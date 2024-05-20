#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

Generate baits for target capture experiments.

Code documentation
------------------
"""


import os
import sys
import csv
import shutil
import argparse
import contextlib
from collections import Counter

import pandas as pd
from Bio import SeqIO
import datapane as dp
import plotly.graph_objs as go

try:
	from __init__ import __version__
	import constants as ct
	import map_utils as mu
	import report_utils as ru
	import cluster_utils as cu
	import general_utils as gu
except ModuleNotFoundError:
	from proBait import __version__
	from proBait import constants as ct
	from proBait import map_utils as mu
	from proBait import report_utils as ru
	from proBait import cluster_utils as cu
	from proBait import general_utils as gu


def incremental_bait_generator(fasta_input, baits_file, output_dir,
							   bait_size, bait_coverage, bait_identity,
							   bait_offset, minimum_region, nr_contigs,
							   short_samples, minimum_exact_match,
							   bait_counts, threads, generate=False, depth=False):
	"""
	"""
	total = 0
	coverage_info = []
	invalid_mappings = []

	short_id = short_samples[fasta_input]
	paf_path = os.path.join(output_dir, short_id+'.paf')

	contigs = gu.import_sequences(fasta_input)
	total_bases = nr_contigs[short_id][2]

	# Run minimap2 to map baits against input
	minimap_std = mu.run_minimap2(fasta_input, baits_file, paf_path, threads)

	# Read PAF file
	# PAF file has variabe number of fields, only the first 12 fields are always present
	paf_lines = gu.read_tabular(paf_path)
	# PAF column headers: Query sequence name, Query sequence length,
	# Query start coordinate (0-based), Query end coordinate (0-based),
	# ‘+’ if query/target on the same strand; ‘-’ if opposite,
	# Target sequence name, Target sequence length,
	# Target start coordinate, Target end coordinate,
	# Number of matching bases,
	# Number bases, including gaps, in the mapping,
	# Mapping quality (0-255 with 255 for missing),
	# NM: Total number of mismatches and gaps in the alignment,
	# ms: DP score of the max scoring segment in the alignment,
	# AS: DP alignment score,
	# nn: Number of ambiguous bases in the alignment,
	# tp: Type of aln: P/primary, S/secondary and I,i/inversion,
	# cm: Number of minimizers on the chain,
	# s1: Chaining score,
	# s2: Chaining score of the best secondary chain,
	# de: Gap-compressed per-base sequence divergence,
	# rl: Length of query regions harboring repetitive seeds,
	# cg: CIGAR string (only in PAF),
	# cs: Difference string.
	# go to https://lh3.github.io/minimap2/minimap2.html#10

	# Filter out alignments with length below bait_coverage*bait_size
	# We use the total number of bases, including gaps
	valid_length = [line
					for line in paf_lines
					if int(line[10]) >= (bait_coverage*bait_size)]

	invalid_mappings.extend([line
							 for line in paf_lines
							 if int(line[10]) < (bait_coverage*bait_size)])

	# Compute alignment identity
	# (number of matching bases divided by total number of bases, including gaps)
	for i in range(len(valid_length)):
		valid_length[i].append(int(valid_length[i][9]) / int(valid_length[i][10]))

	for i in range(len(invalid_mappings)):
		invalid_mappings[i].append(int(invalid_mappings[i][9]) / int(invalid_mappings[i][10]))

	# Filter out alignments below defined identity
	valid_pident = [line
					for line in valid_length
					if line[-1] >= bait_identity]

	invalid_mappings.extend([line
							 for line in valid_length
							 if line[-1] < bait_identity])

	# Match alignment string with regex
	pattern = r':[0-9]+|\*[a-z][a-z]|\+[a-z]+|-[a-z]+'
	for i in range(len(valid_pident)):
		current = valid_pident[i][-2]
		valid_pident[i].append(gu.regex_matcher(current, pattern))

	# Process alignment string for discarded baits
	for i in range(len(invalid_mappings)):
		current = invalid_mappings[i][-2]
		invalid_mappings[i].append(gu.regex_matcher(current, pattern))

	# Filter out alignments that do not have at least
	# N sequential matching bases
	valid_stretch = {i: l[-1]
					 for i, l in enumerate(valid_pident)}
	valid_stretch = {i: [e for e in l if e.startswith(':')]
					 for i, l in valid_stretch.items()}
	valid_stretch = {i: [int(e.lstrip(':')) for e in l]
					 for i, l in valid_stretch.items()}
	valid_stretch = {i: sorted(l, reverse=True)
					 for i, l in valid_stretch.items()}

	invalid_cases = [i
					 for i, l in valid_stretch.items()
					 if l[0] < minimum_exact_match]
	invalid_mappings.extend([valid_pident[i]
							 for i in invalid_cases])

	valid_cases = [i
				   for i, l in valid_stretch.items()
				   if l[0] >= minimum_exact_match]
	valid_pident = [valid_pident[i]
					for i in valid_cases]

	# Get info that matters from discarded baits
	discarded = {}
	for l in invalid_mappings:
		discarded.setdefault(l[5], []).append([int(l[7]), int(l[8]), l[-1]])

	# Save info about discarded baits
	discarded_file = os.path.join(output_dir, '{0}_discarded'.format(short_id))
	gu.pickle_dumper(discarded, discarded_file)

	# Get information about positions that match to determine coverage
	for i in range(len(valid_pident)):
		# Get decomposed cs string
		current = valid_pident[i][-1]
		# Get alignment start position in target
		start = int(valid_pident[i][7])
		# Get bait count to adjust depth of coverage
		if bait_counts is not None:
			current_count = bait_counts.get(valid_pident[i][0], 1)
		else:
			current_count = 1
		# Determine if bait covers all positions in the alignment
		### Check if coverage determination for deletions and insertions is well defined
		valid_pident[i].append(mu.single_position_coverage(current, start, current_count))

	# Identify subsequences that are well covered by baits
	# Covered regions by target sequence
	covered_intervals = {}
	for l in valid_pident:
		# Add target sequence identifier, start and stop in target and dictionary
		# With covered positions in the target
		covered_intervals.setdefault(l[5], []).append([int(l[7]), int(l[8]), l[-1]])

	# Sort covered intervals by start position
	# This groups sequential and overlapping alignments
	covered_intervals_sorted = {k: sorted(v, key=lambda x: x[0])
								for k, v in covered_intervals.items()}

	# Merge overlapping intervals
	# Deepcopy to avoid altering original intervals
	merged_intervals = {k: mu.merge_intervals(v)
						for k, v in covered_intervals_sorted.items()}

	# Determine the number of positions that have depth of coverage of at least 1
	coverage = mu.determine_breadth_coverage(merged_intervals, total_bases)

	# Determine subsequences that are not covered
	missing = [mu.determine_missing_intervals(v, k, len(contigs[k]))
			   for k, v in merged_intervals.items()]

	# Add missing regions for contigs that had 0 baits mapped
	not_mapped = [[{c: [[0, len(contigs[c])]]}, len(contigs[c])]
				  for c in contigs if c not in merged_intervals]
	missing.extend(not_mapped)

	missing_regions = {k: v
					   for i in missing
					   for k, v in i[0].items()}
	# Save missing intervals
	missing_file = os.path.join(output_dir, '{0}_missing'.format(short_id))
	gu.pickle_dumper(missing_regions, missing_file)

	not_covered = sum([i[1] for i in missing])

	coverage_info.extend([*coverage, not_covered])

	# Create baits for missing regions
	if generate is True:
		missing_baits_intervals = {k: mu.cover_intervals(v,
														 len(contigs[k]),
														 bait_size,
														 minimum_region,
														 bait_offset)
								   for k, v in missing_regions.items()}

		# Get sequences of all baits
		baits_seqs = [str(rec.seq)
					  for rec in SeqIO.parse(baits_file, 'fasta')]

		# Create fasta strings
		extra_baits = {}
		for k, v in missing_baits_intervals.items():
			# Check if new baits are not equal to reverse complement of previous baits
			extra_baits[k] = list(set([f'>{short_id}|{k}|{e[0]}\n{contigs[k][e[0]:e[1]]}'
										for e in v
										if contigs[k][e[0]:e[1]] not in baits_seqs
										and gu.reverse_complement(contigs[k][e[0]:e[1]]) not in baits_seqs]))

		new_baits_lines = [v for k, v in extra_baits.items()]
		new_baits_lines = gu.flatten_list(new_baits_lines)

		# Avoid adding newline if there are no baits
		if len(new_baits_lines) > 0:
			gu.write_lines(new_baits_lines, baits_file)

		coverage_info.append(len(new_baits_lines))
		total += len(new_baits_lines)

	if depth is True:
		# Determine depth of coverage
		depth_values = {}
		for k, v in merged_intervals.items():
			depth_values[k] = mu.determine_depth_coverage(v, len(contigs[k]))

		total_counts = {}
		for k, v in depth_values.items():
			for i in v[1]:
				total_counts.setdefault(i[0], []).append(i[1])

		total_counts = {k: sum(v) for k, v in total_counts.items()}

		coverage_info.extend([depth_values, total_counts])

	return [coverage_info, total, discarded_file, missing_file]


def exclude_similar_baits(baits_file, clustering_dir, cluster_identity,
						  cluster_coverage, bait_size, threads):
	"""
	"""
	# Cluster baits and remove based on similarity threshold
	# Create database
	mmseqs_db = os.path.join(clustering_dir, 'mmseqs_db')
	mmseqs_std = cu.create_mmseqs_db(baits_file, mmseqs_db)

	# Output paths
	cluster_db = os.path.join(clustering_dir, 'clusters')
	temp_directory = os.path.join(clustering_dir, 'tmp')
	align_db = os.path.join(clustering_dir, 'alignDB')
	align_out = os.path.join(clustering_dir, 'alignOUT')

	os.mkdir(temp_directory)
	# Clustering
	cluster_std = cu.cluster_baits(mmseqs_db, cluster_db,
								temp_directory, threads)

	# Align clusters
	align_std = cu.align_clusters(mmseqs_db, cluster_db, align_db, threads)

	# Convert alignments
	convert_std = cu.convert_alignmentDB(mmseqs_db, align_db, align_out)

	# Read clustering results
	cluster_lines = gu.read_tabular(align_out)
	clusters = {}
	# Pident at index 2
	for l in cluster_lines:
		clusters.setdefault(l[0], []).append([l[1], l[2], l[3]])

	# Exclude clusters with only the representative
	clusters = {k: v for k, v in clusters.items() if len(v) > 1}

	# Remove representatives from clusters
	clusters = {k: [e for e in v if e[0] != k]
				for k, v in clusters.items()}

	# Get identifiers of baits with identity above threshold
	exclude = [[e for e in v
				if float(e[1]) >= cluster_identity
				and float(e[2]) >= (bait_size*cluster_coverage)]
			   for k, v in clusters.items()]
	exclude = gu.flatten_list(exclude)
	excluded_seqids = [e[0] for e in exclude]

	# Create FASTA without excluded baits
	baits = gu.import_sequences(baits_file)
	baits = {k: v for k, v in baits.items() if k not in excluded_seqids}

	baits_records = ['>{0}\n{1}'.format(k, v) for k, v in baits.items()]
	singular_baits = os.path.join(clustering_dir, 'singular_baits')
	gu.write_lines(baits_records, singular_baits)

	return [singular_baits, excluded_seqids]


def exclude_contaminant(baits_file, exclude, exclude_pident,
						exclude_coverage, bait_size, output_dir, threads):
	"""
	"""
	# Map against target sequences that baits should not be specific for
	gbasename = os.path.basename(exclude).split('.fna')[0]
	paf_path = os.path.join(output_dir, gbasename+'.paf')
	minimap_std = mu.run_minimap2(exclude, baits_file, paf_path, threads)

	# Import mapping results
	mapped_baits = gu.read_tabular(paf_path)
	# Select alignments above defined identity and coverage
	multispecific_baits = [l[0] for l in mapped_baits
							if (int(l[9])/int(l[10])) >= exclude_pident
							and (int(l[3])-int(l[2])) >= (bait_size*exclude_coverage)]

	# Deduplicate baits ientifiers
	multispecific_baits = list(set(multispecific_baits))

	# Remove baits and write final probe set
	baits = gu.import_sequences(baits_file)
	baits = {k: v
			 for k, v in baits.items()
			 if k not in multispecific_baits}

	baits_records = ['>{0}\n{1}'.format(k, v) for k, v in baits.items()]
	decontaminated_baits = os.path.join(output_dir, 'decontaminated_baits.fasta')
	gu.write_lines(baits_records, decontaminated_baits)

	return [decontaminated_baits, multispecific_baits]


def create_report(configs, initial_data, final_data, ref_set,
				  nr_contigs, ordered_contigs,
				  baits_pos, total_baits, bait_counts):
	"""
	"""
	# Create summary text
	summary_text = ru.add_summary_text()

	# Create parameters table
	config_df = pd.DataFrame(configs.items(), columns=['Parameter', 'Value'])
	config_df = config_df.set_index('Parameter')
	parameter_table = dp.Table(config_df)

	# Create coverage table
	coverage_table, coverage_df = ru.coverage_table(initial_data, final_data, ref_set, nr_contigs)

	# Depth of coverage values distribution
	hist_tracers = ru.depth_hists({k: v[4]
								   for k, v in final_data.items()})

	# Missing intervals hist tracers
	missing_intervals_hists_tracers = ru.missing_intervals_hists({k: v[3]
																  for k, v in final_data.items()})

	# Depth of coverage per position
	line_tracers, shapes = ru.depth_lines({k: v[3]
										   for k, v in final_data.items()},
										  ordered_contigs)#, missing_files)

	figs = {}
	for k, v in line_tracers.items():
		fig = go.Figure(data=[*v])
		figs[k] = fig
		figs[k].update_layout(title={'text': 'Depth of coverage per position',
									 'font_size': 20},
							  paper_bgcolor='rgba(255, 255, 255, 0)',
							  plot_bgcolor='#F3F4F6',
							  hovermode='closest',
							  xaxis={'showgrid': False, 'showline': True,
									 'title': {'text': 'Genome position',
											   'font_size': 18,
											   'standoff': 5
											   },
									 'rangeslider': {'visible': True,
													 'range': [1, max(figs[k].data[0].x)]
													 },
									 'range': [1, max(figs[k].data[0].x)]
									 },
							  yaxis={'showgrid': False,
									 'showline': True,
									 'title': {'text': 'Depth of coverage',
											   'font_size': 18,
											   'standoff': 5
											   }
									 },
							  margin={'l': 30, 'r': 30, 't': 30, 'b': 40},
							  template='ggplot2',
							  font_family='sans-serif')

	# Bait start positions per input
	baits_tracers = {k: ru.baits_tracer(baits_pos[k], v)
					 for k, v in ordered_contigs.items() if k in baits_pos}
	for k, v in baits_tracers.items():
		figs[k].add_trace(v)

	# Create shapes for contig boundaries
	for k, v in shapes.items():
		y_value = max(figs[k].data[0].y)
		shapes_tracers, hidden_tracers = ru.create_shapes(list(shapes[k]), y_value)

		# Add shapes for contig boundaries
		figs[k].update_layout(shapes=shapes_tracers, clickmode='event')

		# Add hidden traces to display contig boundaries hover
		for t in hidden_tracers:
			figs[k].add_trace(t)

	# Datapane only accepts Figure objects, cannot pass Bar, Scatter, etc objects
	hist_figs = {}
	for k in hist_tracers:
		hist_figs[k] = go.Figure(data=hist_tracers[k])
		hist_figs[k].update_layout(title={'text': 'Depth of coverage distribution',
										  'font_size': 20
										  },
								   bargap=0.10,
								   paper_bgcolor='rgba(255, 255, 255, 0)',
								   plot_bgcolor='#F3F4F6',
								   hovermode='closest',
								   xaxis={'showgrid': True,
										  'automargin': True,
										  'showline': True,
										  'title': {'text': 'Depth of coverage',
													'font_size': 18,
													'standoff': 5
													}
										  },
								   yaxis={'type': 'log',
										  'showgrid': True,
										  'automargin': True,
										  'showline': True,
										  'title': {'text': 'Count',
													'font_size': 18,
													'standoff': 5
													}
										  },
								   margin={'l': 30, 'r': 30, 't': 30, 'b': 30},
								   template='ggplot2',
								   font_family='sans-serif')

	miss_figs = {}
	for k in missing_intervals_hists_tracers:
		miss_figs[k] = go.Figure(data=missing_intervals_hists_tracers[k])
		miss_figs[k].update_layout(title={'text': 'Uncovered region size distribution',
										  'font_size': 20
										  },
								   bargap=0.10,
								   paper_bgcolor='rgba(255, 255, 255, 0)',
								   plot_bgcolor='#F3F4F6',
								   hovermode='closest',
								   xaxis={'showgrid': True,
										  'automargin': True,
										  'showline': True,
										  'title': {'text': 'Region size',
													'font_size': 18,
													'standoff': 5
													}
										  },
								   yaxis={'type': 'log',
										  'automargin': True,
										  'showgrid': True,
										  'showline': True,
										  'title': {'text': 'Count',
													'font_size': 18,
													'standoff': 5
													}
										  },
								   margin={'l': 30, 'r': 30, 't': 30, 'b': 30},
								   template='ggplot2',
								   font_family='sans-serif')

	# Create big number objects for summary stats page
	total_bn = dp.BigNumber(heading='Total baits',
							value=total_baits)
	mean_coverage_bn = dp.BigNumber(heading='Mean coverage',
									value=round(coverage_df['Breadth of coverage'].mean(), 3))
	mean_depth_bn = dp.BigNumber(heading='Mean depth',
								 value=round(coverage_df['Mean depth of coverage'].mean(), 3))

	# Create bait count histogram
	if bait_counts is not None:
		counted = len(bait_counts)
		counts = Counter(list(bait_counts.values()))
		not_counted = total_baits - counted
		counts.update({1:not_counted})
	else:
		counts = {1: total_baits}

	counts_tracer = ru.bait_count(counts)
	counts_fig = go.Figure(data=counts_tracer)
	counts_fig.update_layout(title={'text': 'Bait count distribution',
									'font_size': 20
									},
							 bargap=0.10,
							 paper_bgcolor='rgba(255, 255, 255, 0)',
							 plot_bgcolor='#F3F4F6',
							 hovermode='closest',
							 xaxis={'showgrid': True,
									'automargin': True,
									'showline': True,
									'title': {'text': 'Count',
											  'font_size': 18,
											  'standoff': 5
											  }
									},
							 yaxis={'type': 'log',
									'showgrid': True,
									'automargin': True,
									'showline': True,
									'title': {'text': 'Number of distinct baits',
											  'font_size': 18,
											  'standoff': 5
											  }
									},
							 margin={'l': 30, 'r': 30, 't': 30, 'b': 30},
							 template='ggplot2',
							 font_family='sans-serif')

	# Create Group objects for Depth page
	depth_groups = []
	for k in figs:
		depth_groups.append(
			dp.Group(
				dp.Group(
					dp.BigNumber(heading='Mean coverage', value=coverage_df[coverage_df['Sample'].str.contains(k)]['Breadth of coverage'].tolist()[0]),
					dp.BigNumber(heading='Mean depth', value=coverage_df[coverage_df['Sample'].str.contains(k)]['Mean depth of coverage'].tolist()[0]),
					dp.BigNumber(heading='Generated baits', value=coverage_df[coverage_df['Sample'].str.contains(k)]['Generated baits'].tolist()[0]),
					columns=3
				),
				dp.Plot(figs[k]),
				dp.Group(
					dp.Plot(hist_figs[k]),
					dp.Plot(miss_figs[k]),
					columns=2
				),
				label=k,
			)
		)

	# Create report object
	report = dp.Report(
		dp.Page(title='Summary', blocks=[
									dp.Group(summary_text, parameter_table, columns=2),
									dp.Group(total_bn, mean_coverage_bn, mean_depth_bn, columns=3),
									counts_fig,
									coverage_table]
				),
		dp.Page(title='Coverage analysis', blocks=[
									dp.Select(blocks=depth_groups, type=dp.SelectType.DROPDOWN)]
				),
		)

	return report


def write_tsv_output(baits_file, output_file):
	"""
	"""
	baits_records = SeqIO.parse(baits_file, 'fasta')
	baits_lines = ['{0}\t{1}'.format(rec.id, str(rec.seq))
				   for rec in baits_records]

	with open(output_file, 'w') as outfile:
		outfile.write('\n'.join(baits_lines)+'\n')

	return len(baits_lines)


# Add option to only determine baits for a set of target loci
# Accept loci annotations and add that info to the coverage per sample pages
def main(input_files, output_directory, generate_baits, baits, bait_proportion,
		 refs, minimum_sequence_length, bait_size, bait_offset,
		 bait_identity, bait_coverage, minimum_region, minimum_exact_match,
		 cluster, cluster_identity, cluster_coverage,
		 exclude, exclude_pident, exclude_coverage, threads,
		 report, report_identities, report_coverages, tsv_output):

	print('proBait version: {0}'.format(__version__))
	print('Authors: {0}'.format(ct.AUTHORS))
	print('Github: {0}'.format(ct.REPOSITORY))
	print('Contacts: {0}'.format(ct.CONTACTS))

	if generate_baits is False and baits is None:
		sys.exit('Please provide the "--generate-baits" argument to '
			  	 'generate baits from input genomes or/and provide a '
			  	 'path to a FASTA file with a set of baits.')
	elif generate_baits is False and baits is not None:
		if report is False:
			sys.exit('User provided a set of baits but does not want to '
				  	 'generate new baits based on the input genomes '
				  	 'nor map the baits to generate a report. Exiting.')

	# Create output directory if it does not exist
	exists = gu.create_directory(output_directory)
	if exists is True:
		sys.exit('Output directory exists. Please provide a path to '
				 'a folder that will be created to save the output files.')

	# Get absolute paths for all input files
	genomes = gu.absolute_paths(input_files)
	# Attribute shorter and unique identifiers to each input
	short_samples = gu.common_suffixes(genomes)
	print(f'Number of inputs: {len(genomes)}')

	initial_baits_file = os.path.join(output_directory, 'initial_baits.fasta')
	if baits is not None:
		shutil.copy(baits, initial_baits_file)
		# Count number of baits provided
		user_baits = len([r.id for r in SeqIO.parse(initial_baits_file, 'fasta')])
	else:
		user_baits = 0

	# Read file with bait proportion
	if bait_proportion is not None:
		bait_counts = gu.read_tabular(bait_proportion)
		bait_counts = {line[0]:int(line[1]) for line in bait_counts}
	else:
		bait_counts = None

	# Determine the number of sequences and total length per input file
	nr_contigs = {short_samples[f]: gu.count_contigs(f, minimum_sequence_length)
				  for f in genomes}

	nr_baits = 0
	if generate_baits is True:
		# Read file with reference basenames
		if refs is not None:
			with open(refs, 'r') as infile:
				ref_set = [g[0] for g in list(csv.reader(infile, delimiter='\t'))]
				ref_set = [os.path.join(input_files, g) for g in ref_set]
		elif refs is None and baits is None:
			ref_set = [genomes[0]]
		else:
			ref_set = []

		# Shred references to generate initial set of baits
		# Does not cover sequence ends
		refs_baits = {}
		if len(ref_set) > 0:
			print('Shredding references to generate initial set of baits...')
			for g in ref_set:
				ref_baits = gu.generate_baits(g, initial_baits_file, bait_size,
											 bait_offset, minimum_sequence_length,
											 short_samples[g])
				refs_baits[g] = ref_baits
				nr_baits += ref_baits
			print(f'Generated {nr_baits} baits based on {len(ref_set)} reference/s.')

		# Identify unique baits
		print('Identifying duplicated baits...')
		baits_file = os.path.join(output_directory, 'baits.fasta')
		total_duplicated, unique_seqids = gu.determine_distinct(initial_baits_file, baits_file)
		print(f'Removed {total_duplicated} duplicated baits.')
		nr_baits -= total_duplicated

		# Map baits against remaining input sequences
		# Maps against references to cover missing regions at contig ends
		bait_creation_dir = os.path.join(output_directory, 'incremental_bait_creation')
		exists = gu.create_directory(bait_creation_dir)

		novel_baits = 0
		coverage_info = {}
		discarded_baits_files = []
		processed = 0
		print('Mapping baits against inputs and generating new baits for uncovered regions...')
		for g in genomes:
			generated = incremental_bait_generator(g, baits_file,
												   bait_creation_dir, bait_size,
												   bait_coverage, bait_identity,
												   bait_offset, minimum_region,
												   nr_contigs, short_samples,
												   minimum_exact_match, None,
												   threads, generate=True, depth=False)
			coverage_info[short_samples[g]] = generated[0]
			if g in ref_set:
				coverage_info[short_samples[g]][3] += refs_baits[g]
			discarded_baits_files.append(generated[2])
			novel_baits += generated[1]
			processed += 1
			print('\r', f'Mapped baits against {processed}/{len(genomes)} inputs.', end='')

		nr_baits += novel_baits

		print(f'\nGenerated {novel_baits} new baits for uncovered regions.')
		print(f'Total of {nr_baits} baits.')
	else:
		baits_file = os.path.join(output_directory, 'baits.fasta')
		shutil.copy(initial_baits_file, baits_file)
		coverage_info = None
		nr_baits = 0
		ref_set = []

	if cluster is True:
		print(f'Clustering baits to remove baits based on identity={cluster_identity} and coverage={cluster_coverage}...')
		clustering_dir = os.path.join(output_directory, 'clustering')
		exists = gu.create_directory(clustering_dir)
		singular_baits, removed = exclude_similar_baits(baits_file,
													    clustering_dir,
													    cluster_identity,
													    cluster_coverage,
													    bait_size,
													    threads)
		print(f'Excluded {len(removed)} baits similar to other baits.')
		# Delete FASTA file with non-decontaminated baits
		os.remove(baits_file)
		# Move and rename FASTA file with decontaminated baits
		shutil.move(singular_baits, baits_file)
		nr_baits -= len(removed)

	exclude_stats = 'None'
	if exclude is not None:
		print(f'Mapping against {exclude} to remove similar baits...')
		exclude_dir = os.path.join(output_directory, 'exclude')
		exists = gu.create_directory(exclude_dir)
		decontaminated_baits, removed = exclude_contaminant(baits_file,
															exclude,
															exclude_pident,
															exclude_coverage,
															bait_size,
															exclude_dir,
															threads)
		print(f'Removed {len(removed)} baits similar to sequences in {exclude}.')
		# Delete FASTA file with non-decontaminated baits
		os.remove(baits_file)
		# Move and rename FASTA file with decontaminated baits
		shutil.move(decontaminated_baits, baits_file)
		# Determine number of exclude regions and total bps
		exclude_stats = [len(rec)
						 for rec in SeqIO.parse(exclude, 'fasta')]
		exclude_stats = len(exclude_stats)
		nr_baits -= len(removed)

	# Create TSV output with bait identifier and bait sequence columns
	if tsv_output is True:
		tsv_output_file = os.path.join(output_directory, 'baits.tsv')
		tsv_baits = write_tsv_output(baits_file, tsv_output_file)
		print(f'Wrote {tsv_baits} baits to {tsv_output_file}')

	if report is True:
		# Create report directory
		report_dir = os.path.join(output_directory, 'report_data')
		exists = gu.create_directory(report_dir)

		# Determine contig order from longest to shortest
		ordered_contigs = gu.order_contigs(genomes, short_samples)

		# Create dict with config values
		configs = {'Number of inputs': len(genomes),
				   'Minimum sequence length': minimum_sequence_length,
				   'Bait size': bait_size,
				   'Bait offset': bait_offset,
				   'Bait identity': bait_identity,
				   'Bait coverage': bait_coverage,
				   'Minimum exact match': minimum_exact_match}

		if refs is not None:
			configs['Number of references'] = len(ref_set)
		else:
			configs['Number of references'] = 'NA'

		configs['Cluster baits'] = str(cluster)
		if cluster is True:
			configs['Cluster identity'] = cluster_identity
			configs['Cluster coverage'] = cluster_coverage

		configs['Sequences to exclude'] = '{0}'.format(exclude_stats)
		if exclude_stats != 'None':
			configs['Exclusion identity'] = exclude_pident
			configs['Exclusion coverage'] = exclude_coverage

		# Get bait positions
		baits_pos = gu.get_baits_pos(baits_file, short_samples)

		# Create a report for each report identity and coverage pair
		print('Evaluating the performance of the set of baits and generating interactive reports...')
		print(f'Identity values: {", ".join(map(str, report_identities))}')
		print(f'Coverage values: {", ".join(map(str, report_coverages))}')
		coverage_dirs = []
		depth_dirs = []
		for i, v in enumerate(report_identities):
			current_identity = v
			current_coverage = report_coverages[i]
			configs['Report bait identity'] = current_identity
			configs['Report bait coverage'] = current_coverage
			# Determine the breadth and depth of coverage for all inputs
			final_coverage_dir = os.path.join(output_directory, 'final_coverage_{0}'.format(i))
			exists = gu.create_directory(final_coverage_dir)
			coverage_dirs.append(final_coverage_dir)

			final_info = {}
			discarded_baits_files = []
			processed = 0
			for g in genomes:
				generated = incremental_bait_generator(g, baits_file,
													   final_coverage_dir,
													   bait_size, current_coverage,
													   current_identity, bait_offset, minimum_region,
													   nr_contigs, short_samples,
													   0, bait_counts,
													   threads, generate=False, depth=True)
				final_info[short_samples[g]] = generated[0]
				discarded_baits_files.append(generated[2])
				processed += 1
				print('\r', f'Processed {processed}/{len(genomes)} inputs for '
					  f'bait_identity={current_identity} and bait_coverage={current_coverage}.', end='')

			# Save depth values
			depth_files_dir = os.path.join(output_directory, 'depth_files_idnt{0}_cov{1}'.format(str(current_identity).replace('.', ''),
																						   str(current_coverage).replace('.', '')))
			exists = gu.create_directory(depth_files_dir)
			depth_dirs.append(depth_files_dir)
			depth_files = [mu.write_depth(k, v[3], depth_files_dir) for k, v in final_info.items()]

			test_report = create_report(configs, coverage_info, final_info, ref_set, nr_contigs,
										ordered_contigs, baits_pos, nr_baits+user_baits,
										bait_counts)

			report_html = os.path.join(output_directory,
									   'proBait_report_idnt{0}_cov{1}.html'.format(configs['Report bait identity'],
																				   configs['Report bait coverage']))
			# Save report as standalone HTML file
			with contextlib.redirect_stdout(None):
				dp.save_report(test_report, path=report_html)

			print('\nCoverage report for bait_identity={0} and '
				  'bait_coverage={1} available in {2}'.format(current_identity, current_coverage, report_dir))

	# Delete intermediate files
	# Initial bait set
	os.remove(initial_baits_file)
	if generate_baits:
		gu.delete_directories([bait_creation_dir])
	if report:
		gu.delete_directories([report_dir]+coverage_dirs+depth_dirs)
	if exclude:
		gu.delete_directories([exclude_dir])
	if cluster:
		gu.delete_directories([clustering_dir])

	print(f'Created a total of {nr_baits} baits to cover {len(genomes)} inputs.')
	print(f'Results available in {output_directory}')


def parse_arguments():

	parser = argparse.ArgumentParser(description=__doc__,
									 formatter_class=argparse.RawDescriptionHelpFormatter)

	parser.add_argument('-i', '--input-files', type=str,
						required=True, dest='input_files',
						help='Path to the directory that contains the '
							 'input FASTA files.')

	parser.add_argument('-o', '--output-directory', type=str,
						required=True, dest='output_directory',
						help='Path to the output directory that will be '
							 'created to store output files. The process '
							 'will exit if the directory already exists.')

	parser.add_argument('-gb', '--generate-baits', action='store_true',
						required=False, dest='generate_baits',
						help='Pass this parameter to generate baits '
							 'based on the sequences in the input files.')

	parser.add_argument('-b', '--baits', type=str,
						required=False, dest='baits',
						help='Path to a FASTA file with a set of baits '
							 'generated on a previous run.')

	parser.add_argument('-bp', '--bait-proportion', type=str,
						required=False, dest='bait_proportion',
						help='Path to a TSV file with data about bait '
							 'proportion to use when evaluating bait '
							 'performance and generating the reports.')

	parser.add_argument('-rf', '--refs', type=str,
						required=False, dest='refs',
						help='Path to a file with the basenames of the input '
							 'FASTA files, one basename per line, that will be '
							 'used as references to create the initial set of '
							 'baits. The references are shredded into baits '
							 'according to the --bait-size and --bait-offset '
							 'values.')

	parser.add_argument('-msl', '--minimum-sequence-length', type=int,
						required=False, default=0,
						dest='minimum_sequence_length',
						help='Do not generate baits for sequences shorter '
							 'than this value.')

	parser.add_argument('-bs', '--bait-size', type=int,
						required=False, default=120,
						dest='bait_size',
						help='The length of the baits in bases. '
							 'All the baits that are generated during the '
							 'process have a size equal to the passed value.')

	parser.add_argument('-bo', '--bait-offset', type=int,
						required=False, default=120,
						dest='bait_offset',
						help='Start position offset between consecutive '
							 'baits.')

	parser.add_argument('-bi', '--bait-identity', type=float,
						required=False, default=0.95,
						dest='bait_identity',
						help='Minimum percent identity to accept an alignment '
							 'between a bait and a region of an input sequence.')

	parser.add_argument('-bc', '--bait-coverage', type=float,
						required=False, default=0.95,
						dest='bait_coverage',
						help='Minimum proportion of a bait that must align to '
							 'accept an alignment.')

	parser.add_argument('-mr', '--minimum-region', type=int,
						required=False, default=0,
						dest='minimum_region',
						help='The process will only generate new baits for '
							 'uncovered regions with length greater than this value.')

	parser.add_argument('-me', '--minimum-exact-match', type=int,
						required=False, default=0,
						dest='minimum_exact_match',
						help='Minimum number of N sequential matching '
							 'bases in an alignment to accept it.')

	parser.add_argument('-c', '--cluster', action='store_true',
						required=False, dest='cluster',
						help='Cluster set of baits to remove similar '
							 'baits and reduce redundancy.')

	parser.add_argument('-ci', '--cluster-identity', type=float,
						required=False, default=1.0,
						dest='cluster_identity',
						help='Exclude baits with an identity value to '
							 'the cluster representative equal or higher '
							 'than this value.')

	parser.add_argument('-cc', '--cluster-coverage', type=float,
						required=False, default=1.0,
						dest='cluster_coverage',
						help='Exclude baits with a coverage value to '
							 'the cluster representative equal or higher '
							 'than this value.')

	parser.add_argument('-e', '--exclude', type=str,
						required=False, default=None,
						dest='exclude',
						help='Path to a FASTA file containing sequences '
							 'to which baits must not be specific.')

	parser.add_argument('-ep', '--exclude-pident', type=float,
						required=False, default=0.7,
						dest='exclude_pident',
						help='Exclude baits with an identity value to '
							 'a region of a sequence to exclude equal '
							 'or higher than this value.')

	parser.add_argument('-ec', '--exclude-coverage', type=float,
						required=False, default=0.5,
						dest='exclude_coverage',
						help='Exclude baits with a coverage value to '
							 'a region of a sequence to exclude equal '
							 'or higher than this value.')

	parser.add_argument('-t', '--threads', type=int,
						required=False, default=1,
						dest='threads',
						help='Number of threads passed to minimap2 and '
							 'MMseqs2.')

	parser.add_argument('-r', '--report', action='store_true',
						required=False, dest='report',
						help='Evaluate bait performance against input '
							 'sequences and generate an interactive report '
							 'with results.')

	parser.add_argument('-ri', '--report-identities', type=float, nargs='+',
						required=False, dest='report_identities',
						help='List of identity values used to evaluate bait '
							 'performance. proBait will generate a report per '
							 'identity value. An equal number of coverage '
							 'values must be provided to the --report-coverages '
							 'parameter to pair with the identity values.')

	parser.add_argument('-rc', '--report-coverages', type=float, nargs='+',
						required=False,
						dest='report_coverages',
						help='List of coverage values used to evaluate bait '
							 'performance. proBait will generate a report per '
							 'coverage value. An equal number of identity '
							 'values must be provided to the --report-identities '
							 'parameter to pair with the coverage values.')

	parser.add_argument('-tsv', '--tsv-output', action='store_true',
						required=False, dest='tsv_output',
						help='Output bait set in TSV format (first column includes '
							 'the bait sequence identifier and the second column '
							 'includes the bait DNA sequence).')

	args = parser.parse_args()

	# Check minimap2 and MMseqs2
	minimap2_path = shutil.which('minimap2')
	mmseqs2_path = shutil.which('mmseqs')
	if minimap2_path is None or mmseqs2_path is None:
		sys.exit('Could not find minimap2 or MMseqs2. Please make sure both are installed.')

	main(**vars(args))


if __name__ == '__main__':

	parse_arguments()
