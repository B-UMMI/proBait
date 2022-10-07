#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""


import os
import sys
import csv
import shutil
import argparse

import pandas as pd
from Bio import SeqIO
import datapane as dp
import plotly.graph_objs as go

#try:
import map_utils as mu
import report_utils as ru
import cluster_utils as cu
import general_utils as gu
#except:
    # import proBait.map_utils as mu
    # import proBait.report_utils as ru
    # import proBait.cluster_utils as cu
    # import proBait.general_utils as gu

# dp.enable_logging()

def incremental_bait_generator(fasta_input, unique_baits, output_dir,
                               bait_size, bait_coverage, bait_identity,
                               bait_offset, minimum_region, nr_contigs,
                               short_samples, minimum_exact_match,
                               generate=False, depth=False):
    """
    """

    total = 0
    coverage_info = []
    invalid_mappings = []

    short_id = short_samples[fasta_input]
    paf_path = os.path.join(output_dir, short_id+'.paf')

    contigs = gu.import_sequences(fasta_input)
    total_bases = nr_contigs[short_id][2]

    # run minimap2 to map baits against input
    minimap_std = mu.run_minimap2(fasta_input, unique_baits, paf_path)

    # read PAF file
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

    # filter out alignments with length below bait_coverage*bait_size
    # we use the total number of bases, including gaps
    valid_length = [line
                    for line in paf_lines
                    if int(line[10]) >= (bait_coverage*bait_size)]

    invalid_mappings.extend([line
                             for line in paf_lines
                             if int(line[10]) < (bait_coverage*bait_size)])

    # compute alignment identity
    # (number of matching bases divided by total number of bases, including gaps)
    for i in range(len(valid_length)):
        valid_length[i].append(int(valid_length[i][9]) / int(valid_length[i][10]))
    
    for i in range(len(invalid_mappings)):
        invalid_mappings[i].append(int(invalid_mappings[i][9]) / int(invalid_mappings[i][10]))

    # filter out alignments below defined identity
    valid_pident = [line
                    for line in valid_length
                    if line[-1] >= bait_identity]

    invalid_mappings.extend([line
                             for line in valid_length
                             if line[-1] < bait_identity])

    # match alignment string with regex
    pattern = r':[0-9]+|\*[a-z][a-z]|\+[a-z]+|-[a-z]+'
    for i in range(len(valid_pident)):
        current = valid_pident[i][-2]
        valid_pident[i].append(gu.regex_matcher(current, pattern))

    # process alignment string for discarded baits
    for i in range(len(invalid_mappings)):
        current = invalid_mappings[i][-2]
        invalid_mappings[i].append(gu.regex_matcher(current, pattern))

    # filter out alignments that do not have at least
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

    # get info that matters from discarded baits
    discarded = {}
    for l in invalid_mappings:
        discarded.setdefault(l[5], []).append([int(l[7]), int(l[8]), l[-1]])

    # save info about discarded baits
    discarded_file = os.path.join(output_dir, '{0}_discarded'.format(short_id))
    gu.pickle_dumper(discarded, discarded_file)

    # get information about positions that match to determine coverage
    for i in range(len(valid_pident)):
        # get decomposed cs string
        current = valid_pident[i][-1]
        # get alignment start position in target
        start = int(valid_pident[i][7])
        # determine if bait covers all positions in the alignment
        ### Check if coverage determination for deletions and insertions is well defined
        valid_pident[i].append(mu.single_position_coverage(current, start))

    # identify subsequences that are well covered by baits
    # covered regions by target sequence
    covered_intervals = {}
    for l in valid_pident:
        # add target sequence identifier, start and stop in target and dictionary
        # with covered positions in the target
        covered_intervals.setdefault(l[5], []).append([int(l[7]), int(l[8]), l[-1]])

    # sort covered intervals by start position
    # this groups sequential and overlapping alignments
    covered_intervals_sorted = {k: sorted(v, key=lambda x: x[0])
                                for k, v in covered_intervals.items()}

    # merge overlapping intervals
    # deepcopy to avoid altering original intervals
    merged_intervals = {k: mu.merge_intervals(v)
                        for k, v in covered_intervals_sorted.items()}

    # determine the number of positions that have depth of coverage of at least 1
    coverage = mu.determine_breadth_coverage(merged_intervals, total_bases)

    # determine subsequences that are not covered
    missing = [mu.determine_missing_intervals(v, k, len(contigs[k]))
               for k, v in merged_intervals.items()]

    # add missing regions for contigs that had 0 baits mapped
    not_mapped = [[{c: [[0, len(contigs[c])]]}, len(contigs[c])]
                  for c in contigs if c not in merged_intervals]
    missing.extend(not_mapped)

    missing_regions = {k: v
                       for i in missing
                       for k, v in i[0].items()}
    # save missing intervals
    missing_file = os.path.join(output_dir, '{0}_missing'.format(short_id))
    gu.pickle_dumper(missing_regions, missing_file)

    not_covered = sum([i[1] for i in missing])

    coverage_info.extend([*coverage, not_covered])

    # create baits for missing regions
    if generate is True:
        missing_baits_intervals = {k: mu.cover_intervals(v, len(contigs[k]), bait_size, minimum_region, bait_offset)
                                   for k, v in missing_regions.items()}

        # get sequences of all probes
        baits_seqs = [str(rec.seq)
                      for rec in SeqIO.parse(unique_baits, 'fasta')]

        # create fasta strings
        extra_probes = {}
        for k, v in missing_baits_intervals.items():
            # check if new baits are not equal to reverse complement of previous baits
            extra_probes[k] = list(set(['>{0}_{1}\n{2}'.format(k, e[0], contigs[k][e[0]:e[1]])
                                        for e in v
                                        if contigs[k][e[0]:e[1]] not in baits_seqs
                                        and gu.reverse_complement(contigs[k][e[0]:e[1]]) not in baits_seqs]))

        new_baits_lines = [v for k, v in extra_probes.items()]
        new_baits_lines = gu.flatten_list(new_baits_lines)

        # avoid adding newline if there are no baits
        if len(new_baits_lines) > 0:
            gu.write_lines(new_baits_lines, unique_baits)

        coverage_info.append(len(new_baits_lines))
        total += len(new_baits_lines)

    if depth is True:
        # determine depth of coverage
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


def exclude_similar_probes(unique_baits, clustering_dir, cluster_identity,
                           cluster_coverage, bait_size, threads):
    """
    """

    print('Clustering probes...')
    # cluster baits and remove based on similarity threshold
    # create database
    mmseqs_db = os.path.join(clustering_dir, 'mmseqs_db')
    mmseqs_std = cu.create_mmseqs_db(unique_baits, mmseqs_db)

    # output paths
    cluster_db = os.path.join(clustering_dir, 'clusters')
    temp_directory = os.path.join(clustering_dir, 'tmp')
    align_db = os.path.join(clustering_dir, 'alignDB')
    align_out = os.path.join(clustering_dir, 'alignOUT')

    os.mkdir(temp_directory)
    # clustering
    cluster_std = cu.cluster_baits(mmseqs_db, cluster_db,
                                temp_directory, threads)

    # align clusters
    align_std = cu.align_clusters(mmseqs_db, cluster_db, align_db, threads)

    # convert alignments
    convert_std = cu.convert_alignmentDB(mmseqs_db, align_db, align_out)

    # read clustering results
    cluster_lines = gu.read_tabular(align_out)
    clusters = {}
    # pident at index 2
    for l in cluster_lines:
        clusters.setdefault(l[0], []).append([l[1], l[2], l[3]])

    # exclude clusters with only the representative
    clusters = {k: v for k, v in clusters.items() if len(v) > 1}

    # remove representatives from clusters
    clusters = {k: [e for e in v if e[0] != k]
                for k, v in clusters.items()}

    # get identifiers of baits with identity above threshold
    exclude = [[e for e in v
                if float(e[1]) >= cluster_identity
                and float(e[2]) >= (bait_size*cluster_coverage)]
               for k, v in clusters.items()]
    exclude = gu.flatten_list(exclude)
    excluded_seqids = [e[0] for e in exclude]
    print('Excluded {0} probes highly similar to other '
          'probes.'.format(len(excluded_seqids)))

    # create FASTA without excluded baits
    baits = gu.import_sequences(unique_baits)
    baits = {k: v for k, v in baits.items() if k not in excluded_seqids}

    baits_records = ['>{0}\n{1}'.format(k, v) for k, v in baits.items()]
    filtered_baits = os.path.join(clustering_dir, 'filtered_baits')
    gu.write_lines(baits_records, filtered_baits)

    return [filtered_baits, excluded_seqids]


def exclude_contaminant(unique_baits, exclude_regions, exclude_pident,
                        exclude_coverage, bait_size, output_dir):
    """
    """

    print('Mapping against and removing similar probes...')
    # map against target genome that baits should not be specific for
    gbasename = os.path.basename(exclude_regions).split('.fna')[0]
    paf_path = os.path.join(output_dir, gbasename+'.paf')
    minimap_std = mu.run_minimap2(exclude_regions, unique_baits, paf_path)

    # import mapping results
    mapped_probes = gu.read_tabular(paf_path)
    # select alignments above defined identity and coverage
    multispecific_probes = [l[0] for l in mapped_probes
                            if (int(l[9])/int(l[10])) >= exclude_pident
                            and (int(l[3])-int(l[2])) >= (bait_size*exclude_coverage)]

    # deduplicate baits ientifiers
    multispecific_probes = list(set(multispecific_probes))

    # remove probes and write final probe set
    baits = gu.import_sequences(unique_baits)
    baits = {k: v
             for k, v in baits.items()
             if k not in multispecific_probes}

    print('Removed {0} probes similar with contaminant '
          'genome.'.format(len(multispecific_probes)))

    baits_records = ['>{0}\n{1}'.format(k, v) for k, v in baits.items()]
    final_baits = os.path.join(output_dir, 'final_baits.fasta')
    gu.write_lines(baits_records, final_baits)

    return [final_baits, multispecific_probes]


# configs = configs
# initial_data = coverage_info
# final_data = final_info
# ref_set = ref_set
# nr_contigs = nr_contigs
# ordered_contigs = ordered_contigs
# missing_files = missing_files
# baits_pos= baits_pos
# total_baits = nr_baits+total
def create_report(configs, initial_data, final_data, ref_set,
                  nr_contigs, ordered_contigs,
                  baits_pos, total_baits):
    """
    """

    # create summary text
    summary_text = ru.add_summary_text()

    # create parameters table
    config_df = pd.DataFrame(configs.items(), columns=['Parameter', 'Value'])
    config_df = config_df.set_index('Parameter')
    parameter_table = dp.Table(config_df)

    # create coverage table
    coverage_table, coverage_df = ru.coverage_table(initial_data, final_data, ref_set, nr_contigs)

    # depth of coverage values distribution
    hist_tracers = ru.depth_hists({k: v[4]
                                   for k, v in final_data.items()})

    # missing intervals hist tracers
    ### Use the same data to create red shapes for missing intervals
    missing_intervals_hists_tracers = ru.missing_intervals_hists({k: v[3]
                                                                  for k, v in final_data.items()})

    # depth of coverage per position
    line_tracers, shapes = ru.depth_lines({k: v[3]
                                           for k, v in final_data.items()},
                                          ordered_contigs)#, missing_files)

    figs = {}
    for k, v in line_tracers.items():
        fig = go.Figure(data=[*v])
        figs[k] = fig
        figs[k].update_layout(title=dict(text='Depth of coverage per position', font_size=20),
                              paper_bgcolor='rgba(255, 255, 255, 0)',
                              plot_bgcolor='#F3F4F6',
                              hovermode='closest',
                              xaxis=dict(showgrid=False, showline=True,
                                         title=dict(text='Genome position', font_size=18, standoff=5),
                                         rangeslider=dict(visible=True, range=[1, max(figs[k].data[0].x)]),
                                         # control plot are that is shown, only show area with data points
                                         range=[1, max(figs[k].data[0].x)]),
                              yaxis=dict(showgrid=False, showline=True,
                                         title=dict(text='Depth of coverage', font_size=18, standoff=5)),
                              margin=dict(l=30, r=30, t=30, b=40),
                              template='ggplot2',
                              font_family='sans-serif')

    # baits start position per input
    baits_tracers = {k: ru.baits_tracer(baits_pos[k], v)
                     for k, v in ordered_contigs.items() if k in baits_pos}
    for k, v in baits_tracers.items():
        figs[k].add_trace(v)

    # create shapes for contig boundaries
    for k, v in shapes.items():
        y_value = max(figs[k].data[0].y)
        shapes_tracers, hidden_tracers = ru.create_shapes(list(shapes[k]), y_value)

        # add shapes for contig boundaries
        figs[k].update_layout(shapes=shapes_tracers, clickmode='event')

        # add hidden tracers to display contig boundaries hover
        for t in hidden_tracers:
            figs[k].add_trace(t)

    # Datapane only accepts Figure objects, cannot pass Bar, Scatter, etc objects
    hist_figs = {}
    for k in hist_tracers:
        hist_figs[k] = go.Figure(data=hist_tracers[k])
        hist_figs[k].update_layout(title=dict(text='Depth of coverage distribution', font_size=20),
                                   bargap=0.10,
                                   paper_bgcolor='rgba(255, 255, 255, 0)',
                                   plot_bgcolor='#F3F4F6',
                                   hovermode='closest',
                                   xaxis=dict(showgrid=True, automargin=True,
                                              showline=True, title=dict(text='Depth of coverage', font_size=18, standoff=5)),
                                   yaxis=dict(type='log', showgrid=True, automargin=True,
                                              showline=True, title=dict(text='Count', font_size=18, standoff=5)),
                                   margin=dict(l=30, r=30, t=30, b=30),
                                   template='ggplot2',
                                   font_family='sans-serif')

    miss_figs = {}
    for k in missing_intervals_hists_tracers:
        miss_figs[k] = go.Figure(data=missing_intervals_hists_tracers[k])
        miss_figs[k].update_layout(title=dict(text='Uncovered region size distribution', font_size=20),
                                   bargap=0.10,
                                   paper_bgcolor='rgba(255, 255, 255, 0)',
                                   plot_bgcolor='#F3F4F6',
                                   hovermode='closest',
                                   xaxis=dict(showgrid=True, automargin=True,
                                              showline=True, title=dict(text='Region size', font_size=18, standoff=5)),
                                   yaxis=dict(type='log', automargin=True,
                                              showgrid=True, showline=True,
                                              title=dict(text='Count', font_size=18, standoff=5)),
                                   margin=dict(l=30, r=30, t=30, b=30),
                                   template='ggplot2',
                                   font_family='sans-serif')

    # create big number objects for summary stats page
    total_bn = dp.BigNumber(heading='Total baits', value=total_baits)
    mean_coverage_bn = dp.BigNumber(heading='Mean coverage', value=round(coverage_df['Breadth of coverage'].mean(), 3))
    mean_depth_bn = dp.BigNumber(heading='Mean depth', value=round(coverage_df['Mean depth of coverage'].mean(), 3))

    # create Group objects for Depth page
    depth_groups = []
    for k in figs:
        depth_groups.append(
            dp.Group(
                dp.Group(
                    dp.BigNumber(heading='Mean coverage', value=coverage_df[coverage_df['Sample'].str.contains(k)]['Breadth of coverage'].tolist()[0]),
                    dp.BigNumber(heading='Mean depth', value=coverage_df[coverage_df['Sample'].str.contains(k)]['Mean depth of coverage'].tolist()[0]),
                    dp.BigNumber(heading='Generated probes', value=coverage_df[coverage_df['Sample'].str.contains(k)]['Generated probes'].tolist()[0]),
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

    # create report object
    report = dp.Report(
        dp.Page(title='Summary', blocks=[
                                    dp.Group(summary_text, parameter_table, columns=2),
                                    dp.Group(total_bn, mean_coverage_bn, mean_depth_bn, columns=3),
                                    coverage_table]
                ),
        dp.Page(title='Coverage analysis', blocks=[
                                    dp.Select(blocks=depth_groups, type=dp.SelectType.DROPDOWN)]
                ),
        layout=dp.PageLayout.SIDE
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


input_files = 'assemblies'
output_directory = 'bait_generation'
generate_baits = False
#baits = None
baits = '/home/rmamede/Desktop/rmamede/Cloned_Repos/proBait/proBait/unique_baits.fasta'
bait_size = 120
bait_offset = 60
refs = None
bait_identity = 0.85
bait_coverage = 0.9
minimum_region = 60
minimum_exact_match = 40
cluster_probes = False
cluster_identity = 0.85
cluster_coverage = 0.9
minimum_sequence_length = 120
exclude_regions = None
exclude_pident = 0.5
exclude_coverage = 0.6
threads = 4
report = True
report_identities = [0.95]
report_coverages = [0.95]
tsv_output = True

# number of SNPs, deletions, insertions and etc
# final table with number of uncovered regions per sample, number of SNPs, number of deletions, insertions, etc.

# add option to determine baits from target regions and then only generate baits
# to capture diversity in those regions in other genomes (this allows to determine baits
# only for targeted regions and will not generate baits for other uncovered loci)
# users can provide annotation labels for the target regions that are included in
# the hovertext

# for the step that maps against the human genome, give 1 chromossome at a time

# option to receive baits as input and add more baits generated based on input genomes

# color baits according to feature that led to bait creation (insertion, deletion, etc)

# plot with markers representing baits should also have the following lines:
# - Lines for the intervals of baits that were accepted (green?)
# - Lines for intervals of baits that were not accepted (red?)
# - Lines for regions that had no baits mapped (grey?)
# The lines for baits that were not accepted should also have information
# about why the baits were not accepted (include SNO, deletion, insertion
# info somehow...). Adding these lines should not increase report file size
# that much because each interval will only have 2 points to represent
# accepted baits, refused baits or uncovered regions.
def main(input_files, output_directory, generate_baits, baits, bait_proportion,
         refs, minimum_sequence_length, bait_size, bait_offset,
         bait_identity, bait_coverage, minimum_region, minimum_exact_match,
         cluster_probes, cluster_identity, cluster_coverage,
         exclude_regions, exclude_pident, exclude_coverage, threads,
         report, report_identities, report_coverages, tsv_output):

    if generate_baits is False and baits is None:
        print('Please provide the "--generate-baits" argument to '
              'generate baits from input genomes or/and provide a '
              'path to a FASTA file with a set of baits.')
        sys.exit()
    elif generate_baits is False and baits is not None:
        if report is False:
            print('Provided set of baits but does not want to '
                  'generate new baits based on the input genomes '
                  'nor map the baits to generate a report. Exiting.')
            sys.exit()

    # create output directory if it does not exist
    exists = gu.create_directory(output_directory)

    baits_file = os.path.join(output_directory, 'baits.fasta')
    if baits is not None:
        shutil.copy(baits, baits_file)
        # count number of baits provided
        user_baits = len([r.id for r in SeqIO.parse(baits_file, 'fasta')])
    else:
        user_baits = 0

    # get absolute paths for all input files
    genomes = gu.absolute_paths(input_files)
    # attribute shorter and unique identifiers to each input
    short_samples = gu.common_suffixes(genomes)
    # determine number of sequences and total length per input file
    nr_contigs = {short_samples[f]: gu.count_contigs(f, minimum_sequence_length)
                  for f in genomes}

    if generate_baits is True:
        # read file with references basenames
        if refs is not None:
            with open(refs, 'r') as infile:
                ref_set = [g[0] for g in list(csv.reader(infile, delimiter='\t'))]
                ref_set = [os.path.join(input_files, g) for g in ref_set]
        elif refs is None and baits is None:
            ref_set = [genomes[0]]
        else:
            ref_set = []

        # shred refs to get initial set of baits
        # does not generate some baits to cover sequence ends
        if len(ref_set) > 0:
            for g in ref_set:
                nr_baits = gu.generate_baits(g, baits_file, bait_size,
                                             bait_offset, minimum_sequence_length)

            print('\nGenerated set of {0} probes based on {1} reference '
                  'inputs.'.format(nr_baits, len(ref_set)))

        # identify unique baits
        unique_baits = os.path.join(output_directory, 'unique_baits.fasta')
        total, unique_seqids = gu.determine_distinct(baits_file, unique_baits)
        print('Removed {0} repeated probes.\n'.format(total))
    
        # start mapping baits against remaining genomes
        # maps against reference inputs to cover missing regions at contig ends
        bait_creation_dir = os.path.join(output_directory, 'incremental_bait_creation')
        exists = gu.create_directory(bait_creation_dir)
    
        total = 0
        coverage_info = {}
        discarded_baits_files = []
        #missing_files = []
        processed = 0
        for g in genomes:
            generated = incremental_bait_generator(g, unique_baits,
                                                   bait_creation_dir, bait_size,
                                                   bait_coverage, bait_identity,
                                                   bait_offset, minimum_region,
                                                   nr_contigs, short_samples,
                                                   minimum_exact_match,
                                                   generate=True, depth=False)
            coverage_info[short_samples[g]] = generated[0]
            if g in ref_set:
                coverage_info[short_samples[g]][3] += nr_baits
            discarded_baits_files.append(generated[2])
            #missing_files.append(generated[3])
            total += generated[1]
            processed += 1
            print('\r', 'Generated baits for {0}/{1} '
                  'inputs.'.format(processed, len(genomes)), end='')
    
        print('\nAdded {0} probes to cover {1} assemblies.\nTotal '
              'of {2} probes.'.format(total, len(genomes),
                                      nr_baits+total))
    else:
        unique_baits = baits_file
        coverage_info = None
        nr_baits = 0
        ref_set = []
        total = 0

    if cluster_probes is True:
        clustering_dir = os.path.join(output_directory, 'clustering')
        exists = gu.create_directory(clustering_dir)
        unique_baits, removed = exclude_similar_probes(unique_baits,
                                                       clustering_dir,
                                                       cluster_identity,
                                                       cluster_coverage,
                                                       bait_size,
                                                       threads)
        total -= len(removed)

    # need to copy the Fasta file without the excluded baits to the main output directory!
    exclude_stats = 'None'
    if exclude_regions is not None:
        exclude_dir = os.path.join(output_directory, 'exclude')
        exists = gu.create_directory(exclude_dir)
        unique_baits, removed = exclude_contaminant(unique_baits,
                                                    exclude_regions,
                                                    exclude_pident,
                                                    exclude_coverage,
                                                    bait_size,
                                                    exclude_dir)
        # determine number of exclude regions and total bps
        exclude_stats = [len(rec)
                         for rec in SeqIO.parse(exclude_regions, 'fasta')]
        exclude_stats = len(exclude_stats)
        total -= len(removed)

    # create TSV output with bait identifier and bait sequence columns
    if tsv_output is True:
        tsv_output_file = os.path.join(output_directory, 'unique_baits.tsv')
        tsv_baits = write_tsv_output(unique_baits, tsv_output_file)
        print('Wrote {0} baits to {1}'.format(tsv_baits, tsv_output_file))

    if report is True:
        # create report directory
        report_dir = os.path.join(output_directory, 'report_data')
        exists = gu.create_directory(report_dir)

        # determine contig order from longest to shortest
        ordered_contigs = gu.order_contigs(genomes, short_samples)

        # create dict with config values
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

        configs['Cluster probes'] = str(cluster_probes)
        if cluster_probes is True:
            configs['Cluster identity'] = cluster_identity
            configs['Cluster coverage'] = cluster_coverage
        
        configs['Sequences to exclude'] = '{0}'.format(exclude_stats)
        if exclude_stats != 'None':
            configs['Exclusion identity'] = exclude_pident
            configs['Exclusion coverage'] = exclude_coverage

        # get baits positions
        baits_pos = gu.get_baits_pos(unique_baits, short_samples)

        # create a report for each bait identity and coverage pair
        for i, v in enumerate(report_identities):
            configs['Report bait identity'] = v
            configs['Report bait coverage'] = report_coverages[i]
            # determine breadth of coverage for all assemblies
            # and depth of coverage for each base
            final_coverage_dir = os.path.join(output_directory, 'final_coverage_{0}'.format(i))
            exists = gu.create_directory(final_coverage_dir)

            final_info = {}
            discarded_baits_files = []
            processed = 0
            for g in genomes:
                generated = incremental_bait_generator(g, unique_baits,
                                                       final_coverage_dir,
                                                       bait_size, report_coverages[i],
                                                       v, bait_offset, minimum_region,
                                                       nr_contigs, short_samples,
                                                       0,
                                                       generate=False, depth=True)
                final_info[short_samples[g]] = generated[0]
                discarded_baits_files.append(generated[2])
                processed += 1
                print('\r', 'Evaluated coverage for {0}/{1} '
                      'inputs with bait_identity={2} and '
                      'bait_coverage={3}.'.format(processed, len(genomes),
                                                  v, report_coverages[i]),
                      end='')

            # save depth values
            depth_files_dir = os.path.join(output_directory, 'depth_files_idnt{0}_cov{1}'.format(str(v).replace('.', ''),
                                                                                           str(report_coverages[i]).replace('.', '')))
            exists = gu.create_directory(depth_files_dir)
            depth_files = [mu.write_depth(k, v[3], depth_files_dir) for k, v in final_info.items()]

            test_report = create_report(configs, coverage_info, final_info, ref_set, nr_contigs,
                                        ordered_contigs, baits_pos, nr_baits+total+user_baits)

            report_html = os.path.join(output_directory,
                                       'proBait_report_idnt{0}_cov{1}.html'.format(configs['Report bait identity'],
                                                                                   configs['Report bait coverage']))
            test_report.save(path=report_html,
                             formatting=dp.ReportFormatting(
                                 font=dp.FontChoice.SANS,
                                 width=dp.ReportWidth.FULL))

            print('\nCoverage report for bait_identity={0} and '
                  'bait_coverage={1} available in {2}'.format(v, report_coverages[i], report_dir))


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input-files', type=str,
                        required=True, dest='input_files',
                        help='Path to the directory with '
                             'input FASTA files.')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=True, dest='output_directory',
                        help='Path to the output directory where '
                             'files will be saved to (must not exist, '
                             'process will create this directory).')

    parser.add_argument('-gb', '--generate-baits', required=False,
                        action='store_true', dest='generate_baits',
                        help='')

    parser.add_argument('-b', '--baits', type=str,
                        required=False, dest='baits',
                        help='')

    parser.add_argument('-bp', '--bait-proportion', type=str,
                        required=False, dest='bait_proportion',
                        help='')

    parser.add_argument('-rf', '--refs', type=str,
                        required=False,
                        dest='refs',
                        help='Path to file with the basename of files '
                             'that will be used as references to create '
                             'the initial set of probes.')

    parser.add_argument('-mcl', '--minimum-sequence-length', type=int,
                        required=False, default=0,
                        dest='minimum_sequence_length',
                        help='Minimum sequence length. Probes will '
                             'not be created for sequences with a '
                             'length value that is smaller than '
                             'this value.')

    parser.add_argument('-bs', '--bait-size', type=int,
                        required=False,
                        default=120,
                        dest='bait_size',
                        help='Length of the probes that the process '
                             'will create.')

    parser.add_argument('-bo', '--bait-offset', type=int,
                        required=False,
                        default=120,
                        dest='bait_offset',
                        help='Start position offset between consecutive '
                             'probes.')

    parser.add_argument('-bi', '--bait-identity', type=float,
                        required=False,
                        default=1.0,
                        dest='bait_identity',
                        help='Minimum percent identity that '
                             'aligned probes need to have when '
                             'aligned against input genomes ('
                             'probes with lower identity are '
                             'not included in the set of probes '
                             'that cover well regions of input '
                             'genomes).')

    parser.add_argument('-bc', '--bait-coverage', type=float,
                        required=False,
                        default=1.0,
                        dest='bait_coverage',
                        help='Minimum percent length of the '
                             'probe that has to align against the '
                             'genome assembly (probes with lower '
                             'coverage are not included in the set '
                             'of probes that cover well regions of '
                             'input genomes).')

    parser.add_argument('-mr', '--minimum-region', type=int,
                        required=False,
                        default=0,
                        dest='minimum_region',
                        help='Uncovered regions must have a length '
                             'value equal or greater than this value. '
                             'If the uncovered region is smaller than '
                             'this value the process will not generate '
                             'new baits to cover that region.')

    parser.add_argument('-me', '--minimum-exact-match', type=int,
                        required=False,
                        default=0,
                        dest='minimum_exact_match',
                        help='Minimum number of N sequential matching '
                             'bases. Alignments are only accepted if they '
                             'have at least N sequential matching bases '
                             'with an input.')

    parser.add_argument('-c', '--cluster-probes', required=False,
                        action='store_true',
                        dest='cluster_probes',
                        help='Cluster set of probes after generating '
                             'probes to cover all input assemblies. '
                             'This clustering step will cluster '
                             'similar probes and remove highly similar '
                             'probes based on percent identity and '
                             'coverage.')

    parser.add_argument('-ci', '--cluster-identity', type=float,
                        required=False,
                        default=1.0,
                        dest='cluster_identity',
                        help='Clustered probes with equal or higher '
                             'percent identity are excluded.')

    parser.add_argument('-cc', '--cluster-coverage', type=float,
                        required=False,
                        default=1.0,
                        dest='cluster_coverage',
                        help='Clustered probes with equal or higher '
                             'coverage may be excluded based on '
                             'percent identity.')

    parser.add_argument('-e', '--exclude-regions', type=str,
                        required=False,
                        default=None,
                        dest='exclude_regions',
                        help='Path to a FASTA file with genomic regions '
                             'that probes must not cover.')

    parser.add_argument('-ep', '--exclude-pident', type=float,
                        required=False,
                        default=0.8,
                        dest='exclude_pident',
                        help='Probes with percent identity equal or '
                             'higher than this value to regions that '
                             'must not be covered will be excluded.')

    parser.add_argument('-ec', '--exclude-coverage', type=float,
                        required=False,
                        default=0.5,
                        dest='exclude_coverage',
                        help='Probes that map against the regions to '
                             'exclude with equal or greater coverage '
                             'may be excluded based on percent identity.')

    parser.add_argument('-t', '--threads', type=int,
                        required=False,
                        default=1, dest='threads',
                        help='')

    parser.add_argument('-r', '--report', required=False,
                        action='store_true',
                        dest='report',
                        help='')

    parser.add_argument('-ri', '--report-identities', type=float, nargs='+',
                        required=False,
                        dest='report_identities',
                        help='')

    parser.add_argument('-rc', '--report-coverages', type=float, nargs='+',
                        required=False,
                        dest='report_coverages',
                        help='')

    parser.add_argument('-tsv', '--tsv-output', required=False,
                        action='store_true',
                        dest='tsv_output',
                        help='')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))
