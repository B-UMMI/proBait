#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""


import os
import sys
import argparse

from Bio import SeqIO

try:
    import map_utils as mu
    import report_utils as ru
    import cluster_utils as cu
    import general_utils as gu
except:
    import proBait.map_utils as mu
    import proBait.report_utils as ru
    import proBait.cluster_utils as cu
    import proBait.general_utils as gu


#fasta_input = genomes[0]
#output_dir = bait_creation_dir
def incremental_bait_generator(fasta_input, unique_baits, output_dir, bait_size,
                               bait_coverage, bait_identity, minimum_region,
                               nr_contigs, short_samples,
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

    # filter out secondary alignments???

    # filter out alignments below bait_coverage*bait_size
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
    valid_pident = [line for line in valid_length if line[-1] >= bait_identity]

    invalid_mappings.extend([line for line in valid_length if line[-1] < bait_identity])

    # match alignment string with regex
    pattern = r':[0-9]+|\*[a-z][a-z]|\+[a-z]+|-[a-z]+'
    for i in range(len(valid_pident)):
        current = valid_pident[i][-2]
        valid_pident[i].append(gu.regex_matcher(current, pattern))

    # process alignment string for discarded baits
    for i in range(len(invalid_mappings)):
        current = invalid_mappings[i][-2]
        invalid_mappings[i].append(gu.regex_matcher(current, pattern))

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
    not_mapped = [[{c: [[0, len(contigs[c])]]}, len(contigs[c])] for c in contigs if c not in merged_intervals]
    missing.extend(not_mapped)

    missing_regions = {k: v for i in missing for k, v in i[0].items()}
    # save missing intervals
    missing_file = os.path.join(output_dir, '{0}_missing'.format(short_id))
    gu.pickle_dumper(missing_regions, missing_file)

    not_covered = sum([i[1] for i in missing])

    coverage_info.extend([*coverage, not_covered])

    # create baits for missing regions
    if generate is True:
        missing_baits_intervals = {k: mu.cover_intervals(v, len(contigs[k]), bait_size, minimum_region)
                                   for k, v in missing_regions.items()}

        # get sequences of all probes
        baits_seqs = [str(rec.seq) for rec in SeqIO.parse(unique_baits, 'fasta')]
        #baits_reverse = [gu.reverse_complement(seq) for seq in baits_seqs]
        #baits_seqs.extend(baits_reverse)

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
    multispecific_probes = [l[0] for l in mapped_probes
                            if (int(l[9])/int(l[10])) >= exclude_pident
                            and (int(l[3])-int(l[2])) >= (bait_size*exclude_coverage)]

    # remove probes and write final probe set
    baits = gu.import_sequences(unique_baits)
    baits = {k: v for k, v in baits.items() if k not in multispecific_probes}

    print('Removed {0} probes similar with contaminant '
          'genome.'.format(len(multispecific_probes)))

    baits_records = ['>{0}\n{1}'.format(k, v) for k, v in baits.items()]
    final_baits = os.path.join(output_dir, 'final_baits.fasta')
    gu.write_lines(baits_records, final_baits)

    return [final_baits, multispecific_probes]


############## Baits positions are wrong! (offset of -1 in each position?)
#initial_data = coverage_info
#final_data = final_info
#short_ids = short_samples
#total_baits = nr_baits+total
#output_dir = report_dir
#initial_baits = nr_baits
#iter_baits = total
def create_report(initial_data, final_data, output_dir, short_ids,
                  ordered_contigs, fixed_xaxis, fixed_yaxis, ref_set,
                  nr_contigs, configs, baits_pos, total_baits,
                  initial_baits, iter_baits):
    """
    """

    # check if user wants equal yaxis ranges for all line plots
    max_x = None
    if fixed_xaxis is True:
        max_x = max([v[2] for v in nr_contigs.values()])

    max_y = None
    coverage_values = {k: max(list(v[4].keys())) for k, v in final_data.items()}
    if fixed_yaxis is True:
        max_y = max(coverage_values.values())

    # create table with run summary
    config_table = ru.create_table_tracer(['Parameter', 'Value'],
                                          dict(size=16),
                                          dict(color='#ffffff', width=2),
                                          dict(color='#9ecae1'),
                                          [list(configs.keys()),
                                           list(configs.values())],
                                          dict(size=14),
                                          dict(color='#ffffff', width=1),
                                          dict(color='#f0f0f0'),
                                          dict(x=[0.5, 1.0], y=[0.2, 1.0]))

    table_tracer = ru.coverage_table(initial_data, final_data,
                                     short_ids, ref_set, nr_contigs)

    # depth of coverage values distribution
    hist_tracers = ru.depth_hists({k: v[4] for k, v in final_data.items()})

    # depth of coverage per position
    line_tracers, shapes = ru.depth_lines({k: v[3]
                                           for k, v in final_data.items()},
                                          ordered_contigs)

    # baits start position per input
    baits_tracers = {k: ru.baits_tracer(baits_pos[k], v)
                     for k, v in ordered_contigs.items() if k in baits_pos}

    nr_rows = len(line_tracers) + 4

    titles = ru.subplot_titles(list(short_ids.values()))

    specs_def = ru.report_specs(len(line_tracers))

    # define figure height
    total_height, row_heights = ru.figure_height(190, 150, 505, len(line_tracers))

    # adding vertical_spacing in combination with row_heights can
    # exceed valid values and lead to error
    fig = ru.create_subplots_fig(nr_rows, 2, titles, specs_def,
                                 shared_yaxes=True, row_heights=row_heights)

    # change subplots titles positions
    fig = ru.adjust_subplot_titles(fig)

    # add traces to fig
    # configuration table
    fig.add_trace(config_table, row=1, col=2)

    # coverage stats table
    fig.add_trace(table_tracer, row=3, col=1)

    # depth of coverage plots traces and baits start position traces
    r = 5
    for k, v in line_tracers.items():
        top_x = nr_contigs[k] if max_x is None else max_x
        top_y = coverage_values[k] if max_y is None else max_y
        traces = [v[0], baits_tracers[k], hist_tracers[k]]
        fig = ru.add_plots_traces(traces, r, 1, top_x, top_y, fig)
        r += 1

    # create shapes for contig boundaries
    ref_axis = 1
    shapes_tracers = []
    hidden_tracers = []
    start_hidden = 5
    for k, v in shapes.items():
        y_value = coverage_values[k] if max_y is None else max_y
        traces = ru.create_shapes(list(shapes[k]), y_value, ref_axis)

        shapes_tracers.extend(traces[0])
        hidden_tracers.append([traces[1], start_hidden, 1])

        ref_axis += 2
        start_hidden += 1

    # add shapes for contig boundaries
    fig.update_layout(shapes=shapes_tracers, clickmode='event')

    # add hidden tracers to display contig boundaries hover
    for h in hidden_tracers:
        traces = h[0]
        row = h[1]
        col = h[2]
        for t in traces:
            fig.add_trace(t, row, col)

    # disable grid
    fig.update_xaxes(showgrid=False)
    fig.update_yaxes(showgrid=False, showline=True, linecolor='black')

    # bars with distribution of depth values in logscale
    for i in range(5, 5+len(line_tracers)):
        fig.update_xaxes(type='log', row=i, col=2)

    fig = ru.add_plots_titles(fig)

    # add summary text
    fig = ru.add_summary_text(fig, total_baits, initial_baits,
                              configs['Bait size'], configs['Bait offset'],
                              iter_baits, total_height)

    # line plots need fixed space
    fig.update_layout(title=dict(text='<b>proBait - Coverage Report</b>',
                                 x=0,
                                 xref='paper',
                                 xanchor='left',
                                 font=dict(size=30)),
                      height=total_height,
                      width=1920,
                      template='ggplot2',
                      plot_bgcolor='rgba(0,0,0,0)')

    output_plot = os.path.join(output_dir, 'proBait_report_idnt{0}_cov{1}.html'
                               ''.format(configs['Report bait identity'],
                                         configs['Report bait coverage']))
    ru.create_html_report(fig, output_plot)

    return fig


#input_files = '/home/rfm/Desktop/rfm/Lab_Analyses/pneumo_baits_design/assemblies'
#input_files = '/home/rfm/Desktop/rfm/Lab_Analyses/pneumo_baits_design/test_assemblies'
#input_files = '/home/rfm/Desktop/rfm/Lab_Analyses/pneumo_baits_design/single_assembly'
input_files = '/home/rfm/Desktop/rfm/Lab_Analyses/pneumo_baits_design/double_assemblies'
output_directory = '/home/rfm/Desktop/rfm/Lab_Analyses/pneumo_baits_design/tmp'
bait_size = 120
bait_offset = 120
number_of_refs = 1
bait_identity = 1.0
bait_coverage = 1.0
minimum_region = 120
cluster_probes = False
cluster_identity = 0.85
cluster_coverage = 0.9
minimum_sequence_length = 120
#exclude_regions = None
#exclude_regions = '/home/rfm/Desktop/rfm/Lab_Analyses/pneumo_baits_design/ncbi-genomes-2020-11-16/GCF_000001405.39_GRCh38.p13_genomic.fna'
exclude_regions = None
exclude_pident = 0.7
exclude_coverage = 0.7
threads = 4
fixed_xaxis = True
fixed_yaxis = True
report = True
report_identities = [1.0]
report_coverages = [1.0]

# determine length of subsequences with 0 coverage and etc
# number of SNPs, deletions, insertions and etc
# final table with number of uncovered regions per sample, number of SNPs, number of deletions, insertions, etc.

# add option to determine baits from target regions and then only generate baits
# to capture diversity in those regions in other genomes (this allows to determine baits
# only for targeted regions and will not generate baits for other uncovered loci)
# users can provide annotation labels for the target regions that are included in
# the hovertext

# for the step that maps against the human genome, give 1 chromossome at a time

# option to receive baits as input and add more baits generated based on input genomes

# add boxplot with distribution of identity values for mapped baits against each genome

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
def main(input_files, output_directory, mode, minimum_sequence_length,
         number_of_refs, bait_size, bait_offset,
         bait_identity, bait_coverage, minimum_region,
         cluster_probes, cluster_identity, cluster_coverage,
         exclude_regions, exclude_pident, exclude_coverage, threads,
         report, report_identities, report_coverages, fixed_xaxis,
         fixed_yaxis):

    # create output directory if it does not exist
    exists = gu.create_directory(output_directory)

    # get absolute paths for all input files
    genomes = gu.absolute_paths(input_files)

    # attribute shorter and unique identifiers to each input
    short_samples = gu.common_suffixes(genomes)
    inv_short = {v: k for k, v in short_samples.items()}

    # determine number of sequences and total length per input file
    nr_contigs = {short_samples[f]: gu.count_contigs(f, minimum_sequence_length)
                  for f in genomes}

    # sort inputs based on number of valid sequences and select assemblies
    # with lowest number of contigs as refs
    sorted_contigs = sorted(list(nr_contigs.items()), key=lambda x: x[1][1])
    ref_set = [t[0] for t in sorted_contigs[0:number_of_refs]]

    # get absolute path for sorted inputs
    genomes = [inv_short[g[0]] for g in sorted_contigs]

    # shred refs to get initial set of baits
    # not generating kmers that cover the end of the sequences!
    baits_file = os.path.join(output_directory, 'baits.fasta')
    for g in ref_set:
        nr_baits = gu.generate_baits(inv_short[g], baits_file, bait_size,
                                     bait_offset, minimum_sequence_length)

    print('\nCreated initial set of {0} probes based on {1} '
          'inputs.'.format(nr_baits, number_of_refs))

    # identify unique baits
    unique_baits = os.path.join(output_directory, 'unique_baits.fasta')
    total, unique_seqids = gu.determine_distinct(baits_file, unique_baits)
    print('Removed {0} repeated probes.\n'.format(total))

    # start mapping baits against remaining genomes
    # mapping against ref_set to cover missing regions
    bait_creation_dir = os.path.join(output_directory, 'incremental_bait_creation')
    exists = gu.create_directory(bait_creation_dir)

    total = 0
    coverage_info = {}
    discarded_baits_files = []
    missing_files = []
    processed = 0
    for g in genomes:
        generated = incremental_bait_generator(g, unique_baits,
                                               bait_creation_dir, bait_size,
                                               bait_coverage, bait_identity,
                                               minimum_region,
                                               nr_contigs, short_samples,
                                               generate=True, depth=False)
        coverage_info[short_samples[g]] = generated[0]
        discarded_baits_files.append(generated[2])
        missing_files.append(generated[3])
        total += generated[1]
        processed += 1
        print('\r', 'Generated baits for {0}/{1} '
              'inputs.'.format(processed, len(genomes)), end='')

    print('\nAdded {0} probes to cover {1} assemblies.\nTotal '
          'of {2} probes.'.format(total, len(genomes),
                                  nr_baits+total))

    if cluster_probes is True:
        clustering_dir = os.path.join(output_directory, 'clustering')
        exists = gu.create_directory(clustering_dir)
        unique_baits, removed = exclude_similar_probes(unique_baits,
                                                       clustering_dir,
                                                       cluster_identity,
                                                       cluster_coverage,
                                                       bait_size,
                                                       threads)

    exclude_stats = 'None'
    if exclude_regions is not None:
        exclude_dir = os.path.join(output_directory, 'exclude')
        os.mkdir(exclude_dir)
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

    if report is True:
        # create report directory
        report_dir = os.path.join(output_directory, 'report_data')
        os.mkdir(report_dir)

        # determine contig order from longest to shortest
        ordered_contigs = gu.order_contigs(genomes, short_samples)

        # create dict with config values
        configs = {'Number of inputs': len(genomes),
                   'Minimum sequence length': minimum_sequence_length,
                   'Number of references': number_of_refs,
                   'Bait size': bait_size,
                   'Bait offset': bait_offset,
                   'Bait identity': bait_identity,
                   'Bait coverage': bait_coverage,
                   'Cluster probes': str(cluster_probes),
                   'Cluster identity': cluster_identity,
                   'Cluster coverage': cluster_coverage,
                   'Exclude regions': '{0} regions'.format(exclude_stats),
                   'Exclude identity': exclude_pident,
                   'Exclude coverage': exclude_coverage,
                   'Report bait identity': None,
                   'Report bait coverage': None}

        # get baits positions
        baits_pos = gu.get_baits_pos(unique_baits, short_samples)

        # create a report for each bait identity and coverage pair
        for i, v in enumerate(report_identities):
            configs['Report bait identity'] = v
            configs['Report bait coverage'] = report_coverages[i]
            # determine breadth of coverage for all assemblies
            # and depth of coverage for each base
            final_coverage_dir = os.path.join(output_directory, 'final_coverage_{0}'.format(i))
            os.mkdir(final_coverage_dir)

            final_info = {}
            discarded_baits_files = []
            missing_files = []
            processed = 0
            for g in genomes:
                generated = incremental_bait_generator(g, unique_baits,
                                                       final_coverage_dir,
                                                       bait_size, report_coverages[i],
                                                       v, minimum_region,
                                                       nr_contigs, short_samples,
                                                       generate=False, depth=True)
                final_info[short_samples[g]] = generated[0]
                discarded_baits_files.append(generated[2])
                missing_files.append(generated[3])
                processed += 1
                print('\r', 'Evaluated coverage for {0}/{1} '
                      'inputs with bait_identity={2} and '
                      'bait_coverage={3}.'.format(processed, len(genomes),
                                                  v, report_coverages[i]),
                      end='')

            # save depth values
            depth_files_dir = os.path.join(output_directory, 'depth_files_idnt{0}_cov{1}'.format(str(v).replace('.', ''),
                                                                                           str(report_coverages[i]).replace('.', '')))
            os.mkdir(depth_files_dir)
            depth_files = [mu.write_depth(k, v[3], depth_files_dir) for k, v in final_info.items()]

            test_fig = create_report(coverage_info, final_info, report_dir,
                                     short_samples, ordered_contigs, fixed_xaxis,
                                     fixed_yaxis, ref_set, nr_contigs,
                                     configs, baits_pos, nr_baits+total, nr_baits,
                                     total)

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

    parser.add_argument('-m', '--mode', type=str, required=False,
                        choices=['inclusive', 'exclusive'],
                        default='inclusive',
                        dest='mode',
                        help='')

    parser.add_argument('-mcl', '--minimum-sequence-length', type=int,
                        required=False, default=0,
                        dest='minimum_sequence_length',
                        help='Minimum sequence length. Probes will '
                             'not be created for sequences with a '
                             'length value that is smaller than '
                             'this value.')

    parser.add_argument('-rf', '--refs', type=str,
                        required=True, default=1,
                        dest='refs',
                        help='Number of genome assemblies that will '
                             'be selected to create the initial set '
                             'of probes (the process selects the '
                             'assemblies with least contigs).')

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

    parser.add_argument('-fx', '--fixed-xaxis', required=False,
                        action='store_true',
                        dest='fixed_xaxis',
                        help='')

    parser.add_argument('-fy', '--fixed-yaxis', required=False,
                        action='store_true',
                        dest='fixed_yaxis',
                        help='')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()

    main(**vars(args))
