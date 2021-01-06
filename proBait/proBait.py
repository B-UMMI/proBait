#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""


import os
import sys
import argparse
from copy import deepcopy

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


def incremental_bait_generator(genomes, unique_baits, output_dir, bait_size,
                               bait_coverage, bait_identity, bait_region,
                               nr_contigs, short_samples, generate=False,
                               depth=False):
    """
    """

    if generate is True:
        header = ('{0:<30}  {1:^10}  {2:^10}  {3:^10}  '
                  '{4:^7}'. format('sample', '%cov', '#cov', '#uncov', '+baits'))
    else:
        header = ('{0:<30}  {1:^10}  {2:^10}  '
                  '{3:^10}'. format('sample', '%cov', '#cov', '#uncov'))

    print('-'*len(header))
    print(header)
    print('-'*len(header))

    total = 0
    discarded_baits_files = []
    missing_files = []
    coverage_info = {}
    for g in genomes:
        invalid_mappings = []
        gbasename = os.path.basename(g).split('.fasta')[0]
        paf_path = os.path.join(output_dir, gbasename+'.paf')

        contigs = gu.import_sequences(g)
        total_bases = nr_contigs[gbasename][2]

        minimap_std = mu.run_minimap2(g, unique_baits, paf_path)

        # read PAF file lines
        paf_lines = gu.read_tabular(paf_path)

        # filter out matches below bait_coverage*bait_size
        # length includes gaps due to deletions and might span more than query length
        valid_length = [line
                        for line in paf_lines
                        if int(line[10]) >= (bait_coverage*bait_size)]

        invalid_mappings.extend([line
                                 for line in paf_lines
                                 if int(line[10]) < (bait_coverage*bait_size)])

        # compute alignment identity
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
        discarded_file = os.path.join(output_dir, '{0}_discarded'.format(gbasename))
        discarded_baits_files.append(discarded_file)
        gu.pickle_dumper(discarded, discarded_file)

        # get information about positions that match to determine coverage
        for i in range(len(valid_pident)):
            current = valid_pident[i][-1]
            start = int(valid_pident[i][7])
            valid_pident[i].append(mu.single_position_coverage(current, start))

        # identify subsequences that are well covered by baits
        covered_intervals = {}
        for l in valid_pident:
            covered_intervals.setdefault(l[5], []).append([int(l[7]), int(l[8]), l[-1]])

        # sort covered intervals
        covered_intervals_sorted = {k: sorted(v, key=lambda x: x[0])
                                    for k, v in covered_intervals.items()}

        # merge overlapping intervals
        # deepcopy to avoid altering original intervals
        merged_intervals = {k: mu.merge_intervals(v)
                            for k, v in covered_intervals_sorted.items()}

        coverage = mu.determine_breadth_coverage(merged_intervals, total_bases)

        # determine subsequences that are not covered
        missing = [mu.determine_missing_intervals(v, k, len(contigs[k]))
                   for k, v in merged_intervals.items()]

        # add missing regions for contigs that had 0 baits mapped
        not_mapped = [[{c: [[0, len(contigs[c])]]}, len(contigs[c])] for c in contigs if c not in merged_intervals]
        missing.extend(not_mapped)

        missing_regions = {k: v for i in missing for k, v in i[0].items()}
        # save missing intervals
        missing_file = os.path.join(output_dir, '{0}_missing'.format(gbasename))
        missing_files.append(missing_file)
        gu.pickle_dumper(missing_regions, missing_file)

        not_covered = sum([i[1] for i in missing])

        coverage_info[gbasename] = [*coverage, not_covered]

        # create baits for missing regions
        if generate is True:
            missing_baits_intervals = {k: mu.cover_intervals(v, len(contigs[k]), bait_size, bait_region)
                                       for k, v in missing_regions.items()}

            extra_probes = {}
            for k, v in missing_baits_intervals.items():
                extra_probes[k] = ['>{0}_{1}\n{2}'.format(k, e[0], contigs[k][e[0]:e[1]]) for e in v]

            new_baits_lines = [v for k, v in extra_probes.items()]
            new_baits_lines = gu.flatten_list(new_baits_lines)

            gu.write_lines(new_baits_lines, unique_baits)

            coverage_info[gbasename].append(len(new_baits_lines))
            total += len(new_baits_lines)

            print('{0:<30}  {1:^10.4f}  {2:^10}  {3:^10}  '
                  '{4:^7}'.format(short_samples[gbasename], coverage[0],
                                  coverage[1], not_covered,
                                  len(new_baits_lines)))

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

            coverage_info[gbasename].extend([depth_values, total_counts])

        if generate is False:
            print('{0:<30}  {1:^10.4f}  {2:^10}  '
                  '{3:^10}'.format(short_samples[gbasename], coverage[0],
                                   coverage[1], not_covered))

    print('-'*len(header))

    return [coverage_info, total, discarded_baits_files, missing_files]


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

    print('Mapping against  and removing similar probes...')
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

##############
############## Baits positions are wrong!!! (offset of -1 in each position?)
initial_data = coverage_info
final_data = final_info
output_dir = report_dir
short_ids = short_samples
# color contig regions that were not covered by probes and that were used
# to generate new probes in different color (add arrows to start and stop)
def create_report(initial_data, final_data, output_dir, short_ids,
                  ordered_contigs, fixed_xaxis, fixed_yaxis, ref_ids,
                  nr_contigs, configs, baits_pos, total_baits,
                  initial_baits, iter_baits):
    """
    """

    # check if user wants equal yaxis ranges for all line plots
    max_x = None
    assemblies_lengths = {k.split('.fasta')[0]: v for k, v in nr_contigs.items()}
    if fixed_xaxis is True:
        max_x = max([v[2] for v in assemblies_lengths.values()])

    max_y = None
    coverage_values = {k: max(list(v[4].keys())) for k, v in final_data[0].items()}
    if fixed_yaxis is True:
        max_y = max(coverage_values.values())

    # create table with run summary
    config_table = ru.create_table_tracer(['Parameter', 'Value'],
                                         dict(size=16),
                                         dict(color='#ffffff', width=2),
                                         dict(color='#9ecae1'),
                                         [list(configs.keys()), list(configs.values())],
                                         dict(size=14),
                                         dict(color='#ffffff', width=1),
                                         dict(color='#f0f0f0'),
                                         dict(x=[0.5, 1.0], y=[0.2, 1.0]))

    table_tracer = ru.coverage_table(initial_data[0], final_data[0], short_ids, ref_ids,
                                     assemblies_lengths)

    # depth of coverage values distribution
    hist_tracers = ru.depth_hists({k: v[4] for k, v in final_data[0].items()})

    # depth of coverage per position
    line_tracers, shapes = ru.depth_lines({k: v[3]
                                           for k, v in final_data[0].items()},
                                          ordered_contigs)

    # line tracers for accepted baits, discarded baits and uncovered regions
#    discarded = {os.path.basename(f).split('_discarded')[0]: f for f in initial_data[2]}
#    merged_discarded = {}
#    for k, v in ordered_contigs.items():
#        start = 0
#        discarded_file = discarded[k]
#        data = gu.pickle_loader(discarded_file)
#        merged_discarded[k] = {}
#        discarded_tracers[k] = []
#        for c in ordered_contigs[k]:
#            if c[0] in data:
#                sorted_intervals = sorted(data[c[0]], key=lambda x: x[0])
#                merged = [deepcopy(sorted_intervals[0])]
#                for i in sorted_intervals[1:]:
#                    previous = merged[-1]
#                    # current and previous intervals intersect
#                    if i[0] <= previous[1]:
#                        # determine top position
#                        previous[1] = max(previous[1], i[1])
#                        previous[2].extend(i[2])
#                    else:
#                        merged.append(deepcopy(i))
#
#                merged_discarded[k][c[0]] = [[i[0]+start, i[1]+start, i[2]] for i in merged]
#                start += c[1]


    # baits start position per input
    baits_tracers = {k: ru.baits_tracer(baits_pos[k.split('.')[0]], v)
                     for k, v in ordered_contigs.items() if k.split('.')[0] in baits_pos}

    nr_rows = len(line_tracers) + 4
    titles = [' ', '<b>Configuration</b>', '<b>Coverage statistics</b>']
    for s in list(short_ids.values()):
        titles += ['<b>{0}</b>'.format(s), '']

    specs_def = [[{'type': 'table', 'rowspan': 2, 'colspan': 1},
                  {'type': 'table', 'rowspan': 2, 'colspan': 1}],
                 [None,
                  None],
                 [{'type': 'table', 'rowspan': 2, 'colspan': 2},
                  None],
                 [None,
                  None]] + \
                [[{'type': 'scatter', 'rowspan': 1, 'colspan': 1},
                  {'type': 'bar', 'rowspan': 1, 'colspan': 1}]]*len(line_tracers)

    # define figure height
    height = int(190*len(line_tracers) + 150*(len(line_tracers)/4) + 505)
    plots_percentage = round((190*len(line_tracers)) / height, 2)
    coverage_table_percentage = round((150*(len(line_tracers)/4)) / height, 2)
    summary_table_percentage = round(1 - (plots_percentage+coverage_table_percentage), 2)

    # determine row heights
    plot_height = plots_percentage / len(line_tracers)

    row_heights = [summary_table_percentage/2]*2 +\
                  [coverage_table_percentage/2]*2 +\
                  [plot_height]*(len(line_tracers))

    # adding vertical_spacing in combination with row_heights can
    # exceed valid values and lead to error
    fig = ru.create_subplots_fig(nr_rows, 2, titles, specs_def, shared_yaxes=True,
                                 row_heights=row_heights)

    # change subplots titles positions
    subplot12_x = fig.get_subplot(1, 2).x[0]
    subplot12_y = fig.get_subplot(1, 2).y[1]
    fig.layout.annotations[1].update(x=subplot12_x, xref='paper',
                                     xanchor='left', y=subplot12_y,
                                     font=dict(size=18))

    # get topmost y coordinate for subplot with coverage stats
    subplot31_x = fig.get_subplot(3, 1).x[0]
    subplot31_y = fig.get_subplot(3, 1).y[1]
    fig.layout.annotations[2].update(x=subplot31_x, xref='paper',
                                     xanchor='left', y=subplot31_y,
                                     font=dict(size=18))

    for a in fig.layout.annotations[3:]:
        a.update(x=0, xref='paper', xanchor='left',
                 font=dict(size=18))

    fig.add_trace(config_table, row=1, col=2)

    # add tracer with coverage stats
    fig.add_trace(table_tracer, row=3, col=1)

    r = 5
    c = 1
    for k, v in line_tracers.items():
        # add tracer with depth per position
        fig.add_trace(v[0], row=r, col=c)
        fig.update_yaxes(title_text='Coverage', title_font_size=16, row=r, col=c)
        fig.update_xaxes(title_text='Position', title_font_size=16, domain=[0, 0.9], row=r, col=c)

        # add tracer with baits start position
        fig.add_trace(baits_tracers[k], row=r, col=c)

        # add tracer with depth values distribution
        fig.add_trace(hist_tracers[k], row=r, col=c+1)
        fig.update_yaxes(showticklabels=False, ticks='', row=r, col=c+1)
        fig.update_xaxes(showticklabels=False, ticks='', zeroline=False, domain=[0.905, 1.0], row=r, col=c+1)

        top_x = assemblies_lengths[k] if max_x is None else max_x
        top_y = coverage_values[k] if max_y is None else max_y

        # adjust axis range
        fig.update_xaxes(range=[-0.2, top_x], row=r, col=c)
        y_tickvals = list(range(0, top_y, int(top_y/4))) + [top_y]
        fig.update_yaxes(range=[0-top_y*0.08, top_y+(top_y*0.08)], tickvals=y_tickvals, row=r, col=c)
        fig.update_yaxes(range=[0-top_y*0.08, top_y+(top_y*0.08)], row=r, col=c+1)

        r += 1

    # create shapes for contig boundaries
    ref_axis = 1
    shapes_tracers = []
    hidden_tracers = []
    start_hidden = 5
    for k, v in shapes.items():
        current_shapes = list(shapes[k])
        y_value = coverage_values[k] if max_y is None else max_y
        for i, s in enumerate(current_shapes):
            axis_str = '' if ref_axis == 1 else ref_axis
            xref = 'x{0}'.format(axis_str)
            yref = 'y{0}'.format(axis_str)
            # do not create line for last contig
            if s != current_shapes[-1]:
                # only create tracer for end position
                # start position is equal to end position of previous contig
                shape_tracer = ru.create_shape(xref, yref, [s[1], s[1]], [0, y_value])
                # create invisible scatter to add hovertext
                hovertext = [s[2], current_shapes[i+1][2]]
                hover_str = '<b><--{0}<b><br><b>{1}--><b>'.format(*hovertext)
                hidden_ticks = list(range(1, top_y, int(top_y/4)))+[top_y]
                hidden_tracer = ru.create_scatter([s[1]]*len(hidden_ticks),
                                                  hidden_ticks,
                                                  mode='lines',
                                                  hovertext=[hover_str]*y_value)
                hidden_tracers.append([hidden_tracer, start_hidden, 1])
                shapes_tracers.append(shape_tracer)

        ref_axis += 2
        start_hidden += 1

    fig.update_layout(shapes=shapes_tracers, clickmode='event')

    # add hidden tracers
    for t in hidden_tracers:
        fig.add_trace(t[0], t[1], t[2])

    # disable grid
    fig.update_xaxes(showgrid=False)
    fig.update_yaxes(showgrid=False, showline=True, linecolor='black')

    # bars with distribution of depth values in logscale
    for i in range(5, 5+len(line_tracers)):
        fig.update_xaxes(type='log', row=i, col=2)

    annotations_topy = fig.get_subplot(5, 1).yaxis.domain[1]
    annotations_boty = fig.get_subplot(5, 1).yaxis.domain[0]
    annotations_y = annotations_topy + (annotations_topy-annotations_boty) / 2.5

    fig.add_annotation(x=0, xref='paper', xanchor='left',
                       y=annotations_y, yref='paper',
                       yanchor='bottom',
                       text='<b>Depth per position</b>',
                       showarrow=False,
                       font=dict(size=18))

    fig.add_annotation(x=0.905, xref='paper', xanchor='left',
                       y=annotations_y, yref='paper',
                       yanchor='middle',
                       text='<b>Depth values<br>distribution (log)</b>',
                       showarrow=False,
                       font=dict(size=18))

    # add summary text
    summary_text = ('Generated a total of <b>{0}</b> baits.<br>'
                    'An initial set of {1} baits was generated by decomposing 1 '
                    'reference(s) into baits of size {2}bps with an offset of {3}.'
                    '<br>An additional set of {4} baits were generated through '
                    'the iterative process of mapping the initial set of '
                    'baits against all<br>inputs and generating new baits for regions that '
                    'did not have any mapped baits or had mapped '
                    'baits without sufficient<br>identity and/or coverage.<br>'
                    '<br>The report has the following sections:<br>'
                    '    <b>- Configuration:</b> values passed to proBait\'s parameters.<br>'
                    '    <b>- Coverage statistics:</b> coverage statistics determined by mapping '
                    'the final set of baits against each input.<br>'
                    '    <b>- Depth per position:</b> depth of coverage per position. Vertical '
                    'dashed lines are contig boundaries and green<br>'
                    '      markers along the xaxis zeroline are the start positions of baits '
                    'that were generated to cover regions not<br>'
                    '      covered by baits. Contigs are ordered based on decreasing length.<br>'
                    '    <b>- Depth values distribution:</b> distribution of depth of '
                    'coverage values for each input (yaxis is shared with '
                    '<br>      "Depth per position" plot in the same line).<br>'
                    '<br>If you have any question or wish to report an '
                    'issue, please go to proBait\'s '
                    '<a href="https://github.com/rfm-targa/'
                    'proBait">Github</a> repo.').format(total_baits,
                                                        initial_baits,
                                                        configs['Bait size'],
                                                        configs['Bait offset'],
                                                        iter_baits)

    run_summary = fig.add_annotation(x=0,
                                     xref='paper',
                                     y=1,
                                     yref='paper',
                                     text=summary_text,
                                     showarrow=False,
                                     font=dict(size=16),
                                     align='left',
                                     bordercolor='#9ecae1',
                                     borderwidth=2,
                                     bgcolor='#f0f0f0')

    fig.update_layout(title=dict(text='<b>proBait - Coverage Report</b>',
                                 x=0,
                                 xref='paper',
                                 xanchor='left',
                                 font=dict(size=30)))

    # line plots need fixed space
    fig.update_layout(height=height, width=1920,
                      template='ggplot2',
                      #paper_bgcolor='#f0f0f0')
                      plot_bgcolor='rgba(0,0,0,0)')

    output_plot = os.path.join(output_dir, 'report.html')
    ru.create_html_report(fig, output_plot)

    return fig


#input_files = '/home/rfm/Desktop/rfm/Lab_Analyses/pneumo_baits_design/assemblies'
input_files = '/home/rfm/Desktop/rfm/Lab_Analyses/pneumo_baits_design/test_assemblies'
#input_files = '/home/rfm/Desktop/rfm/Lab_Analyses/pneumo_baits_design/single_assembly'
#input_files = '/home/rfm/Desktop/rfm/Lab_Analyses/pneumo_baits_design/ref32'
output_dir = '/home/rfm/Desktop/rfm/Lab_Analyses/pneumo_baits_design/tmp'
bait_size = 120
bait_offset = 120
number_refs = 1
bait_identity = 1.0
bait_coverage = 1.0
bait_region = 10
cluster_probes = False
cluster_identity = 0.8
cluster_coverage = 0.9
minlen_contig = 120
#exclude_regions = None
exclude_regions = '/home/rfm/Desktop/rfm/Lab_Analyses/pneumo_baits_design/ncbi-genomes-2020-11-16/GCF_000001405.39_GRCh38.p13_genomic.fna'
exclude_pident = 0.7
exclude_coverage = 0.7
threads = 4
fixed_xaxis = True
fixed_yaxis = True
contig_boundaries = 100
report = True


# Add features to control depth of coverage of regions.
# e.g.: duplicate coverage of regions only covered once.

# determine length of subsequences with 0 coverage and etc
# number of SNPs, deletions, insertions and etc

# add option to avoid generating baits close to contigs boundaries

# add option to determine baits from target regions and then only generate baits
# to capture diversity in those regions in other genomes (this allows to determine baits
# only for targeted regions and will not generate baits for other uncovered loci)
# users can provide annotation labels for the target regions that are included in
# the hovertext
# important to search for options to control mmseqs2 and minimap2 memory usage
# the process might fail because those tools use too much memory
# for the step that maps against the human genome, give 1 chromossome at a time
# possible to change info shown in table with dropdown as filters?
# option to receive baits as input and add more baits generated based on input genomes
# final table with number of uncovered regions per sample, number of SNPs, number of deletions, insertions, etc.
# add boxplot with distribution of identity values for mapped baits against each genome
# color baits according to feature that led to bait creation (insertion, deletion, etc)
# detect missing regions that are close and can be included in same probe.
# plot with markers representing baits should also have the following lines:
# - Lines for the intervals of baits that were accepted (green?)
# - Lines for intervals of baits that were not accepted (red?)
# - Lines for regions that had no baits mapped
# The lines for baits that were not accepted should also have information
# about why the baits were not accepted (include SNO, deletion, insertion
# info somehow...). Adding these lines should not increase report file size
# that much because each interval will only have 2 points to represent
# accepted baits, refused baits or uncovered regions.
def main(input_files, output_dir, mode, minlen_contig, contig_boundaries,
         number_refs, bait_size, bait_offset, bait_identity, bait_coverage,
         bait_region, cluster_probes, cluster_identity, cluster_coverage,
         exclude_regions, exclude_pident, exclude_coverage, threads,
         report, fixed_xaxis, fixed_yaxis):

    if os.path.isdir(output_dir) is False:
        os.mkdir(output_dir)
    else:
        sys.exit('Output directory exists. Please provide a path '
                 'for a directory that will be created to store files.')

    genomes = [os.path.join(input_files, file)
               for file in os.listdir(input_files)]

    # determine number of contigs and total length
    nr_contigs = {f: gu.count_contigs(f, minlen_contig) for f in genomes}

    # select assemblies with lowest number of contigs
    sorted_contigs = sorted(list(nr_contigs.items()), key=lambda x: x[1][1])
    ref_set = [t[0] for t in sorted_contigs[0:number_refs]]

    genomes = ref_set + [g for g in genomes if g not in ref_set]
    nr_contigs = {os.path.basename(f).split('.fasta')[0]: gu.count_contigs(f, minlen_contig) for f in genomes}

    # get short identifiers
    short_samples = gu.common_suffixes([os.path.basename(f).split('.fasta')[0] for f in genomes])

    # shred genomic sequences
    # not generating kmers that cover the end of the sequences!
    baits_file = os.path.join(output_dir, 'baits.fasta')
    for g in ref_set:
        nr_baits = gu.generate_baits(g, baits_file, bait_size,
                                     bait_offset, minlen_contig)

    print('\nCreated initial set of {0} probes based on {1} '
          'assemblies.'.format(nr_baits, number_refs))

    # identify unique baits
    unique_baits = os.path.join(output_dir, 'unique_baits.fasta')
    total, unique_seqids = gu.determine_distinct(baits_file, unique_baits)
    print('Removed {0} repeated probes.\n'.format(total))

    # start mapping baits against remaining genomes
    # mapping against ref_set to cover missing regions
    bait_creation_dir = os.path.join(output_dir, 'incremental_bait_creation')
    os.mkdir(bait_creation_dir)
    coverage_info = incremental_bait_generator(genomes, unique_baits,
                                               bait_creation_dir, bait_size,
                                               bait_coverage, bait_identity,
                                               bait_region, nr_contigs,
                                               short_samples,
                                               generate=True, depth=False)

    print('Added {0} probes to cover {1} assemblies.\nTotal '
          'of {2} probes.'.format(coverage_info[1], len(genomes), nr_baits+coverage_info[1]))

    if cluster_probes is True:
        clustering_dir = os.path.join(output_dir, 'clustering')
        os.mkdir(clustering_dir)
        unique_baits, removed = cu.exclude_similar_probes(unique_baits, clustering_dir,
                                                          cluster_identity, cluster_coverage,
                                                          bait_size, threads)

    if exclude_regions is not None:
        exclude_dir = os.path.join(output_dir, 'exclude')
        os.mkdir(exclude_dir)
        unique_baits, removed = exclude_contaminant(unique_baits, exclude_regions,
                                                    exclude_pident, exclude_coverage,
                                                    bait_size, exclude_dir)

    # determine breadth of coverage for all assemblies
    # and depth of coverage for each base
    final_coverage_dir = os.path.join(output_dir, 'final_coverage')
    os.mkdir(final_coverage_dir)
    final_info = incremental_bait_generator(genomes, unique_baits, final_coverage_dir,
                                            bait_size, bait_coverage,
                                            bait_identity, bait_region,
                                            nr_contigs, short_samples,
                                            generate=False, depth=True)

    # save depth values
    depth_files_dir = os.path.join(output_dir, 'depth_files')
    os.mkdir(depth_files_dir)
    depth_files = [mu.write_depth(k, v[3], depth_files_dir) for k, v in final_info[0].items()]

    if report is True:
        # determine contig order from longest to shortest
        ordered_contigs = gu.order_contigs(genomes)

        # create plots
        report_dir = os.path.join(output_dir, 'report_data')
        os.mkdir(report_dir)

        # determine number of exclude regions and total bps
        exclude_stats = [len(rec) for rec in SeqIO.parse(exclude_regions, 'fasta')]
        total_bps = sum(exclude_stats)

        # get baits positions in genomes
        baits_records = SeqIO.parse(unique_baits, 'fasta')
        baits_pos = {s: {} for s in short_samples.values()}
        total_baits = 0
        for rec in baits_records:
            genome = (rec.id).split('_')[0]
            pos = (rec.id).split('_')[-1]
            contig = (rec.id).split('_{0}'.format(pos))[0]
            baits_pos[genome].setdefault(contig, []).append(pos)
            total_baits += 1

        # create dict with config values
        configs = {'Number of inputs': len(genomes),
                   'Minimum contig length': minlen_contig,
                   'Contig boundaries distance': contig_boundaries,
                   'Number of references': number_refs,
                   'Bait size': bait_size,
                   'Bait offset': bait_offset,
                   'Bait identity': bait_identity,
                   'Bait coverage': bait_coverage,
                   'Bait region': bait_region,
                   'Cluster probes': str(cluster_probes),
                   'Cluster identity': cluster_identity,
                   'Cluster coverage': cluster_coverage,
                   'Exclude regions': '{0} regions'.format(len(exclude_stats)),
                   'Exclude identity': exclude_pident,
                   'Exclude coverage': exclude_coverage}

        ref_ids = [os.path.basename(f).split('.fasta')[0] for f in ref_set]
        test_fig = create_report(coverage_info, final_info, report_dir,
                                 short_samples, ordered_contigs, fixed_xaxis,
                                 fixed_yaxis, ref_ids, nr_contigs,
                                 configs, baits_pos, total_baits, nr_baits,
                                 coverage_info[1])

        print('Coverage report available in {0}'.format(report_dir))

    print('Created a set of {0} probes'.format())


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', type=str, required=True,
                        dest='input_files',
                        help='Path to the directory with '
                             'input FASTA files.')

    parser.add_argument('-o', type=str, required=True,
                        dest='output_dir',
                        help='Path to the output directory where '
                             'files will be saved to (must not exist, '
                             'process will create this directory).')

    parser.add_argument('--m', type=str, required=False,
                        choices=['inclusive', 'exclusive'],
                        default='inclusive',
                        dest='mode',
                        help='')

    parser.add_argument('--mc', type=int, required=False,
                        default=360,
                        dest='minlen_contig',
                        help='Minimum contig length. Probes will '
                             'not be created for contigs with a '
                             'length value that is smaller than '
                             'this value.')

    parser.add_argument('--cb', type=int, required=False,
                        default=0,
                        dest='contig_boundaries',
                        help='Distance to contig boundaries. '
                             'Baits will be determined for uncovered '
                             'regions that are at least this value of '
                             'bases from one of the contig boundaries.')

    parser.add_argument('--nr', type=int, required=False,
                        default=1,
                        dest='number_refs',
                        help='Number of genome assemblies that will '
                             'be selected to create the initial set '
                             'of probes (the process selects the '
                             'assemblies with least contigs).')

    parser.add_argument('--bs', type=int, required=False,
                        default=120,
                        dest='bait_size',
                        help='Length of the probes that the process '
                             'will create.')

    parser.add_argument('--bo', type=int, required=False,
                        default=120,
                        dest='bait_offset',
                        help='Start position offset between consecutive '
                             'probes.')

    parser.add_argument('--bi', type=float, required=False,
                        default=1.0,
                        dest='bait_identity',
                        help='Minimum percent identity that '
                             'aligned probes need to have when '
                             'aligned against input genomes ('
                             'probes with lower identity are '
                             'not included in the set of probes '
                             'that cover well regions of input '
                             'genomes).')

    parser.add_argument('--bc', type=float, required=False,
                        default=1.0,
                        dest='bait_coverage',
                        help='Minimum percent length of the '
                             'probe that has to align against the '
                             'genome assembly (probes with lower '
                             'coverage are not included in the set '
                             'of probes that cover well regions of '
                             'input genomes).')

    parser.add_argument('--br', type=int, required=False,
                        default=0,
                        dest='bait_region',
                        help='Uncovered regions must have a length '
                             'value equal or greater than this value. '
                             'If the uncovered region is smaller than '
                             'this value the process will not generate '
                             'new baits to cover that region.')

    parser.add_argument('--c', required=False, action='store_true',
                        dest='cluster_probes',
                        help='Cluster set of probes after generating '
                             'probes to cover all input assemblies. '
                             'This clustering step will cluster '
                             'similar probes and remove highly similar '
                             'probes based on percent identity and '
                             'coverage.')

    parser.add_argument('--ci', type=float, required=False,
                        default=1.0,
                        dest='cluster_identity',
                        help='Clustered probes with equal or higher '
                             'percent identity are excluded.')

    parser.add_argument('--cc', type=float, required=False,
                        default=1.0,
                        dest='cluster_coverage',
                        help='Clustered probes with equal or higher '
                             'coverage may be excluded based on '
                             'percent identity.')

    parser.add_argument('--e', type=str, required=False,
                        default=None,
                        dest='exclude_regions',
                        help='Path to a FASTA file with genomic regions '
                             'that probes must not cover.')

    parser.add_argument('--ep', type=float, required=False,
                        default=0.8,
                        dest='exclude_pident',
                        help='Probes with percent identity equal or '
                             'higher than this value to regions that '
                             'must not be covered will be excluded.')

    parser.add_argument('--ec', type=float, required=False,
                        default=0.5,
                        dest='exclude_coverage',
                        help='Probes that map against the regions to '
                             'exclude with equal or greater coverage '
                             'may be excluded based on percent identity.')

    parser.add_argument('--t', type=int, required=False,
                        default=1, dest='threads',
                        help='')

    parser.add_argument('--report', required=False, action='store_true',
                        dest='report',
                        help='')

    parser.add_argument('--fx', required=False, action='store_true',
                        dest='fixed_xaxis',
                        help='')

    parser.add_argument('--fy', required=False, action='store_true',
                        dest='fixed_yaxis',
                        help='')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()

    main(**vars(args))
