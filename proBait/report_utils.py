#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""


import plotly.graph_objs as go
from plotly.offline import plot
from plotly.subplots import make_subplots


def depth_hists(depth_values):
    """
    """

    tracers = {}
    for k, v in depth_values.items():
        x_values = list(v.values())
        y_values = list(v.keys())
        tracer = go.Bar(x=x_values,
                        y=y_values,
                        hovertemplate=('<b>Coverage:<b> %{y}'
                                      '<br><b>Number of pos.:<b> %{x}'),
                        marker=dict(color='#67a9cf'),
                        showlegend=False,
                        orientation='h')
        tracers[k] = tracer

    return tracers


def depth_lines(depth_values, ordered_contigs):
    """
    """

    tracers = {}
    shapes = {}
    for k, v in depth_values.items():
        x_values = []
        y_values = []
        hovertext = []
        # start genome at xaxis=1 in plot
        cumulative_pos = 1
        contig_order = {}
        shapes[k] = []
        tracers[k] = []
        for e in ordered_contigs[k]:
            if e[0] in v:
                contig_order[e[0]] = v[e[0]]
            else:
                contig_order[e[0]] = [{i: 0 for i in range(e[1])}]

        for p, c in contig_order.items():
            contig_pos = 1
            values_groups = [list(j) for i, j in groupby(c[0].values())]
            shape_start = cumulative_pos
            for g in values_groups:
                hovertext.append(contig_pos)
                hovertext.append(contig_pos + (len(g) - 1))

                start_x = cumulative_pos
                stop_x = start_x + (len(g) - 1)

                cumulative_pos += len(g)
                contig_pos += len(g)

                x_values.extend([start_x, stop_x])
                y_values.extend([g[0], g[0]])

            shapes[k].append([shape_start, stop_x])
        # use Scattergl to deal with large datasets
        tracer = go.Scattergl(x=x_values,
                              y=y_values,
                              text=hovertext,
                              hovertemplate=('<b>Contig pos.:<b> %{text}'
                                             '<br><b>Cumulative pos.:<b> %{x}'
                                             '<br><b>Coverage:<b> %{y}'),
                              showlegend=False,
                              mode='lines',
                              line=dict(color='#3690c0', width=0.5),
                              fill='tozeroy')
        tracers[k].append(tracer)

    return [tracers, shapes]


def create_table_tracer(header_values, cells_values, domain):
    """
    """

    tracer = go.Table(header=dict(values=header_values,
                                  font=dict(size=12),
                                  align='left'),
                      cells=dict(values=cells_values,
                                 align='left'),
                      domain=domain)

    return tracer


def coverage_table(initial2_data, final2_data, short_samples, ref_ids,
                   assemblies_lengths):
    """
    """

    ids = {os.path.basename(k): v + [short_samples[k+'.fasta']]
           for k, v in assemblies_lengths.items()}

    samples = [v[3]+' (ref)'
               if k in ref_ids
               else v[3]
               for k, v in ids.items()]
    nr_contigs = [v[0] for k, v in ids.items()]
    total_lengths = [v[2] for k, v in ids.items()]

    initial_cov = [round(initial2_data[k][0], 4) for k in ids]
    initial_covered = [initial2_data[k][1] for k in ids]
    initial_uncovered = [initial2_data[k][2] for k in ids]

    generated_probes = [initial2_data[k][3] for k in ids]

    final_cov = [round(final2_data[k][0], 4) for k in ids]
    final_covered = [final2_data[k][1] for k in ids]
    final_uncovered = [final2_data[k][2] for k in ids]

    # determine mean depth of coverage
    mean_depth = []
    for k in ids:
        length = ids[k][2]
        depth_counts = final2_data[k][4]
        depth_sum = sum([d*c for d, c in depth_counts.items()])
        mean = round(depth_sum/length, 4)
        mean_depth.append(mean)

    header_values = ['Sample', 'Number of contigs', 'Total length',
                     'Initial breadth of coverage', 'Covered bases',
                     'Uncovered bases', 'Generated probes',
                     'Final breadth of coverage', 'Covered bases',
                     'Uncovered bases', 'Mean depth of coverage']

    cells_values = [samples, nr_contigs, total_lengths, initial_cov,
                    initial_covered, initial_uncovered, generated_probes,
                    final_cov, final_covered, final_uncovered, mean_depth]

    table_tracer = create_table_tracer(header_values, cells_values, dict(x=[0, 1]))

    return table_tracer


def create_shape(xref, yref, xaxis_pos, yaxis_pos,
                 line_width=1, dash_type='dashdot'):
    """
    """

    shape_tracer = dict(type='line',
                        xref=xref,
                        yref=yref,
                        x0=xaxis_pos[0], x1=xaxis_pos[1],
                        y0=yaxis_pos[0], y1=yaxis_pos[1],
                        line=dict(width=line_width,
                                  dash=dash_type))

    return shape_tracer


# color contig regions that were not covered by probes and that were used
# to generate new probes in different color (add arrows to start and stop)
def create_report(initial_data, final_data, output_dir, short_ids,
                  ordered_contigs, fixed_xaxis, fixed_yaxis, ref_ids,
                  nr_contigs, configs):
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

    table_tracer = coverage_table(initial_data[0], final_data[0], short_ids, ref_ids,
                                  assemblies_lengths)

    # depth of coverage values distribution
    hist_tracers = depth_hists({k: v[4] for k, v in final_data[0].items()})

    # depth of coverage per position
    line_tracers, shapes = depth_lines({k: v[3]
                                        for k, v in final_data[0].items()},
                                       ordered_contigs)

    nr_rows = len(line_tracers) + 4
    titles = ['Configs', 'Scatter', 'Coverage statistics']
    for s in list(short_ids.values()):
        titles += [s, '']

    specs_def = [[{'type': 'table', 'rowspan': 2, 'colspan': 1}, {'type': 'scatter', 'rowspan': 2, 'colspan': 1}],
                 [None, None],
                 [{'type': 'table', 'rowspan': 2, 'colspan': 2}, None],
                 [None, None]]+[[{'type': 'scatter'}, {'type': 'bar'}]]*len(line_tracers)

    fig = make_subplots(rows=nr_rows, cols=2,
                        subplot_titles=titles,
                        horizontal_spacing=0.002,
                        shared_yaxes=True,
                        #column_widths=[0.9, 0.1],
                        specs=specs_def)

    # change subplots titles positions
    # lock/link table subplots titles to xaxis2 to force fixed position
    fig.layout['annotations'][0]['x'] = 0
    fig.layout['annotations'][0]['xref'] = 'x2'
    fig.layout['annotations'][0]['xanchor'] = 'left'

    fig.layout['annotations'][2]['x'] = 0
    fig.layout['annotations'][2]['xref'] = 'x2'
    fig.layout['annotations'][2]['xanchor'] = 'left'

    # change title of first scatter
    fig.layout['annotations'][1]['x'] = 0
    fig.layout['annotations'][1]['xref'] = 'x1'
    fig.layout['annotations'][1]['xanchor'] = 'left'

    x = 2
    for a in fig.layout['annotations'][3:]:
        a['x'] = 0
        a['xref'] = 'x{0}'.format(x)
        a['xanchor'] = 'left'
        x += 2

    print(fig.layout)
    # create table with run summary
    run_summary = create_table_tracer(['Parameter', 'Value'], [list(configs.keys()), list(configs.values())], dict(x=[0, 0.5]))
    fig.add_trace(run_summary, row=1, col=1)

    # add empty scatter
    empty_tracer = go.Scatter(x=[1], y=[1],
                              showlegend=False)
    fig.add_trace(empty_tracer, row=1, col=2)
    # update domain of first scatter
    fig.update_xaxes(domain=[0.53, 1.0], row=1, col=2)

    # add tracer with coverage stats
    fig.add_trace(table_tracer, row=3, col=1)

    r = 5
    c = 1
    for k, v in line_tracers.items():
        fig.add_trace(v[0], row=r, col=c)
        fig.update_yaxes(title_text='Coverage', row=r, col=c)
        fig.update_xaxes(title_text='Position', domain=[0, 0.9], row=r, col=c)

        fig.add_trace(hist_tracers[k], row=r, col=c+1)
        fig.update_yaxes(showticklabels=False, ticks='', row=r, col=c+1)
        fig.update_xaxes(showticklabels=False, ticks='', domain=[0.905, 1.0], row=r, col=c+1)

        top_x = assemblies_lengths[k] if max_x is None else max_x
        top_y = coverage_values[k] if max_y is None else max_y

        # adjust axis range
        fig.update_xaxes(range=[-0.2, top_x], row=r, col=c)
        fig.update_yaxes(range=[0-top_y*0.08, top_y+(top_y*0.08)], row=r, col=c)
        fig.update_yaxes(range=[0-top_y*0.08, top_y+(top_y*0.08)], row=r, col=c+1)

        r += 1

    # create shapes for contig boundaries
    ref_axis = 2
    shapes_tracers = []
    for k, v in shapes.items():
        current_shapes = list(shapes[k])
        y_value = coverage_values[k] if max_y is None else max_y
        for s in current_shapes:
            axis_str = '' if ref_axis == 1 else ref_axis
            xref = 'x{0}'.format(axis_str)
            yref = 'y{0}'.format(axis_str)
            # do not create line for last contig
            if s != current_shapes[-1]:
                # only create tracer for end position
                # start position is equal to end position of previous contig
                shape_tracer = create_shape(xref, yref, [s[1], s[1]], [0, y_value])
                shapes_tracers.append(shape_tracer)

        ref_axis += 2

    fig.update_layout(shapes=shapes_tracers, clickmode='event')

    # disable grid
    fig.update_xaxes(showgrid=False)
    fig.update_yaxes(showgrid=False)

    # change titles position
    for annotation in fig.layout.annotations:
        annotation.update(x=0.03)

    # bars with distribution of depth values in logscale
    for i in range(5, 5+len(line_tracers)):
        fig.update_xaxes(type='log', row=i, col=2)

    # line plots need fixed space
    fig.update_layout(title='proBait - Coverage Report',
                      height=200*len(line_tracers)+400,
                      template='ggplot2')#,
                      #paper_bgcolor='rgba(0,0,0,0)',
                      #plot_bgcolor='rgba(0,0,0,0)')  # plotly_dark, presentation+ggplot2

    output_plot = os.path.join(output_dir, 'report.html')
    plot(fig, filename=output_plot, auto_open=False)