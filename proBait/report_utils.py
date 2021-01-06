#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""


import os
import math
import random
from itertools import groupby

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

            shapes[k].append([shape_start, stop_x, p])
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


def create_table_tracer(header_values, header_font, header_line, header_fill,
                        cells_values, cells_font, cells_line, cells_fill,
                        domain):
    """
    """

    tracer = go.Table(header=dict(values=header_values,
                                  font=header_font,
                                  align='left',
                                  line=header_line,
                                  fill=header_fill),
                      cells=dict(values=cells_values,
                                 font=cells_font,
                                 align='left',
                                 line=cells_line,
                                 fill=cells_fill),
                      domain=domain)

    return tracer


def coverage_table(initial2_data, final2_data, short_samples,
                   ref_ids, assemblies_lengths):
    """
    """

    ids = {k: v + [short_samples[k]]
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

    table_tracer = create_table_tracer(header_values,
                                       dict(size=16),# family='Courier New'),
                                       dict(color='#ffffff', width=2),
                                       dict(color='#9ecae1'),
                                       cells_values,
                                       dict(size=14),# family='Courier New'),
                                       dict(color='#ffffff', width=1),
                                       dict(color='#f0f0f0'),
                                       dict(x=[0, 1.0]))

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


def create_subplots_fig(nr_rows, nr_cols, titles, specs,
                        shared_yaxes, row_heights):
    """
    """

    fig = make_subplots(rows=nr_rows, cols=nr_cols,
                        subplot_titles=titles,
                        #vertical_spacing=vertical_spacing,
                        #horizontal_spacing=horizontal_spacing,
                        shared_yaxes=shared_yaxes,
                        #column_widths=[0.9, 0.1],
                        specs=specs,
                        row_heights=row_heights)

    return fig


def create_html_report(plotly_fig, output_file, plotlyjs=True):
    """
        
        plotlyjs --> option include True, 'cdn'...check
        https://plotly.com/python/interactive-html-export/

    """

    plot(plotly_fig, filename=output_file,
         auto_open=False, include_plotlyjs=plotlyjs,
         config={"displayModeBar": False, "showTips": False})


def baits_tracer(data, ordered_contigs):
    """
    """

    # add baits scatter
    baits_x = []
    baits_y = []
    baits_labels = []
    start = 0
    for contig in ordered_contigs:
        if contig[0] in data:
            current_baits = [start+int(n) for n in data[contig[0]]]
            baits_x.extend(current_baits)

            baits_labels.extend([str(n) for n in data[contig[0]]])

            #y_values = [random.uniform(0.1, 0.5) for i in range(0, len(current_baits))]
            y_values = [0] * len(current_baits)
            baits_y.extend(y_values)

            start += contig[1]

    tracer = go.Scattergl(x=baits_x, y=baits_y,
                          mode='markers',
                          marker=dict(size=2, color='#02818a'),
                          showlegend=False,
                          text=baits_labels,
                          hovertemplate=('<b>Contig pos.:<b> %{text}'
                                         '<br><b>Cumulative pos.:<b> %{x}'),
                          # set to False so that it is not displayed as
                          # default before selecting dropdown
                          visible=True)

    return tracer


def create_scatter(x_values, y_values, mode, hovertext):
  """
  """

  tracer = go.Scattergl(x=x_values, y=y_values,
                        mode=mode,
                        #line=dict(color='black'),
                        line=dict(color='rgba(147,112,219,0.1)'),
                        showlegend=False,
                        text=hovertext,
                        hovertemplate=('%{text}'),
                        visible=True)

  return tracer
