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

# initial data keys are full paths but it expects basenames!!!
def coverage_table(initial2_data, final2_data, short_samples,
                   ref_ids, nr_contigs):
    """
    """

    samples = [k+' (ref)'
               if k in ref_ids
               else k
               for k, v in nr_contigs.items()]
    inputs_contigs = [v[0] for k, v in nr_contigs.items()]
    total_lengths = [v[2] for k, v in nr_contigs.items()]

    initial_cov = [round(initial2_data[k][0], 4) for k in nr_contigs]
    initial_covered = [initial2_data[k][1] for k in nr_contigs]
    initial_uncovered = [initial2_data[k][2] for k in nr_contigs]

    generated_probes = [initial2_data[k][3] for k in nr_contigs]

    final_cov = [round(final2_data[k][0], 4) for k in nr_contigs]
    final_covered = [final2_data[k][1] for k in nr_contigs]
    final_uncovered = [final2_data[k][2] for k in nr_contigs]

    # determine mean depth of coverage
    mean_depth = []
    for k in nr_contigs:
        length = nr_contigs[k][2]
        depth_counts = final2_data[k][4]
        depth_sum = sum([d*c for d, c in depth_counts.items()])
        mean = round(depth_sum/length, 4)
        mean_depth.append(mean)

    header_values = ['Sample', 'Number of contigs', 'Total length',
                     'Initial breadth of coverage', 'Covered bases',
                     'Uncovered bases', 'Generated probes',
                     'Final breadth of coverage', 'Covered bases',
                     'Uncovered bases', 'Mean depth of coverage']

    cells_values = [samples, inputs_contigs, total_lengths, initial_cov,
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


def report_specs(number_of_inputs):
  """
  """

  specs_def = [[{'type': 'table', 'rowspan': 2, 'colspan': 1},
                  {'type': 'table', 'rowspan': 2, 'colspan': 1}],
                 [None,
                  None],
                 [{'type': 'table', 'rowspan': 2, 'colspan': 2},
                  None],
                 [None,
                  None]] + \
                [[{'type': 'scatter', 'rowspan': 1, 'colspan': 1},
                  {'type': 'bar', 'rowspan': 1, 'colspan': 1}]]*number_of_inputs

  return specs_def


def subplot_titles(inputs_ids):
  """
  """

  titles = [' ', '<b>Configuration</b>', '<b>Coverage statistics</b>']
  for s in inputs_ids:
      titles += ['<b>{0}</b>'.format(s), '']

  return titles


def figure_height(plot_height, table_height, config_height, total_plots):
  """
  """

  total_height = int(plot_height*total_plots + table_height*(total_plots/4) + config_height)
  plots_percentage = round((plot_height*total_plots) / total_height, 2)
  coverage_table_percentage = round((table_height*(total_plots/4)) / total_height, 2)
  summary_table_percentage = round(1 - (plots_percentage+coverage_table_percentage), 2)

  # determine row heights
  plot_height = plots_percentage / total_plots

  row_heights = [summary_table_percentage/2]*2 +\
                [coverage_table_percentage/2]*2 +\
                [plot_height]*(total_plots)

  return [total_height, row_heights]


def adjust_subplot_titles(plotly_fig):
  """
  """

  # adjust configuration table title position and style
  # lock table position to subplot x0 and y1 positions
  subplot12_x = plotly_fig.get_subplot(1, 2).x[0]
  subplot12_y = plotly_fig.get_subplot(1, 2).y[1]
  plotly_fig.layout.annotations[1].update(x=subplot12_x, xref='paper',
                                          xanchor='left', y=subplot12_y,
                                          font=dict(size=18))

  # adjust coverage table
  # lock to subplot x0 and y1
  subplot31_x = plotly_fig.get_subplot(3, 1).x[0]
  subplot31_y = plotly_fig.get_subplot(3, 1).y[1]
  plotly_fig.layout.annotations[2].update(x=subplot31_x, xref='paper',
                                          xanchor='left', y=subplot31_y,
                                          font=dict(size=18))

  # lock depth of coverage plots to paper x0
  for a in plotly_fig.layout.annotations[3:]:
      a.update(x=0, xref='paper', xanchor='left',
               font=dict(size=18))

  return plotly_fig


def add_plots_traces(traces, row, col, top_x, top_y, plotly_fig):
  """
  """

  plotly_fig.add_trace(traces[0], row=row, col=col)
  plotly_fig.update_yaxes(title_text='Coverage', title_font_size=16, row=row, col=col)
  plotly_fig.update_xaxes(title_text='Position', title_font_size=16, domain=[0, 0.9], row=row, col=col)

  # scatter trace with baits start positions
  plotly_fig.add_trace(traces[1], row=row, col=col)

  # add tracer with depth values distribution
  plotly_fig.add_trace(traces[2], row=row, col=col+1)
  plotly_fig.update_yaxes(showticklabels=False, ticks='', row=row, col=col+1)
  plotly_fig.update_xaxes(showticklabels=False, ticks='', zeroline=False, domain=[0.905, 1.0], row=row, col=col+1)

  # adjust axis range
  plotly_fig.update_xaxes(range=[-0.2, top_x], row=row, col=col)
  y_tickvals = list(range(0, top_y, int(top_y/4))) + [top_y]
  plotly_fig.update_yaxes(range=[0-top_y*0.08, top_y+(top_y*0.08)], tickvals=y_tickvals, row=row, col=col)
  plotly_fig.update_yaxes(range=[0-top_y*0.08, top_y+(top_y*0.08)], row=row, col=col+1)

  return plotly_fig


def create_shapes(shapes_data, y_value, ref_axis):
  """
  """

  shapes_traces = []
  hidden_traces = []
  for i, s in enumerate(shapes_data):
      axis_str = '' if ref_axis == 1 else ref_axis
      xref = 'x{0}'.format(axis_str)
      yref = 'y{0}'.format(axis_str)
      # do not create line for last contig
      if s != shapes_data[-1]:
          # only create tracer for end position
          # start position is equal to end position of previous contig
          shape_tracer = create_shape(xref, yref, [s[1], s[1]], [0, y_value])
          shapes_traces.append(shape_tracer)
          # create invisible scatter to add hovertext
          hovertext = [s[2], shapes_data[i+1][2]]
          hover_str = '<b><--{0}<b><br><b>{1}--><b>'.format(*hovertext)
          hidden_ticks = list(range(1, y_value, int(y_value/4)))+[y_value]
          hidden_tracer = create_scatter([s[1]]*len(hidden_ticks),
                                            hidden_ticks,
                                            mode='lines',
                                            hovertext=[hover_str]*y_value)
          hidden_traces.append(hidden_tracer)

  return [shapes_traces, hidden_traces]


def add_plots_titles(plotly_fig):
  """
  """

  annotations_topy = plotly_fig.get_subplot(5, 1).yaxis.domain[1]
  annotations_boty = plotly_fig.get_subplot(5, 1).yaxis.domain[0]
  annotations_y = annotations_topy + (annotations_topy-annotations_boty) / 2.5

  plotly_fig.add_annotation(x=0, xref='paper', xanchor='left',
                     y=annotations_y, yref='paper',
                     yanchor='bottom',
                     text='<b>Depth per position</b>',
                     showarrow=False,
                     font=dict(size=18))

  plotly_fig.add_annotation(x=0.905, xref='paper', xanchor='left',
                     y=annotations_y, yref='paper',
                     yanchor='middle',
                     text='<b>Depth values<br>distribution (log)</b>',
                     showarrow=False,
                     font=dict(size=18))

  return plotly_fig


def add_summary_text(plotly_fig, total_baits, initial_baits, bait_size, bait_offset,
                     iter_baits, total_height):
  """
  """

  # text width is not properly adjusted when screen resolution is different than 1920p
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
                  '      markers along the x-axis are the start positions of baits '
                  'that were generated to cover regions not<br>'
                  '      covered by baits. Contigs are ordered based on decreasing length.<br>'
                  '    <b>- Depth values distribution:</b> distribution of depth of '
                  'coverage values for each input (y-axis is shared with '
                  '<br>      "Depth per position" plot in the same line).<br>'
                  '<br>If you have any question or wish to report an '
                  'issue, please go to proBait\'s '
                  '<a href="https://github.com/rfm-targa/'
                  'proBait">Github</a> repo.').format(total_baits,
                                                      initial_baits,
                                                      bait_size,
                                                      bait_offset,
                                                      iter_baits)

  summary_width = (plotly_fig.get_subplot(1, 1).x[1])*1920
  summary_height = (plotly_fig.get_subplot(1, 1).y[1]-plotly_fig.get_subplot(1, 1).y[0]) * total_height
  plotly_fig.add_annotation(x=0,
                     xref='paper',
                     y=1,
                     yref='paper',
                     text=summary_text,
                     showarrow=False,
                     font=dict(size=16),
                     align='left',
                     bordercolor='#9ecae1',
                     borderwidth=2,
                     bgcolor='#f0f0f0',
                     width=summary_width,
                     height=summary_height)

  return plotly_fig
