#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains functions used to create the HTML interactive report.

Code documentation
------------------
"""


from itertools import groupby

import pandas as pd
import datapane as dp
import plotly.graph_objs as go

try:
	import constants as ct
except:
	import proBait.constants as ct


def bait_count(counts):
	"""
	"""
	tracer = go.Bar(x=list(counts.keys()),
					y=list(counts.values()),
					marker={'color': 'rgba(180, 202, 224, 0.7)',
							'line_color': '#045a8d',
							'line_width': 1.5
							},
					showlegend=False)

	return tracer


def depth_hists(depth_values):
	"""
	"""
	tracers = {}
	for k, v in depth_values.items():
		x_values = list(v.values())
		y_values = list(v.keys())
		tracer = go.Bar(x=y_values,
						y=x_values,
						hovertemplate=('<b>Coverage:<b> %{x}'
									   '<br><b>Number of pos.:<b> %{y}'),
						marker={'color': 'rgba(180, 202, 224, 0.7)',
								'line_color': '#045a8d',
								'line_width': 1.5
								},
						showlegend=False)
		tracers[k] = tracer

	return tracers


def missing_intervals_hists(depth_values):
	"""
	"""
	tracers = {}
	for k, v in depth_values.items():
		current_counts = {}
		for c, d in v.items():
			values_groups = [list(j) for i, j in groupby(d[0].values())]
			for n in values_groups:
				if n[0] == 0:
					if len(n) in current_counts:
						current_counts[len(n)] += 1
					else:
						current_counts[len(n)] = 1

		x_values = sorted(list(current_counts.keys()))
		y_values = [current_counts[j] for j in x_values]
		hist_tracer = go.Bar(x=x_values,
							 y=y_values,
							 hovertemplate=('<b>Interval size:<b> %{x}'
											'<br><b>Count:<b> %{y}'),
							 marker={'color': 'rgba(180, 202, 224, 0.7)',
									 'line_color': '#045a8d',
									 'line_width': 1.5
									 },
							 showlegend=False)

		tracers[k] = hist_tracer

	return tracers


def depth_lines(depth_values, ordered_contigs, step=1):
	"""
	"""
	shapes = {}
	tracers = {}
	for k, v in depth_values.items():
		# order contigs based on decreasing length
		contig_order = {}
		for e in ordered_contigs[k]:
			if e[0] in v:
				contig_order[e[0]] = v[e[0]]
			else:
				contig_order[e[0]] = [{i: 0 for i in range(e[1])}]

		x_values = []
		y_values = []
		hovertext = []
		shapes[k] = []
		tracers[k] = []
		# results have 0-based coordinates
		# but plots will be 1-based
		cumulative_pos = 1
		for p, c in contig_order.items():
			# switch to 1-based coordinates for contigs
			contig_pos = 1

			# group depth values into groups of equal sequential values
			values_groups = [list(j) for i, j in groupby(c[0].values())]
			shape_start = cumulative_pos
			j = 0
			for i, g in enumerate(values_groups):
				if i == 0 or i == len(values_groups) - 1 or j >= step:
					# cumulative and contig values already include +1
					# subtract 1 from total sequence length
					hovertext.append(contig_pos)
					hovertext.append(contig_pos + (len(g) - 1))

					start_x = cumulative_pos
					stop_x = start_x + (len(g) - 1)

					x_values.extend([start_x, stop_x])
					y_values.extend([g[0], g[0]])
					if j >= step:
						j = 0

				# Add full length to get start position of next contig
				j += len(g)
				cumulative_pos += len(g)
				contig_pos += len(g)

			shapes[k].append([shape_start, stop_x, p])
		# Use Scattergl to deal with large datasets
		print(len(x_values))
		tracer = go.Scattergl(x=x_values,
							  y=y_values,
							  text=hovertext,
							  hovertemplate=('<b>Contig pos.:<b> %{text}'
											 '<br><b>Cumulative pos.:<b> %{x}'
											 '<br><b>Coverage:<b> %{y}'),
							  showlegend=False,
							  mode='lines',
							  line={'color': '#045a8d',
									'width': 0.5
									},
							  fill='tozeroy')
		tracers[k].append(tracer)

	return [tracers, shapes]


def coverage_table(initial2_data, final2_data, files_info, ref_ids, nr_contigs):
	"""
	"""
	samples = [k+' (ref)'
			   if k in ref_ids
			   else k
			   for k, v in nr_contigs.items()]
	inputs_contigs = [v[0] for k, v in nr_contigs.items()]
	total_lengths = [v[2] for k, v in nr_contigs.items()]

	# initial_cov = [round(initial2_data[k][0], 4) for k in nr_contigs]
	# initial_covered = [initial2_data[k][1] for k in nr_contigs]
	# initial_uncovered = [initial2_data[k][2] for k in nr_contigs]

	# is None when user does not want to generate baits
	if initial2_data is not None:
		generated_baits = [initial2_data[k][3] for k in nr_contigs]
	else:
		generated_baits = [0 for k in nr_contigs]

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

	data = {'Sample': samples,
			'Number of contigs': inputs_contigs,
			'Total length': total_lengths,
			'Breadth of coverage': final_cov,
			'Covered bases': final_covered,
			'Uncovered bases': final_uncovered,
			'Mean depth of coverage': mean_depth,
			# 'Initial breadth of coverage': initial_cov,
			# 'Initial covered bases': initial_covered,
			# 'Initial uncovered bases': initial_uncovered,
			'Generated baits': generated_baits}

	if len(files_info) > 0:
		files_cov = [round(files_info[k][0], 4) for k in nr_contigs]
		files_covered = [files_info[k][1] for k in nr_contigs]
		files_uncovered = [files_info[k][2] for k in nr_contigs]
		data['Files breadth of coverage'] = files_cov
		data['Files covered bases'] = files_covered
		data['Files uncovered bases'] = files_uncovered

		map_mean_depth = []
		for k in nr_contigs:
			length = nr_contigs[k][2]
			depth_counts = files_info[k][4]
			depth_sum = sum([d*c for d, c in depth_counts.items()])
			mean = round(depth_sum/length, 4)
			map_mean_depth.append(mean)

		data['Files mean depth of coverage'] = map_mean_depth

	coverage_df = pd.DataFrame(data)

	table = dp.DataTable(coverage_df)

	return [table, coverage_df]


def create_shape(xref, yref, xaxis_pos, yaxis_pos,
				 line_width=1, dash_type='dashdot'):
	"""
	"""
	shape_tracer = {'type': 'line',
					'xref': xref,
					'yref': yref,
					'x0': xaxis_pos[0],
					'x1': xaxis_pos[1],
					'y0': yaxis_pos[0],
					'y1': yaxis_pos[1],
					'line': {'width': line_width,
							 'dash': dash_type
							 }
					}

	return shape_tracer


def baits_tracer(data, ordered_contigs):
	"""
	"""
	# add baits scatter
	baits_x = []
	baits_y = []
	baits_labels = []
	start = 1
	for contig in ordered_contigs:
		if contig[0] in data:
			# cumulative coordinates
			current_baits = [start+int(n) for n in data[contig[0]]]
			baits_x.extend(current_baits)
			# contig coordinates
			baits_labels.extend([str(int(n)+1) for n in data[contig[0]]])

			y_values = [0] * len(current_baits)
			baits_y.extend(y_values)

		start += contig[1]

	tracer = go.Scattergl(x=baits_x,
						  y=baits_y,
						  mode='markers',
						  marker={'size': 4,
								  'color': '#41ab5d'
								  },
						  showlegend=False,
						  text=baits_labels,
						  hovertemplate=('<b>Contig pos.:<b> %{text}'
										 '<br><b>Cumulative pos.:<b> %{x}'),
						  visible=True)

	return tracer


def create_scatter(x_values, y_values, mode, hovertext):
	"""
	"""
	tracer = go.Scattergl(x=x_values,
						  y=y_values,
						  mode=mode,
						  #line=dict(color='black'),
						  line={'color': 'rgba(147,112,219,0.1)'},
						  showlegend=False,
						  text=hovertext,
						  hovertemplate=('%{text}'),
						  visible=True)

	return tracer


def create_shapes(shapes_data, y_value):
	"""
	"""
	# y_step = int(y_value/4) if int(y_value/4) > 0 else 1
	# hidden_ticks = list(range(1, y_value, y_step))
	hidden_ticks = list(range(1, y_value+1))
	# if y_value not in hidden_ticks:
	# 	hidden_ticks += [y_value]

	shapes_traces = []
	hidden_traces = []
	# Do not create line for last contig or if there is only one contig
	for i, s in enumerate(shapes_data[:-1]):
		# if s != shapes_data[-1]:
		# Only create tracer for end position
		# Start position is equal to end position of previous contig
		shape_tracer = create_shape('x', 'y', [s[1], s[1]], [0, y_value])
		shapes_traces.append(shape_tracer)
		# Create invisible scatter to add hovertext
		hovertext = [s[2], shapes_data[i+1][2]]
		hover_str = '<b><--{0}<b><br><b>{1}--><b>'.format(*hovertext)
		hidden_tracer = create_scatter([s[1]]*len(hidden_ticks),
									hidden_ticks,
									mode='lines',
									hovertext=[hover_str]*y_value)
		hidden_traces.append(hidden_tracer)

	return [shapes_traces, hidden_traces]


def add_summary_text():
	"""
	"""
	summary_text = dp.Text(ct.SUMMARY_TEXT)

	return summary_text
