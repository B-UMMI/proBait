#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains functions used to map baits against sequences with
minimap2 and compute coverage statistics.

Code documentation
------------------
"""


import os
import math
import subprocess
from copy import deepcopy
from itertools import groupby
from collections import Counter

try:
    import general_utils as gu
except:
    import proBait.general_utils as gu


def run_minimap2(reference, map_fasta, output_file, threads):
    """Execute minimap2 to map sequences to a genome.

    Parameters
    ----------
    reference : str
        Path to the FASTA file with the reference.
    map_fasta : str
        Path to FASTA file with the short sequences
        to map against the reference.
    output_file : str
        Path to the output file with mapping results.

    Returns
    -------
    List with following elements:
        stdout : list
            List with stdout from minimap2 in bytes.
        stderr : list
            List with the stderr from minimpa2 in bytes.
    """
    # -I parameter to control number of target bases loaded into memory
    # --secondary=yes to output secondary alignments that might have high score
    minimap_args = ['minimap2 -I 100M --cs -cx sr --secondary=yes -t {3} {0} {1} > '
                    '{2}'.format(reference, map_fasta, output_file, threads)]

    minimap_proc = subprocess.Popen(minimap_args,
                                    shell=True,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)

    stderr = minimap_proc.stderr.readlines()
    stdout = minimap_proc.stdout.readlines()

    return [stdout, stderr]


def determine_breadth_coverage(intervals, total_bases):
    """Determine the percentage and total number of covered bases.

    Parameters
    ----------
    intervals : dict
        Dictionary with sequence identifiers as keys
        and a list of lists as values. Each sublist has
        a start and stop position in the sequence and
        a dictionary with the coverage for every position
        in the sequence interval.
    total_bases : int
        Total number of bases in the reference.

    Returns
    -------
    List with following elements:
        breadth_coverage : float
            Percentage of covered bases.
        covered_bases : int
            Total number of covered bases.
    """
    # Determine breadth of coverage
    covered_bases = 0
    for k, v in intervals.items():
        for e in v:
            covered_bases += sum([1 for p, c in e[2].items() if c > 0])

    breadth_coverage = covered_bases / total_bases

    return [breadth_coverage, covered_bases]


def determine_depth_coverage(intervals, total_len):
    """Determine depth of coverage for a sequence.

    Parameters
    ----------
    intervals : dict
        Dictionary with sequence identifiers as keys
        and a list of lists as values. Each sublist has
        a start and stop position in the sequence and
        a dictionary with the coverage for every position
        in the sequence interval.
    total_len : int
        Total length of the sequence.

    Returns
    -------
    List with following elements:
        positions_depth : dict
            Dictonary with sequence positions and keys
            and coverage for each position as values.
        counts : dict
            Dictionary with coverage values as keys and
            the total number of positions with that coverage
            value as values.
    """
    # Create dictionary to add coverage for all positions
    positions = list(range(0, total_len))
    positions_depth = {p: 0 for p in positions}
    # Increment coverage values based on intervals
    for i in intervals:
        for p, c in i[2].items():
            positions_depth[p] += c

    # Determine coverage distribution
    counts = sorted(Counter(positions_depth.values()).most_common(),
                    key=lambda x: x[0])

    return [positions_depth, counts]


def single_position_coverage(coverage_info, start, count):
    """Determine if positions in a subsequence are covered based on PAF data.

    Parameters
    ----------
    coverage_info : list
        List with subsequent operations extracted
        from the cd field in a PAF file created by
        minimap2.
    start : int
        Subsequence start position in the complete
        sequence.

    Returns
    -------
    coverage : dict
        Dictionary with sequence positions as keys
        and coverage for each position as values.
    """
    coverage = {}
    for m in coverage_info:
        # Subsequence part with exact matches
        if m[0] == ':':
            # Create dctionary entries with coverage = count
            new_cov = {i: count for i in range(start, start+int(m[1:]))}
            coverage = {**coverage, **new_cov}
            # Increment start position
            start = start + int(m[1:])
        # Position with substitution
        elif m[0] == '*':
            coverage[start] = 0
            start += 1
        # Position with deletion
        elif m[0] == '-':
            # Coverage 0 for missing bases
            new_cov = {i: 0 for i in range(start, start+len(m[1:]))}
            coverage = {**coverage, **new_cov}
            start = start + len(m[1:])
        # Insertion
        elif m[0] == '+':
            # Do not add coverage values for positions because
            # insertion does not exist in reference
            pass

    return coverage


def write_depth(identifier, depth_values, output_dir):
    """
    """
    depth_file = os.path.join(output_dir, identifier+'_depth.tsv')
    depth_lines = []
    for k, v in depth_values.items():
        depth_lines.append(k)
        depth_lines.extend(['{0}\t{1}'.format(p, e) for p, e in v[0].items()])

    gu.write_lines(depth_lines, depth_file)

    return depth_file


def determine_small_bait(span, bait_size, start,
                         stop, sequence_length):
    """Determine baits for regions shorter than bait length.

    Parameters
    ----------
    span : int
        Length of the region that is not covered.
    bait_size : int
        Bait size in bases.
    start : int
        Start position of the subsequence that is
        not covered.
    stop : int
        Stop position of the subsequence that is
        not covered.
    sequence_length : int
        Total length of the sequence.

    Returns
    -------
    bait_interval : list
        List with the start and stop position for
        the determined bait.
    """
    # Bait size is larger than uncovered interval
    # Determine number of bait bases in excess
    rest = bait_size - span
    bot = math.floor(rest / 2)
    top = math.ceil(rest / 2)
    # Uncovered region starts at sequence start
    if start == 0:
        bait_interval = [start, start + bait_size]
    # Centering the bait leads to negative start
    elif (start - bot) < 0:
        bait_interval = [0, top+(bot-start)]
    # Centering the bait leads to stop larger than total sequence
    elif (stop + top) > sequence_length:
        bait_interval = [start-(bot+(top-(sequence_length-stop))), sequence_length]
    # Center bait to get same overlap on both sides
    else:
        bait_interval = [start-bot, stop+top]

    return bait_interval


def determine_interval_baits(bait_size, bait_offset, start, stop, total_len):
    """Determine baits for regions with length equal or greater than bait size.

    Parameters
    ----------
    bait_size : int
        Bait size in bases.
    start : int
        Start position of the subsequence that is
        not covered.
    stop : int
        Stop position of the subsequence that is
        not covered.

    Returns
    -------
    baits : list of list
        List with one sublist per determined bait.
        Each sublist has the start and stop position
        for a bait.
    """
    baits = []
    reach = False
    while reach is False:
        # Remaining interval has length equal or greater than bait size
        if (start + bait_size) >= stop:
            span = stop - start
            bait_interval = determine_small_bait(span, bait_size,
                                                 start, stop, total_len)
            reach = True
        elif (start + bait_size) < stop:
            bait_interval = [start, start+bait_size]
            start = start + bait_offset
        baits.append(bait_interval)

    return baits


def merge_intervals(intervals):
    """Merge intersecting intervals.

    Parameters
    ----------
    intervals : dict
        Dictionary with sequence identifiers as keys
        and a list of lists as values. Each sublist has
        a start and stop position in the sequence and
        a dictionary with the coverage for every position
        in the sequence interval.

    Returns
    -------
    merged : list
        Dictionary with the result of merging intervals
        that overlapped (coverage data is updated and
        incremented for positions in common).
    """

    merged = [deepcopy(intervals[0])]
    for current in intervals[1:]:
        previous = merged[-1]
        # Current and previous intervals intersect
        if current[0] <= previous[1]:
            # Determine top position
            previous[1] = max(previous[1], current[1])
            # Merge coverage dictionaries
            previous_cov = previous[2]
            current_cov = current[2]
            for k, v in current_cov.items():
                if k not in previous_cov:
                    previous_cov[k] = v
                else:
                    previous_cov[k] += v
            previous[2] = previous_cov
        # Current and previous intervals do not intersect
        else:
            merged.append(deepcopy(current))

    return merged


def determine_missing_intervals(intervals, identifier, total_len):
    """Determine sequence intervals that are not covered by any baits.

    Parameters
    ----------
    intervals : dict
        Dictionary with sequence identifiers as keys
        and a list of lists as values. Each sublist has
        a start and stop position in the sequence and
        a dictionary with the coverage for every position
        in the sequence interval.
    identifier : str
        Sequence identifier.
    total_len : int
        Total length of the sequence.

    Returns
    -------
    List with following elements:
        missing_regions : dict
            Dictionary with sequence identifiers as keys
            a list of lists as values. Each sublist has
            the start and stop positions for a sequence
            interval that is not covered by baits.
        not_covered : int
            Total number of bases not covered by baits.
    """
    start = 0
    not_covered = 0
    missing_regions = {identifier: []}
    for i in intervals:
        # Determine if there is a coverage gap between previous and current interval
        diff = i[0] - start
        if diff > 0:
            missing_regions[identifier].append([start, start+diff])
            not_covered += diff
            start += diff

        # Create groups of equal values
        values_groups = [list(j) for i, j in groupby(i[2].values())]
        for g in values_groups:
            if g[0] == 0:
                missing_regions[identifier].append([start, start+len(g)])
                not_covered += len(g)

            start += len(g)

    # Add terminal region
    if start != total_len:
        missing_regions[identifier].append([start, total_len])
        not_covered += total_len - start

    return [missing_regions, not_covered]


# add bait_overlap parameter!
def cover_intervals(intervals, total_len, bait_size,
                    minimum_region, bait_offset):
    """Create baits to cover specified sequence regions.

    Parameters
    ----------
    intervals : list
        List of lists. Each sublist has start and stop
        positions for sequence regions with no coverage.
    total_len : int
        Total length of the sequence.
    bait_size : int
        Bait size in bases.
    minimum_region : int
        Minimum length of the region with no coverage.
        Baits will not be determined to cover regions
        that are shorter than this value.

    Returns
    -------
    cover_baits : list
        List of lists. Each sublist has the start and
        stop positions for a bait.
    """
    cover_baits = []
    for i in intervals:
        start = i[0]
        stop = i[1]
        span = stop - start
        # Check if uncovered region is equal or greater than
        # minimum length value defined for uncovered regions
        if span >= minimum_region:
            if span < bait_size:
                bait_interval = determine_small_bait(span, bait_size,
                                                     start, stop,
                                                     total_len)
                cover_baits.append(bait_interval)
            # Will slide and determine baits
            # if in the last iter, uncovered region is very small it will overlap
            # with regions that are already covered and increase depth of coverage
            elif span >= bait_size:
                # Subtract bait offset to start position to ensure overlap
                # Do not subtract offset if it leads to negative position
                if bait_offset < bait_size and (start-bait_offset) >= 0:
                    start -= bait_offset
                baits = determine_interval_baits(bait_size, bait_offset, start, stop, total_len)
                cover_baits.extend(baits)

    return cover_baits
