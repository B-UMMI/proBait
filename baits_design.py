#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""


import os
import re
import sys
import csv
import math
import shutil
import pickle
import hashlib
import argparse
import itertools
import subprocess
from copy import deepcopy
from itertools import groupby
from collections import Counter

from Bio import SeqIO
import plotly.graph_objs as go
from plotly.offline import plot
from plotly.subplots import make_subplots


def pickle_dumper(content, output_file):
    """ Use the Pickle module to serialize an object.

        Parameters
        ----------
        content : type
            Variable that refers to the object that will
            be serialized and written to the output file.
        output_file : str
            Path to the output file.
    """

    with open(output_file, 'wb') as po:
        pickle.dump(content, po)


def pickle_loader(input_file):
    """ Use the Pickle module to de-serialize an object.

        Parameters
        ----------
        input_file : str
            Path to file with byte stream to be de-serialized.

        Returns
        -------
        content : type
            Variable that refers to the de-serialized
            object.
    """

    with open(input_file, 'rb') as pi:
        content = pickle.load(pi)

    return content


def read_tabular(input_file, delimiter='\t'):
    """ Read tabular file.

        Parameters
        ----------
        input_file : str
            Path to a tabular file.
        delimiter : str
            Delimiter used to separate file fields.

        Returns
        -------
        lines : list
            A list with a sublist per line in the input file.
            Each sublist has the fields that were separated by
            the defined delimiter.
    """

    with open(input_file, 'r') as infile:
        reader = csv.reader(infile, delimiter=delimiter)
        lines = [line for line in reader]

    return lines


def import_sequences(fasta_path):
    """ Imports sequences from a FASTA file.

        Parameters
        ----------
        fasta_path : str
            Path to a FASTA file.

        Returns
        -------
        seqs_dict : dict
            Dictionary that has sequences ids as keys and
            sequences as values.
    """

    records = SeqIO.parse(fasta_path, 'fasta')
    seqs_dict = {rec.id: str(rec.seq.upper()) for rec in records}

    return seqs_dict


def flatten_list(list_to_flatten):
    """ Flattens one level of a nested list.

        Parameters
        ----------
        list_to_flatten : list
            List with nested lists.

        Returns
        -------
        flattened_list : str
            Input list flattened by one level.
    """

    flattened_list = list(itertools.chain(*list_to_flatten))

    return flattened_list


def sequence_kmerizer(sequence, k_value, offset=1, position=False):
    """ Decomposes a sequence into kmers.

        Parameters
        ----------
        sequence : str
            Sequence to divide into kmers.
        k_value : int
            Value for the size k of kmers.
        offset : int
            Value to indicate offset of consecutive kmers.
        position : bool
            If the start position of the kmers in the sequence
            should be stored.

        Returns
        -------
        kmers : list
            List with the kmers determined for the input
            sequence. The list will contain strings if
            it is not specified that positions should be
            stored and tuples of kmer and start position
            if the position is stored.
    """

    if position is False:
        kmers = [sequence[i:i+k_value]
                 for i in range(0, len(sequence)-k_value+1, offset)]
    elif position is True:
        kmers = [(sequence[i:i+k_value], i)
                 for i in range(0, len(sequence)-k_value+1, offset)]

    return kmers


def join_list(lst, link):
    """ Joins all elements in a list into a single string.

        Parameters
        ----------
        lst : list
            List with elements to be joined.
        link : str
            Character used to join list elements.

        Returns
        -------
        joined_list : str
            A single string with all elements in the input
            list joined by the character chosen as link.
    """

    joined_list = link.join(lst)

    return joined_list


def write_to_file(text, output_file, write_mode, end_char):
    """ Writes a single string to a file.

        Parameters
        ----------
        text : str
            A single string to write to the output file.
        output_file : str
            Path to the output file.
        write_mode : str
            Write mode can be 'w', writes text and overwrites
            any text in file, or 'a', appends text to text
            already in file.
        end_char : str
            Character added to the end of the file.
    """

    with open(output_file, write_mode) as out:
        out.write(text+end_char)


def write_lines(lines, output_file):
    """ Writes a list of strings to a file. The strings
        are joined with newlines before being written to
        file.

        Parameters
        ----------
        lines : list
            List with the lines/strings to write to the
            output file.
        output_file : str
            Path to the output file.
    """

    joined_lines = join_list(lines, '\n')

    write_to_file(joined_lines, output_file, 'a', '\n')


def concatenate_files(files, output_file, header=None):
    """ Concatenates the contents of a set of files.

        Parameters
        ----------
        files : list
            List with the paths to the files to concatenate.
        output_file : str
            Path to the output file that will store the
            concatenation of input files.
        header : str or NoneType
            Specify a header that should be written as the
            first line in the output file.

        Returns
        -------
        output_file : str
            Path to the output file that was created with
            the concatenation of input files.
    """

    with open(output_file, 'w') as of:
        if header is not None:
            of.write(header)
        for f in files:
            with open(f, 'r') as fd:
                shutil.copyfileobj(fd, of)

    return output_file


def hash_sequence(string):
    """ Compute SHA256 for an input string.

        Parameters
        ----------
        string : str
            Input string to hash.

        Returns
        -------
        sha256 : str
            String representation of the sha256 HASH object.
    """

    sha256 = hashlib.sha256(string.encode('utf-8')).hexdigest()

    return sha256


def fasta_str_record(seqid, sequence):
    """ Creates the string representation of a FASTA record.

        Parameters
        ----------
        seqid : str
            Sequence identifier to include in the header.
        sequence : str
            String representing DNA or Protein sequence.

        Returns
        -------
        record : str
            String representation of the FASTA record.
    """

    record = '>{0}\n{1}'.format(seqid, sequence)

    return record


def determine_distinct(sequences_file, unique_fasta):
    """ Identifies duplicated sequences in a FASTA file.
        Returns a single sequence identifier per distinct
        sequence and saves distinct sequences to a FASTA
        file.

        Parameters
        ----------
        sequences_file : str
            Path to a FASTA file.
        unique_fasta : str
            Path to a FASTA file that will be created to
            store distinct sequences.

        Returns
        -------
        List with following elements:
            total : int
                Total number of times sequences were repeated.
            unique_seqids : list
                List with one sequence identifier per distinct
                sequence. The first identifier observed for a
                distinct sequence is the one stored in the list.
    """

    total = 0
    out_seqs = []
    seqs_dict = {}
    exausted = False
    out_limit = 10000
    seq_generator = SeqIO.parse(sequences_file, 'fasta')
    while exausted is False:
        record = next(seq_generator, None)
        if record is not None:
            # seq object has to be converted to string
            sequence = str(record.seq.upper())
            seqid = record.id
            seq_hash = hash_sequence(sequence)

            # store only the hash for distinct sequences
            if seq_hash not in seqs_dict:
                seqs_dict[seq_hash] = seqid
                recout = fasta_str_record(seqid, sequence)
                out_seqs.append(recout)
            elif seq_hash in seqs_dict:
                total += 1
        else:
            exausted = True

        if len(out_seqs) == out_limit or exausted is True:
            if len(out_seqs) > 0:
                out_seqs = join_list(out_seqs, '\n')
                write_to_file(out_seqs, unique_fasta, 'a', '\n')
                out_seqs = []

    unique_seqids = list(seqs_dict.values())

    return [total, unique_seqids]


def run_minimap2(reference, map_fasta, output_file):
    """ Executes minimap2 to map short sequences
        to a reference genome.

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

    minimap_args = ['minimap2 -I 1G --cs -cx sr {0} {1} > '
                    '{2}'.format(reference, map_fasta, output_file)]

    minimap_proc = subprocess.Popen(minimap_args,
                                    shell=True,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)

    stderr = minimap_proc.stderr.readlines()
    stdout = minimap_proc.stdout.readlines()

    return [stdout, stderr]


def create_mmseqs_db(input_file, output_db):
    """ Executes MMseqs2 to create a database with
        sequences in input FASTA file.

        Parameters
        ----------
        input_file : str
            Path to a FASTA file with DNA sequences.
        output_db : str
            Full path that includes prefix used for
            all database files that are created.

        Returns
        -------
        List with following elements:
            db_stdout : list
                List with the stdout from MMseqs2 in bytes.
            db_stderr : list
                List with the stderr from MMseqs2 in bytes.
    """

    mmseqs_args = ['mmseqs', 'createdb', input_file, output_db]

    mmseqs_proc = subprocess.Popen(mmseqs_args,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)

    db_stderr = mmseqs_proc.stderr.readlines()
    db_stdout = mmseqs_proc.stdout.readlines()

    return [db_stdout, db_stderr]


def cluster_baits(database, cluster_db, temp_directory, threads):
    """ Cluster sequences in a MMseqs2 database.

        Parameters
        ----------
        database : str
            Full path that includes prefix used for
            database files.
        cluster_db : str
            Full path that includes prefix used for
            files with clustering results.
        temp_directory : str
            Path to te temporary directory used to
            store intermediate files.
        threads : int
            Number of threads used in the clustering
            process.

        Returns
        -------
        List with following elements:
            stdout : list
                List with the stdout from MMseqs2 in bytes.
            stderr : list
                List with the stderr from MMseqs2 in bytes.
    """

    # cluster
    mmseqs_args = ['mmseqs', 'cluster', '--cov-mode', '0', '-c',
                   '0.8', '--threads', str(threads), database,
                   cluster_db, temp_directory]

    mmseqs_proc = subprocess.Popen(mmseqs_args,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)

    cluster_stderr = mmseqs_proc.stderr.readlines()
    cluster_stdout = mmseqs_proc.stdout.readlines()

    return [cluster_stdout, cluster_stderr]


def align_clusters(database, cluster_db, align_db, threads):
    """ Aligns sequences in a MMseqs2 database against
        clustering results from MMseqs2.

        Parameters
        ----------
        database : str
            Full path that includes prefix used for
            database files.
        cluster_db : str
            Full path that includes prefix used for
            files with clustering results.
        align_db : str
            Full path that includes prefix used for
            files with alignment results.
        threads : int
            Number of threads used in the clustering
            process.

        Returns
        -------
        List with following elements:
            align_stdout : list
                List with the stdout from MMseqs2 in bytes.
            align_stderr : list
                List with the stderr from MMseqs2 in bytes.
    """

    # align to get identities
    align_args = ['mmseqs', 'align', '--cov-mode', '0', '-c',
                  '0.8', '--threads', str(threads), database,
                  database, cluster_db, align_db, '-a']

    align_proc = subprocess.Popen(align_args,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)

    align_stderr = align_proc.stderr.readlines()
    align_stdout = align_proc.stdout.readlines()

    return [align_stdout, align_stderr]


def convert_alignmentDB(database, align_db, align_out):
    """ Converts MMseqs2 alignment results into tabular format
        to add identities and other alignment information to
        clustering results.

        Parameters
        ----------
        database : str
            Full path that includes prefix used for
            database files.
        align_db : str
            Full path that includes prefix used for
            files with alignment results.
        align_out : str
            Path to the output file with the clustering
            results in tabular format.

        Returns
        -------
        List with following elements:
            convert_stdout : list
                List with the stdout from MMseqs2 in bytes.
            convert_stderr : list
                List with the stderr from MMseqs2 in bytes.
    """

    # convert output
    convert_args = ['mmseqs', 'convertalis', database, database,
                    align_db, align_out]

    convert_proc = subprocess.Popen(convert_args,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)

    convert_stderr = convert_proc.stderr.readlines()
    convert_stdout = convert_proc.stdout.readlines()

    return [convert_stdout, convert_stderr]


def determine_breadth_coverage(intervals, total_bases):
    """ Determines the percentage and total number of covered
        bases according to a set of coverage intervals.

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

    # determine breadth of coverage
    covered_bases = 0
    for k, v in intervals.items():
        for e in v:
            covered_bases += sum([1 for p, c in e[2].items() if c > 0])

    breadth_coverage = covered_bases / total_bases

    return [breadth_coverage, covered_bases]


def determine_small_bait(span, bait_size, start, stop, sequence_length):
    """ Determines baits for regions shorter than bait length.

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

    rest = bait_size - span
    bot = math.floor(rest / 2)
    top = math.ceil(rest / 2)
    if start == 0:
        bait_interval = [start, start + bait_size]
    elif (start - bot) < 0:
        bait_interval = [0, top+(bot-start)]
    elif (stop + top) > sequence_length:
        bait_interval = [start-(bot+(top-(sequence_length-stop))), sequence_length]
    else:
        bait_interval = [start-bot, stop+top]

    return bait_interval


def determine_interval_baits(bait_size, start, stop):
    """ Determines baits for regions with length value
        equal or greater than bait size.

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
        probes : list of list
            List with one sublist per determined bait.
            Each sublist has the start and stop position
            for a bait.
    """

    probes = []
    reach = False
    while reach is False:
        if (start + bait_size) == stop:
            bait_interval = [start, stop]
            reach = True
        elif (start + bait_size) > stop:
            diff = (start + bait_size) - stop
            bot_plus = start - diff
            bait_interval = [bot_plus, stop]
            reach = True
        else:
            bait_interval = [start, start + bait_size]
            start = start + bait_size
        probes.append(bait_interval)

    return probes


def count_contigs(fasta, min_len):
    """ Counts the number of records in a FASTA file.
        Only counts sequences that are longer than a
        specified minimum length.

        Parameters
        ----------
        fasta : str
            Path to a FASTA file.
        min_len : int
            Minimum sequence length. Sequences shorter
            than this value are not counted.

        Returns
        -------
        nr_contigs : int
            Number of sequences in input FASTA file
            (longer than specified minimum length).
    """

    contigs = [str(rec.seq) for rec in SeqIO.parse(fasta, 'fasta')]
    nr_contigs = len(contigs)
    valid_contigs = len([seq for seq in contigs if len(seq) >= min_len])
    total_length = sum([len(seq) for seq in contigs])

    return [nr_contigs, valid_contigs, total_length]


def generate_baits(fasta, output_file, bait_size, bait_offset, min_len):
    """ Generates baits for sequences in a FASTA file.

        Parameters
        ----------
        fasta : str
            Path to a FASTA file.
        output_file : str
            Path to the output FASTA file.
        bait_size : int
            Bait size in bases.
        bait_offset : int
            Position offset between start positions of
            subsequent baits.
        min_len : int
            Sequence minimum length. Baits will not be
            determined for sequences shorter than this
            value.

        Returns
        -------
        None
    """

    sequences = import_sequences(fasta)
    kmers = {cid: sequence_kmerizer(seq, bait_size, offset=bait_offset, position=True)
             for cid, seq in sequences.items() if len(seq) >= min_len}

    baits = [['>{0}_{1}\n{2}'.format(k, e[1], e[0])
              for e in v] for k, v in kmers.items()]
    baits = flatten_list(baits)

    write_lines(baits, output_file)

    return len(baits)


def merge_intervals(intervals):
    """ Merges intersecting intervals.

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
        # current and previous intervals intersect
        if current[0] <= previous[1]:
            # determine top position
            previous[1] = max(previous[1], current[1])
            # merge coverage dictionaries
            previous_cov = previous[2]
            current_cov = current[2]
            for k, v in current_cov.items():
                if k not in previous_cov:
                    previous_cov[k] = v
                else:
                    previous_cov[k] += v
            previous[2] = previous_cov
        # current and previous intervals do not intersect
        else:
            merged.append(deepcopy(current))

    return merged


#intervals = v
#identifier = k
#total_len = len(contigs[k])
def determine_missing_intervals(intervals, identifier, total_len):
    """ Determines sequence intervals that are not covered by any
        probes.

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
                interval that is not covered by probes.
            not_covered : int
                Total number of bases not covered by probes.
    """

    start = 0
    not_covered = 0
    missing_regions = {identifier: []}
    for i in intervals:
        diff = i[0] - start
        if diff > 0:
            missing_regions[identifier].append([start, start+diff])
            not_covered += diff
            start += diff

        # create groups of equal values
        values_groups = [list(j) for i, j in groupby(i[2].values())]
        for g in values_groups:
            if g[0] == 0:
                missing_regions[identifier].append([start, start+len(g)])
                not_covered += len(g)
                start += len(g)
            else:
                start += len(g)

    # add terminal region
    if start != total_len:
        missing_regions[identifier].append([start, total_len])
        not_covered += total_len - start

    return [missing_regions, not_covered]


def cover_intervals(intervals, total_len, bait_size, bait_region):
    """ Determines baits to cover specified sequence regions.

        Parameters
        ----------
        intervals : list
            List of lists. Each sublist has start and stop
            positions for sequence regions with no coverage.
        total_len : int
            Total length of the sequence.
        bait_size : int
            Bait size in bases.
        bait_region : int
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
        span = i[1] - i[0]
        if span >= bait_region:
            if span < bait_size:
                bait_interval = determine_small_bait(span, bait_size,
                                                     i[0], i[1],
                                                     total_len)
                cover_baits.append(bait_interval)
            # will slide and determine baits
            # if in the last iter, uncovered region is very small it will overlap
            # with regions that are already covered and increase depth of coverage
            # pass bait_region as arg to stop determining when region is too small?
            elif span >= bait_size:
                probes = determine_interval_baits(bait_size, i[0], i[1])
                cover_baits.extend(probes)

    return cover_baits


def determine_depth_coverage(intervals, total_len):
    """ Determine depth of coverage for a sequence.

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

    # create dictionary to add coverage for all positions
    positions = list(range(0, total_len))
    positions_depth = {p: 0 for p in positions}
    # increment coverage values based on intervals
    for i in intervals:
        for p, c in i[2].items():
            positions_depth[p] += c

    # determine coverage distribution
    counts = sorted(Counter(positions_depth.values()).most_common(),
                    key=lambda x: x[0])

    return [positions_depth, counts]


def regex_matcher(string, pattern):
    """ Finds substrings that match a regex pattern.

        Parameters
        ----------
        string : str
            Input string.
        pattern : str
            Pattern to match. Patterns should
            include 'r' before string to match
            or characters after backslashes will
            be escaped.

        Returns
        -------
        matches : list
            List with substrings that matched the pattern.
    """

    matches = re.findall(pattern, string)

    return matches


def single_position_coverage(coverage_info, start):
    """ Determine if positions in a subsequence are
        covered based on information in the cs field
        in a PAF file created by minimpa2.

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
        # subsequence part with exact matches
        if m[0] == ':':
            # create dctionary entries with coverage = 1
            new_cov = {i: 1 for i in range(start, start+int(m[1:]))}
            coverage = {**coverage, **new_cov}
            # increment start position
            start = start + int(m[1:])
        # position with substitution
        elif m[0] == '*':
            coverage[start] = 0
            start += 1
        # position with deletion
        elif m[0] == '-':
            # coverage 0 for missing bases
            new_cov = {i: 0 for i in range(start, start+len(m[1:]))}
            coverage = {**coverage, **new_cov}
            start = start + len(m[1:])
        # insertion
        elif m[0] == '+':
            # do not add coverage values for positions because
            # insertion does not exist in reference
            pass

    return coverage


def common_suffixes(strings):
    """
    """

    splitted = [os.path.basename(s).split('.') for s in strings]
    common = set(splitted[0]).intersection(*splitted[1:])
    filtered = [[e for e in s if e not in common] for s in splitted]
    joined = {strings[i]: '.'.join(filtered[i]) for i in range(len(filtered))}

    return joined


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
    coverage_info = {}
    for g in genomes:
        gbasename = os.path.basename(g).split('.fasta')[0]
        paf_path = os.path.join(output_dir, gbasename+'.paf')

        contigs = import_sequences(g)
        total_bases = nr_contigs[g][2]

        minimap_std = run_minimap2(g, unique_baits, paf_path)

        # read PAF file lines
        paf_lines = read_tabular(paf_path)

        # filter out matches below bait_coverage*bait_size
        # length includes gaps due to deletions and might span more than query length
        valid_length = [line
                        for line in paf_lines
                        if int(line[10]) >= (bait_coverage*bait_size)]

        # compute alignment identity
        for i in range(len(valid_length)):
            valid_length[i].append(int(valid_length[i][9]) / int(valid_length[i][10]))

        # filter out alignments below defined identity
        valid_pident = [line for line in valid_length if line[-1] >= bait_identity]

        # match alignment string with regex
        pattern = r':[0-9]+|\*[a-z][a-z]|\+[a-z]+|-[a-z]+'
        for i in range(len(valid_pident)):
            current = valid_pident[i][-2]
            valid_pident[i].append(regex_matcher(current, pattern))

        # get information about positions that match to determine coverage
        for i in range(len(valid_pident)):
            current = valid_pident[i][-1]
            start = int(valid_pident[i][7])
            valid_pident[i].append(single_position_coverage(current, start))

        # identify subsequences that are well covered by baits
        covered_intervals = {}
        for l in valid_pident:
            covered_intervals.setdefault(l[5], []).append([int(l[7]), int(l[8]), l[-1]])

        # sort covered intervals
        covered_intervals_sorted = {k: sorted(v, key=lambda x: x[0])
                                    for k, v in covered_intervals.items()}

        # merge overlapping intervals
        # deepcopy to avoid altering original intervals
        merged_intervals = {k: merge_intervals(v)
                            for k, v in covered_intervals_sorted.items()}

        coverage = determine_breadth_coverage(merged_intervals, total_bases)

        # determine subsequences that are not covered
        missing = [determine_missing_intervals(v, k, len(contigs[k]))
                   for k, v in merged_intervals.items()]

        # add missing regions for contigs that had 0 baits mapped
        not_mapped = [[{c: [[0, len(contigs[c])]]}, len(contigs[c])] for c in contigs if c not in merged_intervals]
        missing.extend(not_mapped)

        missing_regions = {k: v for i in missing for k, v in i[0].items()}
        not_covered = sum([i[1] for i in missing])

        coverage_info[gbasename] = [*coverage, not_covered]

        # create baits for missing regions
        if generate is True:
            missing_baits_intervals = {k: cover_intervals(v, len(contigs[k]), bait_size, bait_region)
                                       for k, v in missing_regions.items()}

            extra_probes = {}
            for k, v in missing_baits_intervals.items():
                extra_probes[k] = ['>{0}_{1}\n{2}'.format(k, e[0], contigs[k][e[0]:e[1]]) for e in v]

            new_baits_lines = [v for k, v in extra_probes.items()]
            new_baits_lines = flatten_list(new_baits_lines)

            write_lines(new_baits_lines, unique_baits)

            coverage_info[gbasename].append(len(new_baits_lines))
            total += len(new_baits_lines)

            print('{0:<30}  {1:^10.4f}  {2:^10}  {3:^10}  '
                  '{4:^7}'.format(short_samples[g], coverage[0],
                                  coverage[1], not_covered,
                                  len(new_baits_lines)))

        if depth is True:
            # determine depth of coverage
            depth_values = {}
            for k, v in merged_intervals.items():
                depth_values[k] = determine_depth_coverage(v, len(contigs[k]))

            total_counts = {}
            for k, v in depth_values.items():
                for i in v[1]:
                    total_counts.setdefault(i[0], []).append(i[1])

            total_counts = {k: sum(v) for k, v in total_counts.items()}

            coverage_info[gbasename].extend([depth_values, total_counts])

        if generate is False:
            print('{0:<30}  {1:^10.4f}  {2:^10}  '
                  '{3:^10}'.format(short_samples[g], coverage[0],
                                   coverage[1], not_covered))

    print('-'*len(header))

    return [coverage_info, total]


def exclude_similar_probes(unique_baits, clustering_dir, cluster_identity,
                           cluster_coverage, bait_size, threads):
    """
    """

    print('Clustering probes...')
    # cluster baits and remove based on similarity threshold
    # create database
    mmseqs_db = os.path.join(clustering_dir, 'mmseqs_db')
    mmseqs_std = create_mmseqs_db(unique_baits, mmseqs_db)

    # output paths
    cluster_db = os.path.join(clustering_dir, 'clusters')
    temp_directory = os.path.join(clustering_dir, 'tmp')
    align_db = os.path.join(clustering_dir, 'alignDB')
    align_out = os.path.join(clustering_dir, 'alignOUT')

    os.mkdir(temp_directory)
    # clustering
    cluster_std = cluster_baits(mmseqs_db, cluster_db,
                                temp_directory, threads)

    # align clusters
    align_std = align_clusters(mmseqs_db, cluster_db, align_db, threads)

    # convert alignments
    convert_std = convert_alignmentDB(mmseqs_db, align_db, align_out)

    # read clustering results
    cluster_lines = read_tabular(align_out)
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
    exclude = flatten_list(exclude)
    excluded_seqids = [e[0] for e in exclude]
    print('Excluded {0} probes highly similar to other '
          'probes.'.format(len(excluded_seqids)))

    # create FASTA without excluded baits
    baits = import_sequences(unique_baits)
    baits = {k: v for k, v in baits.items() if k not in excluded_seqids}

    baits_records = ['>{0}\n{1}'.format(k, v) for k, v in baits.items()]
    filtered_baits = os.path.join(clustering_dir, 'filtered_baits')
    write_lines(baits_records, filtered_baits)

    return [filtered_baits, excluded_seqids]


def exclude_contaminant(unique_baits, exclude_regions, exclude_pident,
                        exclude_coverage, bait_size, output_dir):
    """
    """

    print('Mapping against  and removing similar probes...')
    # map against target genome that baits should not be specific for
    gbasename = os.path.basename(exclude_regions).split('.fna')[0]
    paf_path = os.path.join(output_dir, gbasename+'.paf')
    minimap_std = run_minimap2(exclude_regions, unique_baits, paf_path)

    # import mapping results
    mapped_probes = read_tabular(paf_path)
    multispecific_probes = [l[0] for l in mapped_probes
                            if (int(l[9])/int(l[10])) >= exclude_pident
                            and (int(l[3])-int(l[2])) >= (bait_size*exclude_coverage)]

    # remove probes and write final probe set
    baits = import_sequences(unique_baits)
    baits = {k: v for k, v in baits.items() if k not in multispecific_probes}

    print('Removed {0} probes similar with contaminant '
          'genome.'.format(len(multispecific_probes)))

    baits_records = ['>{0}\n{1}'.format(k, v) for k, v in baits.items()]
    final_baits = os.path.join(output_dir, 'final_baits.fasta')
    write_lines(baits_records, final_baits)

    return [final_baits, multispecific_probes]


def write_depth(identifier, depth_values, output_dir):
    """
    """

    depth_file = os.path.join(output_dir, identifier+'_depth.tsv')
    depth_lines = []
    for k, v in depth_values.items():
        depth_lines.append(k)
        depth_lines.extend(['{0}\t{1}'.format(p, e) for p, e in v[0].items()])

    write_lines(depth_lines, depth_file)

    return depth_file


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


def order_contigs(input_files):
    """
    """

    ordered_contigs = {}
    for g in input_files:
        basename = os.path.basename(g)
        contigs = [[rec.id, str(rec.seq), len(str(rec.seq))]
                   for rec in SeqIO.parse(g, 'fasta')]
        contigs = sorted(contigs, key=lambda x: len(x[1]), reverse=True)
        ordered_contigs[basename.split('.fasta')[0]] = [[c[0], c[2]] for c in contigs]

    return ordered_contigs


input_files = '/home/rfm/Desktop/rfm/Lab_Analyses/pneumo_baits_design/assemblies'
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
# add dropdown to change data that is displayed, it will show only regions uncovered before mapping baits!
# possible to change info shown in table with dropdown as filters?
# option to receive baits as input and add more baits generated based on input genomes
# try to use Scattergl and other gl variants to reduce plot size!
# add graph with baits as nodes to empty scatter plot? (scatterplotgl to deal with many points)
# or add tree? stats for genomes that are not related with coverage?
# dropdown to change data displayed on graph based on identity percentage at several thresholds?
# final table with number of uncovered regions per sample, number of SNPs, number of deletions, insertions, etc.
def main(input_files, output_dir, minlen_contig, contig_boundaries,
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

    # get short identifiers
    short_samples = common_suffixes(genomes)

    # determine number of contigs and total length
    nr_contigs = {f: count_contigs(f, minlen_contig) for f in genomes}

    # select assemblies with lowest number of contigs
    sorted_contigs = sorted(list(nr_contigs.items()), key=lambda x: x[1][1])
    ref_set = [t[0] for t in sorted_contigs[0:number_refs]]

    # shred genomic sequences
    # not generating kmers that cover the end of the sequences!
    baits_file = os.path.join(output_dir, 'baits.fasta')
    for g in ref_set:
        nr_baits = generate_baits(g, baits_file, bait_size,
                                  bait_offset, minlen_contig)

    print('\nCreated initial set of {0} probes based on {1} '
          'assemblies.'.format(nr_baits, number_refs))

    # identify unique baits
    unique_baits = os.path.join(output_dir, 'unique_baits.fasta')
    total, unique_seqids = determine_distinct(baits_file, unique_baits)
    print('Removed {0} repeated probes.\n'.format(total))

    # start mapping baits against remaining genomes
    # mapping against ref_set to cover missing regions
    coverage_info = incremental_bait_generator(genomes, unique_baits,
                                               output_dir, bait_size,
                                               bait_coverage, bait_identity,
                                               bait_region, nr_contigs,
                                               short_samples,
                                               generate=True, depth=False)

    print('Added {0} probes to cover {1} assemblies.\nTotal '
          'of {2} probes.'.format(coverage_info[1], len(genomes), nr_baits+coverage_info[1]))

    if cluster_probes is True:
        clustering_dir = os.path.join(output_dir, 'clustering')
        os.mkdir(clustering_dir)
        unique_baits, removed = exclude_similar_probes(unique_baits, clustering_dir,
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
    final_info = incremental_bait_generator(genomes, unique_baits, output_dir,
                                            bait_size, bait_coverage,
                                            bait_identity, bait_region,
                                            nr_contigs, short_samples,
                                            generate=False, depth=True)

    # save depth values
    depth_files_dir = os.path.join(output_dir, 'depth_files')
    os.mkdir(depth_files_dir)
    depth_files = [write_depth(k, v[3], depth_files_dir) for k, v in final_info[0].items()]

    if report is True:
        # determine contig order from longest to shortest
        ordered_contigs = order_contigs(genomes)

        # create plots
        report_dir = os.path.join(output_dir, 'plots')
        os.mkdir(report_dir)

        # determine number of exclude regions and total bps
        exclude_stats = [len(rec) for rec in SeqIO.parse(exclude_regions, 'fasta')]
        total_bps = sum(exclude_stats)

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
                   'Exclude regions': '{0} regions ({1}bps)'.format(len(exclude_stats), total_bps),
                   'Exclude identity': exclude_pident,
                   'Exclude coverage': exclude_coverage,
                   'Create report': str(report)}

        ref_ids = [os.path.basename(f).split('.fasta')[0] for f in ref_set]
        create_report(coverage_info, final_info, report_dir, short_samples,
                      ordered_contigs, fixed_xaxis, fixed_yaxis, ref_ids,
                      nr_contigs, configs)

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
