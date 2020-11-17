#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 22:26:18 2020

@author: rfm
"""


import os
import csv
import math
import shutil
import pickle
import hashlib
import argparse
import itertools
import subprocess
from copy import deepcopy

from Bio import SeqIO


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
    seqs_dict = {}
    out_limit = 10000
    out_seqs = []
    exausted = False
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
    """
    """

    minimap_args = ['minimap2 -I 1G --cs -cx sr {0} {1} > '
                    '{2}'.format(reference, map_fasta, output_file)]

    minimap_proc = subprocess.Popen(minimap_args,
                                    shell=True,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)

    stderr = minimap_proc.stderr.readlines()

    return stderr


def create_mmseqs_db(input_file, output_db):
    """
    """

    mmseqs_args = ['mmseqs', 'createdb', input_file, output_db]

    mmseqs_proc = subprocess.Popen(mmseqs_args,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)

    stderr = mmseqs_proc.stderr.readlines()

    return stderr


def cluster_baits(database, cluster_db, temp_directory):
    """
    """

    # cluster
    mmseqs_args = ['mmseqs', 'cluster', '--cov-mode', '0', '-c', '0.8',
                   '--threads', '4', database, cluster_db, temp_directory]

    mmseqs_proc = subprocess.Popen(mmseqs_args,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)

    cluster_stderr = mmseqs_proc.stderr.readlines()
    cluster_stdout = mmseqs_proc.stdout.readlines()

    return [cluster_stdout, cluster_stderr]


def align_clusters(database, cluster_db, align_db):
    """
    """

    # align to get identities
    align_args = ['mmseqs', 'align', '--cov-mode', '0', '-c', '0.8',
                  '--threads', '4', database, database, cluster_db,
                  align_db, '-a']

    align_proc = subprocess.Popen(align_args,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)

    align_stderr = align_proc.stderr.readlines()
    align_stdout = align_proc.stdout.readlines()

    return [align_stdout, align_stderr]


def convert_alignmentDB(database, align_db, align_out):
    """
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


input_files = '/home/rfm/Desktop/rfm/Lab_Analyses/pneumo_baits_design/assemblies'
output_dir = '/home/rfm/Desktop/rfm/Lab_Analyses/pneumo_baits_design'
bait_size = 120
bait_offset = 120
number_refs = 1
bait_identity = 1.0
cluster_identity = 0.99
minlen_contig = bait_size * 2
contaminant_genome = '/home/rfm/Desktop/rfm/Lab_Analyses/pneumo_baits_design/ncbi-genomes-2020-11-16/GCF_000001405.39_GRCh38.p13_genomic.fna'


def main(input_files, bait_size, output_dir, number_refs, bait_identity,
         cluster_identity, minlen_contig, contaminant_genome):

    genomes = [os.path.join(input_files, file)
               for file in os.listdir(input_files)]

    # select genomes to shred
    nr_contigs = []
    for g in genomes:
        contigs = [rec for rec in SeqIO.parse(g, 'fasta')
                   if len(str(rec.seq)) >= minlen_contig]
        nr_contigs.append((g, len(contigs)))

    # select assemblies with lowest number of contigs
    nr_contigs = sorted(nr_contigs, key=lambda x: x[1])

    ref_set = [t[0] for t in nr_contigs[0:number_refs]]
    genomes = list(set(genomes) - set(ref_set))
    genomes.sort()
    # shred genomic sequences
    # not generating kmers that cover the end of the sequences!
    baits_files = []
    for g in ref_set:
        sequences = import_sequences(g)
        kmers = {cid: sequence_kmerizer(seq, bait_size, offset=bait_offset, position=True)
                 for cid, seq in sequences.items() if len(seq) >= bait_size}

        baits = [['>{0}_{1}\n{2}'.format(k, e[1], e[0])
                  for e in v] for k, v in kmers.items()]
        baits = flatten_list(baits)

        gbasename = os.path.basename(g).split('.fasta')[0]
        output_file = os.path.join(output_dir, gbasename+'_baits.fasta')
        write_lines(baits, output_file)

        baits_files.append(output_file)

    # concatenate files with baits
    baits_file = os.path.join(output_dir, 'baits.fasta')
    concatenate_files(baits_files, baits_file)

    # identify unique baits
    unique_baits = os.path.join(output_dir, 'unique_baits.fasta')
    total, unique_seqids = determine_distinct(baits_file, unique_baits)

    # start mapping baits against remaining genomes
    for g in genomes:
        gbasename = os.path.basename(g).split('.fasta')[0]
        paf_path = os.path.join(output_dir, gbasename+'.paf')

        stderr = run_minimap2(g, unique_baits, paf_path)

        # read PAF file lines
        paf_lines = read_tabular(paf_path)

        # filter out matches below bait length
        valid_length = [line for line in paf_lines if int(line[10]) == 120]

        # compute alignment identity
        for i in range(len(valid_length)):
            valid_length[i].append(int(valid_length[i][9]) / int(valid_length[i][10]))

        # filter out alignments below defined identity
        valid_length = [line for line in valid_length if line[-1] >= bait_identity]

        # identify subsequences that are well covered by baits
        contigs = {}
        for l in valid_length:
            contigs.setdefault(l[5], []).append([int(l[7]), int(l[8])])

        # sort covered intervals
        contigs = {k: sorted(v, key=lambda x: x[0]) for k, v in contigs.items()}

        # merge overlapping intervals
        # deepcopy to avoid altering original intervals
        merged_intervals = {}
        for k, v in contigs.items():
            merged = [deepcopy(v[0])]
            for current in v:
                previous = merged[-1]
                if current[0] <= previous[1]:
                    previous[1] = max(previous[1], current[1])
                else:
                    merged.append(deepcopy(current))

            merged_intervals[k] = merged

        # determine breath of coverage
        covered_bases = 0
        for k, v in merged_intervals.items():
            for e in v:
                covered_bases += e[1] - e[0]

        # determine subsequences that are not covered
        missing_regions = {}
        not_covered = 0
        sequences = {rec.id: str(rec.seq) for rec in SeqIO.parse(g, 'fasta')}
        for k, v in merged_intervals.items():
            total_len = len(sequences[k])
            start = 0
            missing_regions[k] = []
            for i in v:
                missing_regions[k].append([start, i[0]])
                not_covered += i[0] - start
                start = i[1]
            # add terminal region
            missing_regions[k].append([start, total_len])
            not_covered += total_len - start

        # create baits for missing regions
        missing_baits_intervals = {}
        for k, v in missing_regions.items():
            total_len = len(sequences[k])
            missing_baits_intervals[k] = []
            for e in v:
                span = e[1] - e[0]
                rest = bait_size - span
                bot = math.floor(rest / 2)
                top = math.ceil(rest / 2)
                if span < bait_size:
                    if e[0] == 0:
                        bait_interval = [e[0], e[0] + bait_size]
                    elif (e[0] - bot) < 0:
                        bait_interval = [0, top+(bot-e[0])]
                    elif (e[1] + top) > total_len:
                        bait_interval = [e[0]-(bot+(top-(total_len-e[1]))), total_len]
                    else:
                        bait_interval = [e[0]-bot, e[1]+top]
                    missing_baits_intervals[k].append(bait_interval)
                elif span >= bait_size:
                    new_baits = []
                    reach = False
                    bot2 = e[0]
                    top2 = e[1]
                    while reach is False:
                        if (bot2 + bait_size) == top2:
                            bait_interval = [bot2, top2]
                            reach = True
                        elif (bot2 + bait_size) > top2:
                            diff = (bot2 + bait_size) - top2
                            bot_plus = bot2 - diff
                            bait_interval = [bot_plus, top2]
                            reach = True
                        else:
                            bait_interval = [bot2, bot2+bait_size]
                            bot2 = bot2+bait_size
                        new_baits.append(bait_interval)
                    missing_baits_intervals[k].extend(new_baits)

        extra_baits = {}
        for k, v in missing_baits_intervals.items():
            extra_baits[k] = ['>{0}_{1}\n{2}'.format(k, e[0], sequences[k][e[0]:e[1]]) for e in v]

        new_baits_lines = [v for k, v in extra_baits.items()]
        new_baits_lines = flatten_list(new_baits_lines)

        write_lines(new_baits_lines, unique_baits)
        print('Added {0} baits from genome {1}'.format(len(new_baits_lines), g))

    # cluster baits and remove based on similarity threshold
    # create database
    mmseqs_db = os.path.join(output_dir, 'mmseqs_db')
    stderr = create_mmseqs_db(unique_baits, mmseqs_db)

    # output paths
    cluster_db = os.path.join(output_dir, 'clusters')
    temp_directory = os.path.join(output_dir, 'tmp')
    align_db = os.path.join(output_dir, 'alignDB')
    align_out = os.path.join(output_dir, 'alignOUT')

    os.mkdir(temp_directory)
    # clustering
    cluster_stderr = cluster_baits(mmseqs_db, cluster_db, temp_directory)
    # align clusters
    align_stderr = align_clusters(mmseqs_db, cluster_db, align_db)
    # convert alignments
    convert_stderr = convert_alignmentDB(mmseqs_db, align_db, align_out)

    # read clustering results
    cluster_lines = read_tabular(align_out)
    clusters = {}
    # pident at index 2
    for l in cluster_lines:
        clusters.setdefault(l[0], []).append([l[1], l[2]])

    # exclude clusters with only the representative
    clusters = {k: v for k, v in clusters.items() if len(v) > 1}

    # remove representatives from clusters
    clusters = {k: [e for e in v if e[0] != k]
                for k, v in clusters.items()}

    # get identifiers of baits with identity above threshold
    exclude = [[e for e in v if e[0] != k and float(e[1]) >= cluster_identity]
               for k, v in clusters.items()]
    exclude = flatten_list(exclude)
    excluded_seqids = [e[0] for e in exclude]

    # create FASTA without excluded baits
    baits = import_sequences(unique_baits)
    baits = {k: v for k, v in baits.items() if k not in excluded_seqids}

    baits_records = ['>{0}\n{1}'.format(k, v) for k, v in baits.items()]
    filtered_baits = os.path.join(output_dir, 'filtered_baits')
    write_lines(baits_records, filtered_baits)

    # map against target genome that baits should not be specific for
    gbasename = os.path.basename(contaminant_genome).split('.fna')[0]
    paf_path = os.path.join(output_dir, gbasename+'.paf')
    stderr = run_minimap2(contaminant_genome, filtered_baits, paf_path)


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', type=str, required=True,
                        dest='input_files',
                        help='')

    parser.add_argument('-o', type=str, required=True,
                        dest='output_dir',
                        help='')

    parser.add_argument('-bt', type=int, required=False,
                        default=120,
                        dest='bait_size',
                        help='')

    parser.add_argument('-bo', type=int, required=False,
                        default=60,
                        dest='bait_offset',
                        help='')

    parser.add_argument('--nr', type=int, required=False,
                        default=2,
                        dest='number_refs',
                        help='')

    parser.add_argument('--bi', type=float, required=False,
                        default=0.98,
                        dest='bait_identity',
                        help='')

    parser.add_argument('--ci', type=float, required=False,
                        default=0.98,
                        dest='cluster_identity',
                        help='')
    
    parser.add_argument('--mc', type=int, required=False,
                        default=None,
                        dest='minlen_contig',
                        help='')

    parser.add_argument('--cg', type=str, required=False,
                        default=None,
                        dest='contaminant_genome',
                        help='')

    args = parser.parse_args()

    args.minlen_contig = (args.minlen_contig
                          if args.minlen_contig is not None
                          else args.bait_size)

    return [args.input_files, args.output_dir, args.bait_size,
            args.bait_overlap, args.number_refs, args.bait_identity,
            args.cluster_identity, args.minlen_contig, args.contaminant_genome]


if __name__ == '__main__':

    args = parse_arguments()

    main(*args)
