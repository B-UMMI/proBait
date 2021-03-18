#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""


import os
import re
import csv
import pickle
import shutil
import hashlib
import itertools

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


def order_contigs(input_files, short_ids):
    """
    """

    ordered_contigs = {}
    for g in input_files:
        contigs = [[rec.id, str(rec.seq), len(str(rec.seq))]
                   for rec in SeqIO.parse(g, 'fasta')]
        contigs = sorted(contigs, key=lambda x: len(x[1]), reverse=True)
        ordered_contigs[short_ids[g]] = [[c[0], c[2]] for c in contigs]

    return ordered_contigs


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


def common_suffixes(strings):
    """
    """

    splitted = [os.path.basename(s).split('.') for s in strings]
    common = set(splitted[0]).intersection(*splitted[1:])
    filtered = [[e for e in s if e not in common] for s in splitted]
    joined = {strings[i]: '.'.join(filtered[i]) for i in range(len(filtered))}

    return joined


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


def get_baits_pos(baits_fasta, short_ids):
    """
    """

    baits_records = SeqIO.parse(baits_fasta, 'fasta')
    baits_pos = {s: {} for s in short_ids.values()}
    for rec in baits_records:
        genome = (rec.id).split('_')[0]
        pos = (rec.id).split('_')[-1]
        contig = (rec.id).split('_{0}'.format(pos))[0]
        baits_pos[genome].setdefault(contig, []).append(pos)

    return baits_pos


def reverse_str(string):
    """ Reverse character order in input string.

        Parameters
        ----------
        string : str
         String to be reversed.

        Returns
        -------
        revstr : str
            Reverse of input string.
    """

    revstr = string[::-1]

    return revstr


def reverse_complement(dna_sequence):
    """ Determines the reverse complement of given DNA strand.

        Parameters
        ----------
        dna_sequence : str
            String representing a DNA sequence.

        Returns
        -------
        reverse_complement_strand : str
            The reverse complement of the DNA sequence (lowercase
            is converted to uppercase).
    """

    base_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                       'a': 'T', 'c': 'G', 'g': 'C', 't': 'A',
                       'n': 'N', 'N': 'N'}

    # convert string into list with each character as a separate element
    bases = list(dna_sequence)

    # determine complement strand
    # default to 'N' if nucleotide is not in base_complement dictionary
    bases = [base_complement.get(base, 'N') for base in bases]

    complement_strand = ''.join(bases)

    # reverse strand
    reverse_complement_strand = reverse_str(complement_strand)

    return reverse_complement_strand
