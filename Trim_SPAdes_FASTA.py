#!/usr/bin/env python 
# Created by: Lee Bergstrand
# Modified by: Matt McInnes
# Descript: Trims FASTAs from the assembler SPAdes by coverage to remove low coverage contigs.
#
# Usage: Trim_SPAdes_FASTA_By_Cov.py <sequences.fasta> <length threshold> <coverage threshold>
# Example: Trim_SPAdes_FASTA_By_Cov.py mySeqs.fasta 1000 20
# ----------------------------------------------------------------------------------------
# ===========================================================================================================
# Imports:

import argparse
import functools
import multiprocessing
import os
import re
import sys


class FASTAFormat(object):
    def __init__(self, node, length, coverage, identifier, component, sequence):
        self.node = node
        self.length = length
        self.coverage = coverage
        self.identifier = identifier
        self.component = component
        self.sequence = sequence

    def __repr__(self):
        marker = '\n============================================================\n'
        outstr = marker + "Node: %d \nLength: %d \nCoverage: %f \nIdentifier: %d \nComponent: %d \nSequence:\n\n%s" % (
            self.node, self.length, self.coverage, self.identifier, self.component, self.sequence)

        return outstr

    def __str__(self):
        return self.__repr__()


# ---------------------------------
# Command line interface controller
# ---------------------------------
def main(args):
    if args.input is None or args.input == []:
        parser.print_help()
        print("\nError: At least one input file is needed")
        sys.exit(1)
    else:
        input_paths = args.input

    if args.length is None or args.length == []:
        length_cut = None
        is_length_cut = False
    else:
        length_cut = int(args.length[0])
        is_length_cut = True

    if args.coverage is None or args.coverage == []:
        coverage_cut = None
        is_coverage_cut = False
    else:
        coverage_cut = float(args.coverage[0])
        is_coverage_cut = True

    if not is_length_cut and not is_coverage_cut:
        parser.print_help()
        print("\nError: A length or coverage cutoff must be specified")
        sys.exit(1)

    if args.processes is None or args.processes == []:
        try:
            process_count = multiprocessing.cpu_count()
        except NotImplementedError:
            process_count = 2
    else:
        process_count = int(args.processes[0])

    out_dir = args.outdir[0] if args.outdir else os.getcwd()

    # Generates multiprocessing pool for parallel processing of input files.
    pool = multiprocessing.Pool(processes=process_count)
    pool.map(functools.partial(filter_fasta, out_dir=out_dir, length_cut=length_cut, coverage_cut=coverage_cut),
             input_paths)

    print(">> Done")
    sys.exit(0)


# ----------------------------------------------
# Filters each FASTA file by length and coverage
# ----------------------------------------------
def filter_fasta(input_path, out_dir, length_cut, coverage_cut):
    in_path, in_file = os.path.split(input_path)
    out_path = os.path.join(out_dir, in_file) + ".out"

    print(">> Filtering " + in_file)

    seq_objects = generate_sequence_objects(input_path)
    seq_objects_filtered = [obj for obj in seq_objects if not filter_seq_object(obj, length_cut, coverage_cut)]

    fasta = generate_fasta_string(seq_objects_filtered)

    with open(out_path, "w") as new_fasta:
        new_fasta.write(fasta)

    print(">> Finished Filtering " + in_file)


# ----------------------------------
# Generates list of sequence objects
# ----------------------------------
def generate_sequence_objects(input_file):
    seq_objects = []

    with open(input_file, 'rU') as fasta_file:
        fasta_data = fasta_file.read()
        sequence_blocks = re.findall(">[^>]+", fasta_data)

        for block in sequence_blocks:
            header = None
            sequence = None

            try:
                header, sequence = block.split("\n", 1)  # Split each FASTA into header and sequence.
            except ValueError:
                print("Error: Failed to split header and sequence")

            if header and sequence:
                fasta_obj = create_fasta(header, sequence)

                if fasta_obj:
                    seq_objects.append(fasta_obj)
                else:
                    print("Error: Missing FASTA header value")
                    continue
            else:
                print("Error: FASTA sequence block is missing a header or sequence")
                continue

    return seq_objects


# ---------------------------------------------------
# Creates fasta object from input header and sequence
# ---------------------------------------------------
def create_fasta(header, sequence):
    header_split = header.lstrip('>').split('_')
    header_dict = {}

    for i in range(0, len(header_split), 2):
        header_dict[str(header_split[i]).lower()] = header_split[i + 1]

    sequence_node, sequence_length, sequence_coverage = (None, None, None)
    try:
        sequence_node = int(header_dict['node'])
        sequence_length = int(header_dict['length'])
        sequence_coverage = float(header_dict['cov'])
    except KeyError:
        print("Error: FASTA header is missing node, length or coverage")

    if header_dict.get('id'):
        sequence_id = int(header_dict.get('id'))
    else:
        sequence_id = None

    if header_dict.get('component'):
        component = int(header_dict.get('component'))
    else:
        component = None

    valid_header = True
    if None in (sequence_node, sequence_length, sequence_coverage):
        valid_header = False

    if sequence_id is None and component is None:
        valid_header = False

    if valid_header:
        return FASTAFormat(sequence_node, sequence_length, sequence_coverage, sequence_id, component, sequence)
    else:
        return False


# ----------------------------------------------------
# Filters each sequence objects by length and coverage
# ----------------------------------------------------
def filter_seq_object(seq_object, length_cut, coverage_cut):
    remove = False
    length = seq_object.length
    coverage = seq_object.coverage

    if length_cut is not None and (length < length_cut):
        remove = True
    elif coverage_cut is not None and (coverage < coverage_cut):
        remove = True

    if remove is True:
        if seq_object.identifier:
            print(">> Removing sequence with ID: %d" % seq_object.identifier)
        else:
            print(">> Removing sequence with NODE: %d and COMPONENT: %d" % (seq_object.node, seq_object.component))

    return remove


# ----------------------------------------
# Generates final FASTA string for writing
# ----------------------------------------
def generate_fasta_string(seq_objects):
    output_fasta = []
    for seq in seq_objects:
        out_fasta_header= ">NODE_%d_length_%d_cov_%f" % (
            seq.node, seq.length, seq.coverage)

        if seq.identifier:
            out_fasta_header += "_ID_%d\n" % seq.identifier
        else:
            out_fasta_header += "_component_%d\n" % seq.component

        sequence = seq.sequence

        output_fasta.append(out_fasta_header + sequence)

    output_fasta = ''.join(output_fasta)

    return output_fasta


if __name__ == '__main__':
    # -------------------------------
    # Command line interface options.
    # -------------------------------
    parser = argparse.ArgumentParser(description='Trims FASTAs from the assembler SPAdes by length and coverage')
    parser.add_argument('-i', '--input', metavar='INPATH', nargs='+', help='''Input file paths for the output files''')
    parser.add_argument('-o', '--outdir', metavar='OUTPATH', nargs=1, help='''The directory for the output files''')
    parser.add_argument('-l', '--length', metavar='LENGTH', nargs=1, help='''The length cutoff value''')
    parser.add_argument('-c', '--coverage', metavar='COVERAGE', nargs=1, help='''The coverage cutoff value''')
    parser.add_argument('-p', '--processes', metavar='PROCESSES', nargs=1, help='''The coverage cutoff value''')

    cli_args = parser.parse_args()
    main(cli_args)
