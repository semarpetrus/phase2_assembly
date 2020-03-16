#!/usr/bin/env python3
from os.path import isfile, splitext
from Bio import SeqIO
import gzip
import argparse
import sys


def main():
    parser = argparse.ArgumentParser(
        add_help=True,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        '--spades_ids', metavar="File", type=str, action='store',
        dest='spades_ids', required=True,
        help=("List of spades ids to ignore."))
    parser.add_argument(
        '--ONT_ids', metavar="File", type=str, action='store',
        dest='ONT_ids', required=True,
        help=("List of Oxford contigs to keep."))
    parser.add_argument(
        '--spades_file', metavar="File", type=str, action='store',
        dest='spades_file', required=True,
        help=("Spades Assembly in fasta format."))
    parser.add_argument(
        '--ONT_file', metavar="File", type=str, action='store',
        dest='ONT_file', required=True,
        help=("Oxford reads file in fastq format."))
    parser.add_argument(
        '--spades_output', metavar="File", type=str, action='store',
        dest='spades_output', required=True, help=(
            "Output file in format of fasta."
            " This contains the spades contigs not found in ONT reads."))
    parser.add_argument(
        '--ont_output', metavar="File", type=str, action='store',
        dest='ont_output', required=True, help=(
            "Output file in format of fasta."
            " This contains the ONT reads that contain the spades reads."))

    args = parser.parse_args()
    spades_ids_file = args.spades_ids
    ONT_ids_file = args.ONT_ids
    spades_file = args.spades_file
    ONT_file = args.ONT_file
    spades_out = args.spades_output
    ONT_out = args.ont_output
    # Check that input files exist
    for file in [spades_ids_file, ONT_ids_file, spades_file, ONT_file]:
        if not isfile(file):
            print(f"{file} not found!", file=sys.stderr)
            return 1
    # Read ids files
    spades_ids = set(open(spades_ids_file, "r").read().splitlines())
    ONT_ids = set(open(ONT_ids_file, "r").read().splitlines())
    # Open spades output file
    with open(spades_out, "wt") as output_handle:
        # Write spades contigs
        file_prefix, file_suffix = splitext(spades_file)
        if file_suffix == '.gz':
            with gzip.open(spades_file, 'rt') as input_ob:
                for record in SeqIO.parse(input_ob, 'fasta'):
                    if record.id not in spades_ids:
                        SeqIO.write(record, output_handle, "fasta")
        else:
            with open(spades_file, 'r') as input_ob:
                for record in SeqIO.parse(input_ob, 'fasta'):
                    if record.id not in spades_ids:
                        SeqIO.write(record, output_handle, "fasta")
        file_prefix, file_suffix = splitext(spades_file)
    # Open ONT output file
    with open(ONT_out, "wt") as output_handle:
        # Write spades contigs
        file_prefix, file_suffix = splitext(ONT_file)
        if file_suffix == '.gz':
            with gzip.open(ONT_file, 'rt') as input_ob:
                for record in SeqIO.parse(input_ob, 'fastq'):
                    if record.id in ONT_ids:
                        SeqIO.write(record, output_handle, "fasta")
        else:
            with open(ONT_file, 'r') as input_ob:
                for record in SeqIO.parse(input_ob, 'fastq'):
                    if record.id in ONT_ids:
                        SeqIO.write(record, output_handle, "fasta")
        file_prefix, file_suffix = splitext(spades_file)


if __name__ == "__main__":
    main()
