#!/usr/bin/env python3
from os.path import isfile, splitext
import sys
from Bio import SeqIO
import gzip
import argparse


def main():
    parser = argparse.ArgumentParser(
        add_help=True,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        '--sample', metavar="String", type=str, action='store',
        dest='sample', required=True, help=("Sample name."))
    parser.add_argument(
        '--unmapped', metavar="File", type=str, action='store',
        dest='unmapped', required=True, help=("Unmapped spades contigs."))
    parser.add_argument(
        '--context', metavar="File", type=str, action='store',
        dest='context', required=True,
        help=("ONT context contigs."))
    parser.add_argument(
        '--flye', metavar="File", type=str, action='store',
        dest='flye', required=True,
        help=("Flye assembly."))
    parser.add_argument(
        '--output', metavar="File", type=str, action='store',
        dest='output', required=True, help=(
            "Output file in format of fasta. "
            "This contains all contigs in input files "
            "with reformatted names."))
    args = parser.parse_args()
    sample = args.sample
    unmapped = args.unmapped
    context = args.context
    flye = args.flye
    output = args.output
    for file in [unmapped, context, flye]:
        if not isfile(file):
            print(f"{file} not found!", file=sys.stderr)
            return 1
    # Open output file
    with open(output, "wt") as output_handle:
        # Loop Through each input file
        # Write sequences with reformatted headers
        for source, input_file in zip(
            ['spades', 'recovered_reads', 'flye'],
                [unmapped, context, flye]):
            file_prefix, file_suffix = splitext(input_file)
            if file_suffix == '.gz':
                with gzip.open(input_file, 'rt') as input_ob:
                    for count, record in enumerate(
                            SeqIO.parse(input_ob, 'fasta')):
                        record.id = (
                            f'{sample}_{source}_{count}_{len(record.seq)}')
                        record.description = (
                            f'{sample}_{source}_{count}_{len(record.seq)}')
                        SeqIO.write(record, output_handle, "fasta")
            else:
                with open(input_file, 'r') as input_ob:
                    for count, record in enumerate(
                            SeqIO.parse(input_ob, 'fasta')):
                        record.id = (
                            f'{sample}_{source}_{count}_{len(record.seq)}')
                        record.description = (
                            f'{sample}_{source}_{count}_{len(record.seq)}')
                        SeqIO.write(record, output_handle, "fasta")


if __name__ == "__main__":
    main()
