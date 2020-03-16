#!/usr/bin/env python3
from os.path import isfile, splitext
import sys
from Bio import SeqIO
import gzip
import pandas as pd
from collections import deque
import argparse


def main():
    parser = argparse.ArgumentParser(
        add_help=True,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        '-b', '--blast', metavar="File", type=str, action='store',
        dest='blast', required=True,
        help=("Blast xml file in format of '6 std qlen slen'."))
    parser.add_argument(
        '-q', '--query', metavar="File", type=str, action='store',
        dest='query', required=True, help=(
            "Query file in format of fasta. "
            "This will the file that will be reduced."))
    parser.add_argument(
        '-o', '--output', metavar="File", type=str, action='store',
        dest='output', required=True, help=(
            "Output file in format of fasta."
            " This contains the contigs in the query file "
            "not covered by the blast results."))
    args = parser.parse_args()
    blast_file = args.blast
    query_file = args.query
    output_file = args.output
    # Check if input files exist
    for file in [query_file, blast_file]:
        if not isfile(file):
            print(f"{file} not found!", file=sys.stderr)
            return 1
    # Read blast file into dataframe
    blast_df = pd.read_csv(blast_file, sep='\t', header=None)
    blast_df.columns = [
        'qid', 'sid', 'pident', 'len', 'mismatch',
        'gap', 'qstart', 'qend', 'sstart', 'send',
        'evalue', 'bitscore', 'qlen', 'slen']
    # Filter blast results based on percent identity >= 90 percs
    blast_df = blast_df[blast_df['pident'] >= 90]
    # Adjust blast alignment for all to be in forward strand of query
    blast_df['qsn_end'] = blast_df[['qstart', 'qend']].apply(
        lambda x: max(x), axis=1)
    blast_df['qsn_start'] = blast_df[['qstart', 'qend']].apply(
        lambda x: min(x), axis=1)
    # Group dataframe based on query ids
    qids_group = blast_df.groupby(['qid'])
    patched_deques = deque([])
    # Get all unique query ids
    total_keys = len(qids_group.groups.keys())
    # Loop through each query id and patch the blast hits
    for count, query in enumerate(qids_group.groups.keys(), start=1):
        print(f'{count}/{total_keys}', end='\r')
        sids_group = qids_group.get_group(query).groupby(['sid'])
        for subject in sids_group.groups.keys():
            sid_group = sids_group.get_group(subject).sort_values(
                ['qsn_start', 'qsn_end'])
            patched_deque = deque()
            for qid, sid, start, end, qlen in sid_group[[
                    'qid', 'sid', 'qsn_start', 'qsn_end', 'qlen']].to_numpy():
                if not patched_deque:
                    patched_deque.append([qid, sid, start, end, qlen])
                    continue
                qid, sid, lstart, lend, qlen = patched_deque.pop()
                if start > lend+1:
                    patched_deque.append([qid, sid, lstart, lend, qlen])
                    patched_deque.append([qid, sid, start, end, qlen])
                else:
                    patched_deque.append(
                        [qid, sid, lstart, max(end, lend), qlen])
            patched_deques.extend(patched_deque)
    patched_df = pd.DataFrame(
        patched_deques, columns=['qid', 'sid', 'start', 'end', 'qlen'])
    # Get length of coverage
    patched_df['len'] = patched_df.apply(
        lambda x: x[3]-x[2]+1, axis=1)
    sumdf = patched_df.groupby(['qid', 'sid'])['len'].sum().reset_index()
    spades_dict = patched_df.set_index(['qid'])['qlen'].to_dict()
    sumdf['perc'] = sumdf[['qid', 'len']].apply(
        lambda x: x[1]/spades_dict[x[0]], axis=1)
    # Select contigs with percentage coverage greater than 90
    covered_df = sumdf.loc[sumdf['perc'] >= 0.9]
    covered_contigs = set(covered_df['qid'].unique())
    with open(output_file, "wt") as output_handle:
        file_prefix, file_suffix = splitext(query_file)
        if file_suffix == '.gz':
            with gzip.open(query_file, 'rt') as input_ob:
                for record in SeqIO.parse(input_ob, 'fasta'):
                    if record.id not in covered_contigs:
                        SeqIO.write(record, output_handle, "fasta")
        else:
            with open(query_file, 'r') as input_ob:
                for record in SeqIO.parse(input_ob, 'fasta'):
                    if record.id not in covered_contigs:
                        SeqIO.write(record, output_handle, "fasta")


if __name__ == "__main__":
    main()
