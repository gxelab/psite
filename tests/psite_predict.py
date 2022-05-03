from operator import index
import sys
import argparse
import re
import gzip
import pickle
from itertools import chain, islice

import pandas as pd
import numpy as np
import pysam
from Bio import SeqIO


# helper functions
def smart_open(path):
    """open plain text or gzipped file"""
    if path[-2:] == 'gz':
        return gzip.open(path, 'rt')
    else:
        return open(path, 'rt')


def strip_version(tx_name):
    """Strip transcript version"""
    return re.sub('\.\d+$', '', tx_name)


def read_fasta(path, ignore_version=False):
    """Construct a dict of input fasta"""
    fa = dict()
    with smart_open(path) as f:
        for record in SeqIO.parse(f, 'fasta'):
            if ignore_version:
                record.id = strip_version(record.id)
            fa[record.id] = str(record.seq)
    return fa


WC_PAIRING = str.maketrans('ACGTN', 'TGCAN')
def rev_comp(s):
    """reverse complementation DNA string"""
    return s[::-1].translate(WC_PAIRING)


def chunk_iter(iterable, size=1024):
    """
    iterate by fixed block size
    
    Yield an iterator for each block. If the last chunk do not have enough
    items, all the remaining items is return as the last chunk.
    reference: https://stackoverflow.com/a/8998040/3926543
    """
    it = iter(iterable)
    while True:
        chunk = islice(it, size)
        try:
            fist_element = next(chunk)
        except StopIteration:
            return
        yield chain((fist_element,), chunk)


# main function
def predict(path_ref, path_bam, path_model, path_out,
            ignore_txversion=False, nts=3, chunk_size=1024,
            rlen_min=25, rlen_max=35, threads=1):
    """
    load pre-trained model and predict P-site offsets
    TODO: need to consider both bam and sam
    """
    print('...load ref and model', file=sys.stderr)
    ref = read_fasta(path_ref, ignore_version=ignore_txversion)
    model = pickle.load(open(path_model, 'rb'))
    model.n_jobs = threads

    # prediction
    print('...Parse bam by chunk and predict', file=sys.stderr)
    batch_cols = ['qwidth'] + [f'{i}{j}' for i in ['s', 'e'] for j in range(2*nts)]
    features = ['qwidth'] + [f'{i}{j}_{k}' for i in ['s', 'e'] for j in range(2*nts) for k in 'ACGT']

    with pysam.AlignmentFile(path_bam) as bam:
        output = pysam.AlignmentFile(path_out, "wb", header=bam.header)
        for chunk in chunk_iter(bam, size=chunk_size):
            aligns = list(chunk)
            # filter missing tx
            aligns = [i for i in aligns if i.reference_name in ref]
            if len(aligns) == 0:
                continue
            qlen = np.zeros(len(aligns), dtype=int)
            seqs_start = [None] * len(aligns)
            seqs_end = [None] * len(aligns)
            for i, align in enumerate(aligns):
                seq = ref[align.reference_name]
                qlen[i] = align.query_alignment_length
                if align.reference_start - nts >= 0:
                    seq_left = seq[(align.reference_start - nts):(align.reference_start + nts)]
                else:
                    seq_left = seq[:(align.reference_start + nts)].rjust(2*nts, '-')
                if align.reference_end + nts <= len(seq):
                    seq_right = seq[(align.reference_end - nts):(align.reference_end + nts)]
                else:
                    seq_right = seq[(align.reference_end - nts):].ljust(2*nts, '-')
                if align.is_reverse:
                    seqs_start[i] = rev_comp(seq_right)
                    seqs_end[i] = rev_comp(seq_left)
                else:
                    seqs_start[i] = seq_left
                    seqs_end[i] = seq_right
            seqs_start = pd.DataFrame([list(i) for i in seqs_start])
            seqs_end = pd.DataFrame([list(i) for i in seqs_end])
            # index needed for edge case of one-row data frame
            batch = pd.DataFrame({'qwidth': qlen}, index=range(len(qlen)))
            batch = pd.concat([batch, seqs_start, seqs_end], axis = 1)
            batch.set_axis(batch_cols, axis=1, inplace=True)
            batch = pd.get_dummies(batch)
            X = np.zeros((batch.shape[0], 1 + 16 * nts), dtype=int)
            for k, column in enumerate(features):
                if column in list(batch.columns.values):
                    X[:,k] = batch.loc[:, column]
            p_sites = model.predict(X)
            for p, align in zip(p_sites, aligns):
                if align.query_alignment_length >= rlen_min and align.query_alignment_length <= rlen_max:
                    align.set_tag('PS', p, 'i')
                    output.write(align)
        output.close()
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog = 'python psite_predict.py',
        description = 'Predict P-site offset with the given model',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    # required arguments (input and output)
    parser.add_argument('ref', metavar='ref',
                        help='Reference transcriptome (fasta) matching the bam')
    parser.add_argument('bam', metavar='bam',
                        help='alignments of RPFs to reference transcriptome')
    parser.add_argument('model', metavar='model',
                        help='pre-trained RFC model')
    parser.add_argument('bam_out', metavar='bam_out',
                       help='output path of bam with PS (for P-site) tag')
    # optional arguments
    parser.add_argument('-c', '--chunk_size', type=int, default=65536,
                        help='chunk size for prediction batch')
    parser.add_argument('-l', '--lower', type=int, default=25,
                        help='lower bound for RPF mapped length')
    parser.add_argument('-u', '--upper', type=int, default=35,
                        help='upper bound for mapped read length')
    parser.add_argument('-i', '--ignore_txversion', action='store_true',
                        help='either to ignore trasncript version in ".\d+" format')
    parser.add_argument('-n', '--nts', type=int, default=3,
                        help='fanking nucleotides to consider at each side')
    parser.add_argument('-p', '--threads', type=int, default=1,
                        help='Number of threads used for prediction')
    args = parser.parse_args()
    predict(path_ref=args.ref,
            path_bam=args.bam,
            path_model=args.model,
            path_out=args.bam_out,
            ignore_txversion=args.ignore_txversion,
            nts=args.nts,
            chunk_size=args.chunk_size,
            rlen_min=args.lower,
            rlen_max=args.upper,
            threads=args.threads)
