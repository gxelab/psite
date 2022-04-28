import sys
import pickle

import click
import pandas as pd
import numpy as np
import pysam
from utils import read_fasta, chunk_iter, rev_comp


CLICK_CS = dict(help_option_names=['-h', '--help'], show_default=True)
@click.command(context_settings=CLICK_CS)
@click.argument('path_ref', type=click.STRING)
@click.argument('path_bam', type=click.STRING)
@click.argument('path_model', type=click.STRING)
@click.argument('path_out', type=click.STRING)
@click.option('-f', '--out_format', default='bam',
              type=click.Choice(['bam', 'sam']),
              help='P-site alignemnt output format')
@click.option('-c', '--chunk_size', type=click.INT, default=65536,
              help='chunk size for prediction batch')
@click.option('-i', '--ignore_txversion', is_flag=True, default=False,
              help='either to ignore trasncript version in ".\d+" format')
@click.option('-l', '--rlen_min', type=click.INT, default=25,
              help='lower bound for RPF mapped length')
@click.option('-u', '--rlen_max', type=click.INT, default=35,
              help='upper bound for mapped read length')
@click.option('-n', '--nts', type=click.INT, default=3,
              help='fanking nucleotides to consider at each side')
@click.option('-p', '--threads', type=click.INT, default=1,
              help='Number of threads used for prediction')
def pbam(path_ref, path_bam, path_model, path_out, out_format='bam',
            rlen_min=25, rlen_max=35, ignore_txversion=False, nts=3,
            chunk_size=1024, threads=1):
    """
    generate bam with only P-site regions
    
    \b
    path_ref   : reference transcriptome (fasta) matching the bam
    path_bam   : alignments of RPFs to reference transcriptome
    path_model : path to save the fitted model
    path_out   : output path of bam with P-site regions only
    """
    print('...load ref and model', file=sys.stderr)
    # gzipped fasta is detected automatically
    ref = read_fasta(path_ref, ignore_version=ignore_txversion)
    model = pickle.load(open(path_model, 'rb'))
    model.n_jobs = threads
    
    # prediction
    print('...Parse bam by chunk and predict', file=sys.stderr)
    batch_cols = ['qwidth'] + [f'{i}{j}' for i in ['s', 'e'] for j in range(2*nts)]
    features = ['qwidth'] + [f'{i}{j}_{k}' for i in ['s', 'e'] for j in range(2*nts) for k in 'ACGT']

    with pysam.AlignmentFile(path_bam) as bam:
        if out_format == 'bam':
            output = pysam.AlignmentFile(path_out, "wb", header=bam.header)
        else:
            output =  pysam.AlignmentFile(path_out, "w", header=bam.header)
        for chunk in chunk_iter(bam, chunk_size):
            # print('...reading new chunk: ' + time(), file=sys.stderr)
            aligns = list(chunk)
            # filter missing tx
            aligns = [i for i in aligns if i.reference_name in ref]
            if len(aligns) == 0:
                continue
            qlen = np.zeros(len(aligns), dtype=int)
            seqs_start = [None] * len(aligns)
            seqs_end = [None] * len(aligns)
            for i, align in enumerate(aligns):
                qlen[i] = align.query_alignment_length
                if align.reference_start - nts >= 0:
                    seq_left = ref[align.reference_name][(align.reference_start - nts):(align.reference_start + nts)]
                else:
                    seq_left = ref[align.reference_name][:(align.reference_start + nts)].rjust(2*nts, '-')
                if align.reference_end + nts <= len(ref[align.reference_name]):
                    seq_right = ref[align.reference_name][(align.reference_end - nts):(align.reference_end + nts)]
                else:
                    seq_right = ref[align.reference_name][(align.reference_end - nts):].ljust(2*nts, '-')
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
                a = pysam.AlignedSegment()
                a.query_name = align.query_name
                a.flag = align.flag
                a.reference_id = align.reference_id
                a.mapping_quality = align.mapping_quality
                a.tags = align.tags
                if align.is_reverse:
                    if align.query_alignment_length < rlen_min or align.query_alignment_length > rlen_max:
                        continue
                    # note: both reference_start and reference_positions are 0-based
                    a.query_sequence  = align.query_alignment_sequence[-(p + 1)]
                    a.reference_start = align.get_reference_positions()[-(p + 1)]
                    a.cigar = ((0, 1),)
                    a.query_qualities = [ align.query_alignment_qualities[-(p + 1)] ]
                else:
                    a.query_sequence  = align.query_alignment_sequence[p]
                    a.reference_start = align.get_reference_positions()[p]
                    a.cigar = ((0, 1),)
                    a.query_qualities  = [ align.query_alignment_qualities[p] ]
                output.write(a)
        output.close()
    return