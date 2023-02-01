import sys
from pickle import dump
from random import random
from re import sub
import click
import numpy as np
import pandas as pd
import seaborn as sns
import pysam
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.utils.validation import check_X_y, check_array, check_is_fitted
from sklearn.utils.multiclass import unique_labels
from psite.utils import CLICK_CS, read_fasta, read_txinfo, strip_version
# from utils import CLICK_CS, read_fasta, read_txinfo, strip_version


class DistArgMax(BaseEstimator, ClassifierMixin):
    """
    a learner that assigns offset based on the most frequent distance
    to start codon in the training data.
    """
    def __init__(self):
        self.offsets = dict()

    def fit(self, X, y):
        X, y = check_X_y(X, y)
        self.classes_ = unique_labels(y)
        X_qwidth = X[:, 0]  # the first column should be qwidth
        for qw in unique_labels(X_qwidth):
            val, cnt = np.unique(y[X_qwidth == qw], return_counts=True)
            self.offsets[qw] = val[cnt.argmax()]
        return self

    def predict(self, X):
        check_is_fitted(self)
        X = check_array(X)
        return np.array([self.offsets[i] for i in X[:, 0]])


def frame_test(test, psite):
    """
    calculate proportions of psites in each frame
    """
    frame = (test.tstart + psite - test.start_codon) % 3
    return frame.value_counts() / frame.size


def get_txrep(txinfo, type_rep='longest', path_exp=None, ignore_version=False):
    """Get representative transcript per gene"""
    if 'gene_biotype' in txinfo.columns:
        txinfo = txinfo.loc[txinfo['gene_biotype'] == 'protein_coding']
    if 'transcript_biotype' in txinfo.columns:
        txinfo = txinfo.loc[txinfo['transcript_biotype'] == 'protein_coding']
    if type_rep == 'longest':
        txrep = txinfo.sort_values(['gene_id', 'tx_len'], ascending=[True, False])
        txrep = txrep.groupby('gene_id').first().reset_index()
    elif type_rep == 'principal':
        txrep = txinfo.loc[(txinfo['txtype'] == 'principal')]
    elif type_rep == 'salmon':   # salmon_output_dir/quant.sf
        try:
            tx_quant = pd.read_csv(path_exp, sep='\t')
            if ignore_version:
                tx_quant['Name'] = tx_quant['Name'].map(strip_version)
            tx_quant = tx_quant.merge(
                txinfo[['gene_id', 'tx_name']], how='inner', left_on='Name', right_on='tx_name')
            tx_quant = tx_quant.sort_values(by=['gene_id', 'TPM'], ascending=[True, False])
            tx_quant = tx_quant.groupby('gene_id').first().reset_index()
            txrep = txinfo.loc[(txinfo['tx_name'].isin(tx_quant['Name']))]
        except:
            print('input salmon quant results incorrect', file=sys.stderr)
            exit()
    elif type_rep == 'kallisto':  # kallisto_output_dir/abundance.tsv
        try:
            tx_quant = pd.read_csv(path_exp, sep='\t')
            if ignore_version:
                tx_quant['target_id'] = tx_quant['target_id'].map(strip_version)
            tx_quant = tx_quant.merge(
                txinfo[['gene_id', 'tx_name']], how='inner', left_on='target_id', right_on='tx_name')
            tx_quant = tx_quant.sort_values(by=['gene_id', 'tpm'], ascending=[True, False])
            tx_quant = tx_quant.groupby('gene_id').first().reset_index()
            txrep = txinfo.loc[(txinfo['tx_name'].isin(tx_quant['target_id']))]
        except:
            print('input kallisto quant results incorrect', file=sys.stderr)
            exit()
    else:
        print('Incorrect txinfo_rep option!', file=sys.stderr)
        exit()
    return txrep


def extract_features(path_bam, ref, nts=3, frac=1.0):
    """
    Extract and save features for model training and testing

    path_bam: alignments of RPFs to reference transcriptome
    ref: reference transcriptome (fasta) matching the bam
    nts: number of nucleotides to include from each side of the 5' end
    frac: frac of alignments to read from bam
    """
    alignments = []
    with pysam.AlignmentFile(path_bam) as bam:
        for align in bam:
            if frac < 1 and random() > frac:
                continue
            if align.reference_name not in ref:
                continue
            if align.is_reverse:
                # alignments in transcriptome bam should be in forward strand
                continue
            # .reference_start: 0-based leftmost coordinate
            # .reference_end: reference_end points to one past the last aligned residue.
            # so that reference_end - reference_start = reference_length
            seq = ref[align.reference_name]
            if align.reference_start - nts >= 0:
                flank_5p = seq[(align.reference_start - nts):(align.reference_start + nts)]
            else:
                flank_5p = seq[:(align.reference_start + nts)].rjust(2*nts, '-')
            
            out = [align.query_name, align.reference_name, align.reference_start + 1,
                   align.query_alignment_length] + list(flank_5p)
            alignments.append(out)
    cols = ['rname', 'tx_name', 'tstart', 'qwidth'] + [f's{i}' for i in range(2*nts)]
    return pd.DataFrame(alignments, columns=cols)


@click.command(context_settings=CLICK_CS)
@click.argument('path_ref', type=click.STRING)
@click.argument('path_bam', type=click.STRING)
@click.argument('output_prefix', type=click.STRING)
@click.argument('path_txinfo', type=click.STRING)
@click.option('-t', '--type_rep', default='longest',
              type=click.Choice(['longest', 'principal', 'kallisto', 'salmon']),
              help='type of representative transcripts')
@click.option('-e', '--path_exp', type=click.STRING, default=None,
              help='path of transcript expression quant results')
@click.option('-i', '--ignore_txversion', is_flag=True, default=False,
              help='ignore transcript version in ".\d+" format')
@click.option('-n', '--nts', type=click.INT, default=3,
              help='flanking nucleotides to consider at each side')
@click.option('-f', '--frac', type=click.FLOAT, default=1.0,
              help='fraction of alignments for training (for large datasets)')
@click.option('--offset_min', type=click.INT, default=10,
              help='lower bound of distance between RPF 5p and start codon')
@click.option('--offset_max', type=click.INT, default=14,
              help='upper bound of distance between RPF 5p and start codon')
@click.option('-d', '--max_depth', type=click.INT, default=3,
              help='max depth of trees')
@click.option('-m', '--min_samples_split', type=click.INT, default=6,
              help='min number of alignments required to split an internal node')
@click.option('-k', '--keep', is_flag=True, default=False,
              help='whether to keep intermediate results')
def train(path_ref, path_bam, output_prefix, path_txinfo,
          type_rep='longest', path_exp=None, ignore_txversion=True,
          nts=3, keep=False, frac=1,  offset_min=11, offset_max=14,
          max_depth=3, min_samples_split=6):
    """
    train a model for P-site offset prediction

    \b
    path_ref     : reference transcriptome (fasta) matching the bam
    path_bam     : alignments of RPFs to reference transcriptome
    output_prefix: output prefix of fitted models and logs
    path_txinfo  : transcriptome annotation
    """

    log = open(f'{output_prefix}.log', 'wt')
    print((
        f'# TxInfo: {path_txinfo}\n'
        f'# Transcript quant: {path_exp}\n'
        f'# Alignment file: {path_bam}\n'
        f'# nts                     = 3\n'
        f'# train_offset_min        = {offset_min}\n'
        f'# train_offset_max        = {offset_max}\n'
        f'# train_max_depth         = {max_depth}\n'
        f'# train_min_samples_split = {min_samples_split}'), file=log)

    # load data ===============================================================
    print('# ...load reference fasta', file=log)
    ref = read_fasta(path=path_ref, ignore_version=ignore_txversion)

    print('# ...load gene information', file=log)
    tx_info = read_txinfo(path_txinfo, '\t')
    
    # exclude MT and fly transposable elements
    tx_info = tx_info[~tx_info.chrom.isin(['MT', 'mitochondrion_genome'])]
    if (tx_info.gene_id.str.startswith('FBgn')).sum() > 0:
        tx_info = tx_info[tx_info.gene_id.str.startswith('FBgn')]
    
    print('# ...get representative transcript per gene', file=log)
    txrep = get_txrep(tx_info, type_rep, path_exp, ignore_txversion)

    print('# ...keep items common in ref fasta and txinfo', file=log)
    txrep = txrep[txrep['tx_name'].isin(ref.keys())]
    ref_len = np.array([len(ref[i]) for i in txrep['tx_name']])
    txrep = txrep[txrep['tx_len'] == ref_len]
    ref = {k: ref[k] for k in txrep['tx_name']}
    
    # get the position of start codon and stop codon of each transcript
    txrep['start_codon'] = txrep['utr5_len'] + 1
    txrep['stop_codon'] = txrep['utr5_len'] + txrep['cds_len'] - 2
    
    print('# ...check metada', file=log)
    print(sub(r'(^|\n)', r'\1# ', str(txrep.head())), file=log)
    print('# ', list(ref.keys())[:6], file = log)
    
    print('# ...parse alignments', file=log)
    alignments = extract_features(path_bam, ref, nts=nts, frac=frac)
    alignments = alignments[alignments.tx_name.isin(txrep.tx_name)]
    alignments.drop_duplicates(subset='rname', inplace=True)
    alignments.drop(columns='rname', inplace=True)

    # change nucleotide columns to categorical, so that there are always four columns when run
    # pd.get_dummies for each column, even of one or more of A/C/G/T didn't appear in this column.
    nuc_cols = [f's{i}' for i in range(2*nts)]
    for i in nuc_cols:
        alignments.loc[:, i] = pd.Categorical(alignments[i], categories=['A', 'C', 'G', 'T'])

    data = pd.merge(alignments, txrep[['tx_name', 'start_codon', 'stop_codon']], on='tx_name')
    data = data.assign(
        dist_start=lambda x: x.start_codon - x.tstart,
        dist_stop=lambda x: x.stop_codon - x.tstart - 3)
    del alignments

    if keep:
        data.to_csv(f'', sep='\t')

    # QC plot =================================================================
    print('# ...draw QC plots of read length and offset distribution', file=log)
    qwidth_freq = data.qwidth.value_counts().sort_index()
    qwidth_freq = qwidth_freq[(qwidth_freq.index >= 25) & (qwidth_freq.index <= 40)]
    qwidth_freq = qwidth_freq / qwidth_freq.sum()
    ax = qwidth_freq.plot(kind='bar')
    ax.figure.savefig(f'{output_prefix}.qwidth_distribution.pdf')

    # filter qwidth
    qwidth_range = qwidth_freq.index[qwidth_freq > 0.05]
    print(f'# Read lengths used: {list(qwidth_range)}', file=log)
    data = data[data.qwidth.isin(qwidth_range)]

    # plot 5' end distribution around start and stop codonW
    tmp = data[data.dist_start.between(8, 20)].copy()
    g = sns.FacetGrid(tmp, col="qwidth", col_wrap=5, sharex=False, sharey=False, height=3)
    g.map(sns.histplot, 'dist_start', binwidth=0.25)
    g.set(xticks=[8, 10, 12, 14, 16, 18, 20])
    g.savefig(f'{output_prefix}.offset_distribution.pdf')
    del ax, g, tmp

    # use alignments overlapping start or stop codons for training
    print('# ...filter aligments for training', file=log)
    train_start = data[data.dist_start.between(offset_min, offset_max)].copy()
    train_start['dist'] = train_start.dist_start
    train_stop = data[data.dist_stop.between(offset_min, offset_max)].copy()
    train_stop['dist'] = train_stop.dist_stop
    train = pd.concat([train_start, train_stop], axis=0)
    del train_start, train_stop

    # train different learners ================================================
    print('# ...train DAM and GBT', file=log)
    train_features = ['qwidth'] + nuc_cols
    X = pd.get_dummies(train[train_features], drop_first=True)
    y = train.dist.to_numpy()

    dam = DistArgMax()
    dam.fit(X, y)
    dam.qwidth_range = qwidth_range.to_numpy()
    dam.nts = nts
    dump(dam, open(f'{output_prefix}.dam.pickle', 'wb'))

    gbt = GradientBoostingClassifier(n_estimators=100, learning_rate=0.1,
        min_samples_split=min_samples_split, max_depth=max_depth)
    gbt.fit(X, y)
    gbt.qwidth_range = qwidth_range.to_numpy()
    gbt.nts = nts
    dump(gbt, open(f'{output_prefix}.gbt.pickle', 'wb'))

    # reads within CDS are used for training
    cds_aligns = data[(data.dist_start < 9) & (data.dist_stop > 15)]
    cds_X = pd.get_dummies(cds_aligns[train_features], drop_first=True)

    dam_cds_yhat = dam.predict(cds_X)
    gbt_cds_yhat = gbt.predict(cds_X)

    print(pd.DataFrame({
        'none': frame_test(cds_aligns, 0),
        'base': frame_test(cds_aligns, dam_cds_yhat),
        'GBT': frame_test(cds_aligns, gbt_cds_yhat),
    }), file=log)

    print(pd.DataFrame({
        'feat': X.columns,
        'GBT': gbt.feature_importances_
    }), file=log)

    print('# Done!', file=log)
    log.close()
    return


if __name__ == '__main__':
    train()
