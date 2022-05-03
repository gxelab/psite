import sys
import argparse
import re
import gzip
import pickle
import numpy as np
import pandas as pd
from Bio import SeqIO
import pysam
from sklearn.ensemble import RandomForestClassifier


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


def read_txinfo(path, sep='auto'):
    """Read transcript infomation"""
    if sep == 'auto':
        return pd.read_csv(path, sep=None, engine='python')
    else:
        return pd.read_csv(path, sep=sep)


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


def fit_model(path_ref, path_bam, path_model, path_txinfo,
              sep_txinfo='auto', type_ref='longest', path_exp=None, ignore_txversion=True,
              rlen_min=25, rlen_max=35, nts=3, threads=1):
    """
    Fit the random forest model

    note: soft-clips are not allowed in alignments, which will interfere
          with the extraction of flank sequences at two ends of an alignment.
    """
    # reference transcriptome
    print('...Load reference fasta', file=sys.stderr)
    ref = read_fasta(path=path_ref, ignore_version=ignore_txversion)
    # gene info
    print('...Load gene information', file=sys.stderr)
    tx_info = read_txinfo(path_txinfo, sep_txinfo)
    
    # get representative transcripts
    print('...Get representative transcript per gene', file=sys.stderr)
    txrep = get_txrep(tx_info, type_ref, path_exp, ignore_txversion)

    print('...keep items common in ref fasta and txinfo', file=sys.stderr)
    # check whether tx length in `ref` is consistent with that in `txrep`
    txrep = txrep[txrep['tx_name'].isin(ref.keys())]
    ref_len = np.array([len(ref[i]) for i in txrep['tx_name']])
    txrep = txrep[txrep['tx_len'] == ref_len]
    ref = {k: ref[k] for k in txrep['tx_name']}
    
    # get the position of start codon and stop codon of each transcript
    txrep['start_codon'] = txrep['utr5_len'] + 1
    txrep['stop_codon'] = txrep['utr5_len'] + txrep['cds_len'] - 2
    txrep = txrep.set_index('tx_name', drop=False)  # index for faster value access
    
    # check
    print('...check metada', file=sys.stderr)
    print(txrep.head(), file=sys.stderr)
    print(list(ref.keys())[:6], file = sys.stderr)
    
    # parse alignments
    print('...parse alignments', file=sys.stderr)
    reads = set()
    col_qwidth = list()
    col_sqleft = list()
    col_sqright = list()
    col_label = list()

    with pysam.AlignmentFile(path_bam) as bam:
        for align in bam:
            if align.query_name in reads:
                continue
            # alignment qc
            if align.is_reverse:
                continue
            if align.query_alignment_length < rlen_min or align.query_alignment_length > rlen_max:
                continue
            if align.reference_name not in txrep['tx_name']:
                continue
            tx_info = txrep.loc[align.reference_name]
            seq = ref[align.reference_name]
            if align.reference_start < nts or align.reference_end > tx_info['tx_len'] - nts:
                continue
            reads.add(align.query_name)
            # keep alignments overlapping with start codon or stop codon (1-based)
            # reference_start: 0-based; reference_end: 1-based rightmost coordinate
            dist_start = tx_info['start_codon'] - align.reference_start - 1
            dist_stop = tx_info['stop_codon'] - align.reference_start - 4
            if dist_start >= 11 and dist_start <= 14:
                label = dist_start
            elif dist_stop >= 10 and dist_stop <= 13:
                label = dist_stop
            else: # not overlapping with start/stop in given offset range
                continue
            col_qwidth.append(align.query_alignment_length)
            col_sqleft.append(seq[(align.reference_start - nts):(align.reference_start + nts)])
            col_sqright.append(seq[(align.reference_end - nts):(align.reference_end + nts)])
            col_label.append(label)

    # prepare training data
    print('...prepare training data', file=sys.stderr)
    sqleft = pd.DataFrame([list(i) for i in col_sqleft]).add_prefix('s')
    sqright = pd.DataFrame([list(i) for i in col_sqright]).add_prefix('e')
    training_data = pd.DataFrame({'label': col_label, 'qwidth': col_qwidth})
    training_data = pd.concat([training_data, sqleft, sqright], axis = 1)
    # prevent ambiguous base 'N' being used in training data
    features = pd.get_dummies(training_data)
    scols = ['label', 'qwidth'] 
    scols += [f'{i}{j}_{k}' for i in ['s', 'e'] for j in range(2*nts) for k in 'ACGT']
    cols_to_remove = [ col for col in features.columns if col not in scols]
    features.drop(columns=cols_to_remove, inplace=True)

    X = features.drop(columns='label').to_numpy()
    y = features[['label']].to_numpy().flatten()
    
    # fit model
    print('...fit model', file=sys.stderr)
    rfc = RandomForestClassifier(n_estimators=200, max_features='sqrt', n_jobs=threads)
    rfc.fit(X, y)
    pickle.dump(rfc, open(path_model, 'wb'))
    print('...Done!', file=sys.stderr)
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog = 'python psite_train.py',
        description = 'Fit Random Forest Model for P-site offset estimation',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    # required arguments (input and output)
    parser.add_argument('ref', metavar='ref',
                        help='reference transcriptome (fasta) matching the bam')
    parser.add_argument('path_bam', metavar='bam',
                        help='alignments of RPFs to reference transcriptome')
    parser.add_argument('txinfo', metavar='txinfo',
                        help='basic info of transcripts')
    parser.add_argument('model', metavar='model',
                       help='path to save the fitted model')
    # optional arguments
    parser.add_argument('-s', '--sep_txinfo', type=str, default='auto',
                        help='field delimiter of the txinfo file')
    parser.add_argument('-e', '--expression', type=str, default=None,
                        help='lower bound for RPF mapped length')
    parser.add_argument('-l', '--lower', type=int, default=25,
                        help='lower bound for RPF mapped length')
    parser.add_argument('-u', '--upper', type=int, default=35,
                        help='upper bound for mapped read length')
    parser.add_argument('-t', '--type_ref', type=str, default='longest',
                        choices=('longest', 'principal', 'kallisto', 'salmon'),
                        help='type of representative transcripts')
    parser.add_argument('-i', '--ignore_txversion', action='store_true',
                        help='either to ignore trasncript version in ".\d+" format')
    parser.add_argument('-n', '--nts', type=int, default=3,
                        help='fanking nucleotides to consider at each side')
    parser.add_argument('-p', '--threads', type=int, default=1,
                        help='Number of threads used for model fitting')
    args = parser.parse_args()
    fit_model(
        path_ref=args.ref,
        path_bam=args.bam,
        path_txinfo=args.txinfo,
        path_model=args.model,
        sep_txinfo=args.sep_txinfo,
        type_ref=args.type_ref,
        ignore_txversion=args.ignore_txversion,
        path_exp=args.expression,
        rlen_min=args.lower,
        rlen_max=args.upper,
        nts=args.nts,
        threads=args.threads)