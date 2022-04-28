import sys
import pickle

import click
import numpy as np
import pandas as pd
import pysam
from sklearn.ensemble import RandomForestClassifier
from utils import read_fasta, read_txinfo, get_txrep


@click.command()
@click.argument('path_ref', type=click.STRING)
@click.argument('path_bam', type=click.STRING)
@click.argument('path_model', type=click.STRING)
@click.argument('path_txinfo', type=click.STRING)
def train(path_ref, path_bam, path_model, path_txinfo,
          sep_txinfo='auto', type_ref='longest', path_exp=None, ignore_txversion=True,
          rlen_min=25, rlen_max=35, nts=3, threads=1):
    """
    Fit the random forest model

    PATH_REF: reference transcriptome (fasta) matching the bam
    PATH_BAM: alignments of RPFs to reference transcriptome
    PATH_model: path to save the fitted model
    PATH_TXINFO: basic info of transcripts
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