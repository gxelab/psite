import sys
import pickle
from multiprocessing import Pool

import click
import numpy as np
import pandas as pd
import pysam
from sklearn.ensemble import RandomForestClassifier
from utils import read_fasta, read_txinfo, strip_version, CLICK_CS


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


@click.command(context_settings=CLICK_CS)
@click.argument('path_ref', type=click.STRING)
@click.argument('path_bam', type=click.STRING)
@click.argument('path_model', type=click.STRING)
@click.argument('path_txinfo', type=click.STRING)
@click.option('-s', 'sep_txinfo', type=click.STRING, default='auto',
              help='field delimiter of the txinfo file')
@click.option('-t', '--type_ref', default='longest',
              type=click.Choice(['longest', 'principal', 'kallisto', 'salmon']),
              help='type of representative transcripts')
@click.option('-e', '--path_exp', type=click.STRING, default=None,
              help='lower bound for RPF mapped length')
@click.option('-i', '--ignore_txversion', is_flag=True, default=False,
              help='whether to ignore trasncript version in ".\d+" format')
@click.option('-l', '--rlen_min', type=click.INT, default=25,
              help='lower bound for RPF mapped length')
@click.option('-u', '--rlen_max', type=click.INT, default=35,
              help='upper bound for mapped read length')
@click.option('-n', '--nts', type=click.INT, default=3,
              help='fanking nucleotides to consider at each side')
@click.option('-p', '--threads', type=click.INT, default=1,
              help='Number of threads used for model fitting')
def train(path_ref, path_bam, path_model, path_txinfo,
          sep_txinfo='auto', type_ref='longest', path_exp=None, ignore_txversion=True,
          rlen_min=25, rlen_max=35, nts=3, threads=1):
    """
    train a random forest model of p-site offsets

    \b
    path_ref   : reference transcriptome (fasta) matching the bam
    path_bam   : alignments of RPFs to reference transcriptome
    path_model : path to save the fitted model
    path_txinfo: transcriptome annotation
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
    
    def process_align(align):
        """helper function to process an alignment"""
        # print(align.query_name, file=sys.stderr)
        # alignment qc
        if align.is_reverse:
            return ()
        if align.query_alignment_length < rlen_min or align.query_alignment_length > rlen_max:
            return ()
        if align.reference_name not in txrep['tx_name']:
            return ()
        tx_info = txrep.loc[align.reference_name]
        if align.reference_start < nts or align.reference_end > tx_info['tx_len'] - nts:
            return ()
        # keep alignments overlapping with start codon or stop codon (1-based)
        # reference_start: 0-based; reference_end: 1-based rightmost coordinate
        dist_start = tx_info['start_codon'] - align.reference_start - 1
        dist_stop = tx_info['stop_codon'] - align.reference_start - 4
        if dist_start >= 11 and dist_start <= 14:
            label = dist_start
        elif dist_stop >= 10 and dist_stop <= 13:
            label = dist_stop
        else: # not overlapping with start/stop in given offset range
            return ()
        seq = ref[align.reference_name]
        sqleft = seq[(align.reference_start - nts):(align.reference_start + nts)]
        sqright = seq[(align.reference_end - nts):(align.reference_end + nts)]
        return (align.query_name, label, align.query_alignment_length, sqleft, sqright)


    # parse alignments (consider using chunksize)
    print(f'...parse alignments using {threads} threads', file=sys.stderr)
    with pysam.AlignmentFile(path_bam, 'rb') as bam:
        pool = Pool(processes=threads)
        aligns = pool.imap(process_align, bam, chunksize=1)
        print(next(aligns))
        aligns = [i for i in aligns if i]  # keep non-empty entries
        pool.close()
        pool.join()

    # prepare training data
    print('...prepare training data', file=sys.stderr)
    aligns = pd.DataFrame(aligns)
    aligns.set_axis(['read_name', 'label', 'qwidth', 'sqleft', 'sqright'], axis=1, inplace=True)
    aligns.drop_duplicates(subset=['read_name'], inplace=True)

    sqleft = pd.DataFrame([list(i) for i in aligns['sqleft']]).add_prefix('s')
    sqright = pd.DataFrame([list(i) for i in aligns['sqright']]).add_prefix('e')
    aligns.drop(['read_name', 'sqleft', 'sqright'], inplace=True)
    aligns = pd.concat([aligns, sqleft, sqright], axis = 1)

    # prevent ambiguous base 'N' being used in training data
    features = pd.get_dummies(aligns)
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
    train()