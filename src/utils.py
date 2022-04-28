from operator import index
import sys
import argparse
import re
import gzip
import pickle
from itertools import chain, islice

import click
import pandas as pd
import numpy as np
import pysam
from Bio import SeqIO


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


def chunk_iter(iterable, chunk_size=1024):
    """
    iterate by fixed block size
    
    Yield an iterator for each block. If the last chunk do not have enough
    items, all the remaining items is return as the last chunk.
    reference: https://stackoverflow.com/a/8998040/3926543
    """
    it = iter(iterable)
    while True:
        chunk = islice(it, chunk_size)
        try:
            fist_element = next(chunk)
        except StopIteration:
            return
        yield chain((fist_element,), chunk)

        
def time():
    return str(datetime.now())