from datetime import datetime
import re
import gzip
from itertools import chain, islice

import pandas as pd
from Bio import SeqIO

# click setting ###############################################################
CLICK_CS = dict(help_option_names=['-h', '--help'], show_default=True)

# reading #####################################################################
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


def chunk_iter(iterable, size=1024):
    """
    iterate by fixed block size
    
    Yield an iterator for each block. If the last chunk do not have enough
    items, all the remaining items is return as the last chunk.
    reference: reclosedev, 2012, https://stackoverflow.com/a/8998040/3926543
    """
    it = iter(iterable)
    while True:
        chunk = islice(it, size)
        try:
            fist_element = next(chunk)
        except StopIteration:
            return
        yield chain((fist_element,), chunk)


# helper functions ############################################################
WC_PAIRING = str.maketrans('ACGTN', 'TGCAN')
def rev_comp(s):
    """reverse complementation DNA string"""
    return s[::-1].translate(WC_PAIRING)

        
def time():
    return str(datetime.now())
