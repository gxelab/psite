import sys
import argparse
import numpy as np
import pysam
import pyBigWig


def bw_write_chrom(bw, chrom, cov):
    """write the coverage of a chromosome into the bw file"""
    # note: np.array supports slice/modify by index vector
    run_end = np.append(np.where(np.diff(cov) != 0), len(cov)-1) + 1  # 1-based ends
    run_start = run_end - np.diff(np.append(0, run_end))  # 0-based starts
    run_value = cov[run_start]
    # ignore 0 values
    non_zero  = run_value > 0
    run_start = run_start[non_zero]
    run_end   = run_end[non_zero]
    run_value = run_value[non_zero]

    if len(run_value) > 0:
        bw.addEntries([chrom] * np.sum(non_zero),
                      run_start, ends = run_end, values = run_value)
    return


def coverage(path_bam, prefix, rlen_min=25, rlen_max=35, mapq_min=10, ignore_supp=True):
    """
    calculate the coverage for plus strand and minus strand seperately

    
    """
    # try to open bam file
    bam = pysam.AlignmentFile(path_bam)
    # create bw file and add header

    bw_fw = pyBigWig.open(f'{prefix}_fw.bw', 'w')
    bw_rc = pyBigWig.open(f'{prefix}_rc.bw', 'w')
    bw_header = [(i, bam.get_reference_length(i)) for i in bam.references]
    bw_fw.addHeader(bw_header)
    bw_rc.addHeader(bw_header)

    # processing the bam file
    refname_curr = None  # name of current reference
    for align in bam:
        # mapping quality and read length QC
        if align.mapping_quality < mapq_min:
            continue
        if ignore_supp and align.is_supplementary:
            continue
        if align.query_alignment_length < rlen_min or align.query_alignment_length > rlen_max:
            continue

        # write bw before processing a new chromosome
        if align.reference_name != refname_curr:
            if refname_curr != None:
                bw_write_chrom(bw_fw, refname_curr, cov_plus)
                bw_write_chrom(bw_rc, refname_curr, cov_minus)
            refname_curr = align.reference_name
            cov_plus = np.zeros(bam.get_reference_length(align.reference_name))
            cov_minus = np.zeros(bam.get_reference_length(align.reference_name))
        # get P-site offset for read of this length
        # offset = 0 means no offset and count the 5' most position
        offset = align.get_tag('PS')
        # store coverage of current read to the coverage vectors
        # coverage vectors are 0-based
        # get_reference_positions() return 0-based positions on the reference
        # (https://pysam.readthedocs.io/en/latest/api.html#introduction)
        if align.is_reverse:
            cov_minus[ align.get_reference_positions()[-(offset + 1)] ] += 1
        else:
            cov_plus[ align.get_reference_positions()[offset] ] += 1
    # write last chromosome if not empty (note: at least one of them is not empty)
    if np.sum(cov_plus) > 0:
        bw_write_chrom(bw_fw, align.reference_name, cov_plus)
    if np.sum(cov_minus) > 0:
        bw_write_chrom(bw_rc, align.reference_name, cov_minus)

    # close files
    bam.close()
    bw_fw.close()
    bw_rc.close()
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog = 'python psite_cov.py',
        description = 'Calculate P-site coverage',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    # required arguments (input and output)
    parser.add_argument('bam', metavar='bam',
                        help='aligment bam file with the PS tag (for P-site offset)')
    parser.add_argument('prefix', metavar='prefix',
                        help='output prefix of P-site coverage tracks in bigWig format')
    # optional arguments
    parser.add_argument('-l', '--lower', type=int, default=25,
                        help='lower bound for RPF mapped length')
    parser.add_argument('-u', '--upper', type=int, default=35,
                        help='upper bound for RPF mapped length')
    parser.add_argument('-q', '--mapq_min', type=int, default=35,
                        help='minimum mapping quality')
    parser.add_argument('-i', '--ignore_supp', action='store_true',
                        help='Whether to ignore supplementary alignments')
    args = parser.parse_args()
    coverage(path_bam=args.bam, prefix=args.prefix,
             rlen_min=args.lower, rlen_max=args.upper,
             mapq_min=args.mapq_min, ignore_supp=args.ignore_supp)
