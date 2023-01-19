import sys
import click
import pysam
from psite.utils import read_fasta, chunk_iter, rev_comp, CLICK_CS
# from utils import read_fasta, chunk_iter, rev_comp, CLICK_CS


@click.command(context_settings=CLICK_CS)
@click.argument('path_bam', type=click.STRING)
@click.argument('path_out', type=click.STRING)
@click.option('-l', '--rlen_min', type=click.INT, default=27,
              help='lower bound for mapped read length')
@click.option('-u', '--rlen_max', type=click.INT, default=35,
              help='upper bound for mapped read length')
@click.option('-n', '--nucleotides', type=click.INT, default=12,
              help='fixed global offset value')
def setp(path_bam, path_out, rlen_min=None, rlen_max=None, nucleotides=12):
    """
    set global fixed P-site offset tag
    
    \b
    path_bam   : alignments of RPFs to reference transcriptome
    path_out   : output path of bam with PS (for P-site) tag 
    """
    with pysam.AlignmentFile(path_bam) as bam:
        output = pysam.AlignmentFile(path_out, "wb", header=bam.header)
        for align in bam:
            if align.query_alignment_length >= rlen_min and align.query_alignment_length <= rlen_max:
                align.set_tag('PS', nucleotides, 'i')
                output.write(align)
        output.close()
    return


if __name__ == '__main__':
    setp()
