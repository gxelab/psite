# PSite


[![DOI](https://zenodo.org/badge/474568909.svg)](https://zenodo.org/badge/latestdoi/474568909)

`PSite` is a python package that predicts P-site offsets for footprints generated in ribosome profiling using a Gradient Boosting Trees (GBT) model trained with footprints around both annotated start and stop codons. `PSite` can report estimated P-site offsets in two manners:

- append a `PS` tag to each original alignment in `SAM` or `BAM` format, without any other modifications;
- output a new `BAM` file of the alignments of P-site locations only;

To demonstrate the usage of the `PS` tag, `PSite` also has a `coverage` module that performs genome-wide calculation of P-site coverage of ribosome footprints at nucleotide-resolution.

### Dependency
- `numpy` >= 1.21.2
- `pandas` >= 1.3.4
- `biopython` >= 1.79
- `sklearn` >= 1.1.1
- `pysam` >= 0.17.0
- `pyBigWig` >= 0.3.18
- `click` >= 8.1.2
- `seaborn` >= 0.11.0


### Install and uninstall
To install `PSite`, download the package [tarball](https://github.com/gxelab/psite/releases) from the release page and run
```bash
pip install psite-version-py3-none-any.whl
```

To uninstall it, simply run
```
pip uninstall psite
```

### Build distributions from source
Run the following command in source directory
```bash
python3 -m build
```

---------------------------------------

### Usage
`PSite` is designed to be used from command line on Unix-like operating systems such as Linux or macOS.

```bash
$ psite -h
Usage: psite [OPTIONS] COMMAND [ARGS]...

  main interface

Options:
  -h, --help  Show this message and exit.

Commands:
  coverage  calculate the coverage for plus strand and minus strand...
  pbam      generate bam with only P-site regions
  predict   load pre-trained model and predict P-site offsets
  train     train a model for P-site offset prediction
```

#### `train`
This is the core module that trains the GBT model for P-site offset prediction. It requires transcriptome alignments (`PATH_BAM`) and the corresponding sequences of all transcripts (`PATH_REF`). The required bam can be generated by mapping footprints to the reference genome using [STAR](https://github.com/alexdobin/STAR) and output transcriptome alignments with parameter `--quantMode TranscriptomeSAM`. The trained model is saved in `pickle` format for later use.

```bash
$ psite train -h
Usage: psite train [OPTIONS] PATH_REF PATH_BAM OUTPUT_PREFIX PATH_TXINFO

  train a model for P-site offset prediction

  path_ref     : reference transcriptome (fasta) matching the bam
  path_bam     : alignments of RPFs to reference transcriptome
  output_prefix: output prefix of fitted models and logs
  path_txinfo  : transcriptome annotation

Options:
  -t, --type_rep [longest|principal|kallisto|salmon]
                                  type of representative transcripts
                                  [default: longest]
  -e, --path_exp TEXT             path of transcript expression quant results
  -i, --ignore_txversion          whether to ignore trasncript version in
                                  ".\d+" format  [default: False]
  -n, --nts INTEGER               fanking nucleotides to consider at each side
                                  [default: 3]
  -f, --frac FLOAT                fraction of alignments for training (for
                                  large datasets)  [default: 1.0]
  --offset_min INTEGER            lower bound of distance between RPF 5p and
                                  start codon  [default: 10]
  --offset_max INTEGER            upper bound of distance between RPF 5p and
                                  start codon  [default: 14]
  -d, --max_depth INTEGER         number of threads used for model fitting
                                  [default: 3]
  -m, --min_samples_split INTEGER
                                  min number of alignments required to split
                                  an internal node  [default: 6]
  -k, --keep                      whether to to keep intermediate results
                                  [default: False]
  -h, --help                      Show this message and exit.  [default:
                                  False]
```

#### `predict`
This module predict P-site for each alignment using a pre-trained model and append a `PS` tag (for "P-site") to the original alignment. The input can be either genomic alignments or transcriptomic alignments.

```bash
$ psite predict -h
Usage: psite predict [OPTIONS] PATH_REF PATH_BAM PATH_MODEL PATH_OUT

  load pre-trained model and predict P-site offsets

  path_ref   : reference transcriptome (fasta) matching the bam
  path_bam   : alignments of RPFs to reference transcriptome
  path_model : path to save the fitted model
  path_out   : output path of bam with PS (for P-site) tag

Options:
  -i, --ignore_txversion    whether to ignore trasncript version in ".\d+"
                            format  [default: False]
  -l, --rlen_min INTEGER    lower bound for mapped mapped length
  -u, --rlen_max INTEGER    upper bound for mapped read length
  -c, --chunk_size INTEGER  chunk size for prediction batch  [default: 100000]
  -h, --help                Show this message and exit.  [default: False]
```

#### `pbam`
This module predict P-site for each alignment and keeps only the first nucleotide next to the P-site offset in the alignment. Thus, each aligment in the output contains a single site. The input can be either genomic alignments or transcriptomic alignments.

```bash
$ psite pbam -h
Usage: psite pbam [OPTIONS] PATH_REF PATH_BAM PATH_MODEL PATH_OUT

  generate bam with only P-site regions

  path_ref   : reference transcriptome (fasta) matching the bam
  path_bam   : alignments of RPFs to reference transcriptome
  path_model : path to save the fitted model
  path_out   : output path of bam with P-site regions only

Options:
  -f, --out_format [bam|sam]  P-site alignemnt output format  [default: bam]
  -i, --ignore_txversion      whether to ignore trasncript version in ".\d+"
                              format  [default: False]
  -l, --rlen_min INTEGER      lower bound for mapped mapped length
  -u, --rlen_max INTEGER      upper bound for mapped read length
  -c, --chunk_size INTEGER    chunk size for prediction batch  [default:
                              100000]
  -h, --help                  Show this message and exit.  [default: False]
```

#### `coverage`
This module calculates the genome or transcriptome-wide coverage of RPF P-sites using by employing the `PS` tag generatd by `predict` module.

```bash
$ psite coverage -h
Usage: psite coverage [OPTIONS] PATH_BAM PREFIX

  calculate the coverage for plus strand and minus strand seperately

  path_bam: sorted aligment bam file with the PS tag (for P-site offset)
  prefix  : output prefix of P-site coverage tracks in bigWig format

Options:
  -l, --rlen_min INTEGER  lower bound for RPF mapped length  [default: 25]
  -u, --rlen_max INTEGER  upper bound for mapped read length  [default: 40]
  -q, --mapq_min INTEGER  minimum mapping quality  [default: 10]
  -i, --ignore_supp       whether to ignore supplementary alignments
                          [default: False]
  -h, --help              Show this message and exit.  [default: False]
```

---------------------------------------

### An example workflow to use PSite

##### Prepare input files
After trimming adapters and optionally removing reads derived from rRNAs or tRNAs, map ribosomal footprints to the reference genome with STAR:

```bash
STAR --runThreadN 16 --outFilterType BySJout --outFilterMismatchNmax 2 --genomeDir genome_index --readFilesIn sample_RPF.fq.gz  --outFileNamePrefix sample_RPF --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --outFilterMultimapNmax 1 --outFilterMatchNmin 16 --alignEndsType EndToEnd --outSAMattributes NH HI AS nM NM MD
```

The parameter `--quantMode TranscriptomeSAM` will instruct STAR to translate the genomic alignments into transcript alignments, which will be used to train the GBT model. Since many uniquely mapped reads in genomic alignments will become multi-mapping reads in transcriptome alignment due to the presence of alternative transcript isoforms, `--outFilterMultimapNmax 1` parameter is included to exlude only multi-mapping reads in genomic alignments.

PSite need to know the position of annotated start codons and stop codons of all protein-coding transcripts, which can be obtained with the helper scirpts located in the `scripts` directory:
```bash
Rscript --vanilla scripts/extract_txinfo_ensembl.R gene_annotations.gtf txinfo.tsv
```

For a gene with multiple transcript isoforms, only a represent isoform is used in analysis. By default, PSite use the longest transcript. However, a more reasonable choice is the most abundant isoform. Therefore, if the information of transcript abundance as calculated by [kallisto](https://pachterlab.github.io/kallisto/) or [salmon](https://salmon.readthedocs.io/en/latest/salmon.html) is provided, PSite can automatically determine the most abundant transcript isoform for later use:
```bash
salmon quant -p4 --seqBias --gcBias --posBias -l A -i salmon_index -r sample_RNA.fq.gz -o salmon_results
```

##### Run PSite
The first step is to train a GBT model with `train` module with the transcriptome bam. The fitted model will be saved in [pickle](https://docs.python.org/3/library/pickle.html) format.

```bash
psite train -i -t salmon -e salmon_results/quant.sf \
    cdna.all.fa.gz sample_RPF.Aligned.toTranscriptome.out.bam output_prefix txinfo.tsv
```

Once the model is successfully trained, it can be used to predict P-site off sites for robosome footprints that are mapped to the reference genome or reference trancriptomes. It should be noted that if you use genome bam for prediction, genomic fasta sould be used as input, and vice versa.

Model training is slow for large dataset. `-f` parameter can be used to select only a subset of alignments for training. In our experience, this can significantly improve speed and reduce memory usage, while maintaining similar accuracy.

```bash
# with transcriptomic bam
psite predict -i all_transcripts.fa sample_RPF.Aligned.toTranscriptome.out.bam output_prefix.gbt.pickle sample_RPF.transcriptome.tag.bam

# with genomic bam
psite predict -i genome.fa sample_RPF.Aligned.sortedByCoord.out.bam output_prefix.gbt.pickle sample_RPF.genome.tag.bam
```
example output:

```
r1	16	2L	10716	255	29M	*	0	0	TACAATTTATTAAATGGGGACGGACCAAT	IIIIIIIIIIIIIIIIIIIIIHIIDDDDD	NH:i:1	HI:i:1	AS:i:28	nM:i:0	NM:i:0	MD:Z:29	jM:B:c,-1	jI:B:i,-1	PS:i:10
r2	16	2L	10836	255	33M	*	0	0	TGTCAACTTTTATCCTTTGTACCTTTCTACAAA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIDDDDD	NH:i:1	HI:i:1	AS:i:32	nM:i:0	NM:i:0	MD:Z:33	jM:B:c,-1	jI:B:i,-1	PS:i:12
r3	16	2L	10891	255	30M	*	0	0	CGGGTAAAGGGTATAAAGTCACTACGCGAA	GGD1HEC?HHIHHGF<1?IHHDHHF0@DD0	NH:i:1	HI:i:1	AS:i:29	nM:i:0	NM:i:0	MD:Z:30	jM:B:c,-1	jI:B:i,-1	PS:i:12
r4	16	2L	11027	255	31M	*	0	0	TTTCTGTTTGTATGTAAATCGCGTTTAATTT	IIIIIIIIIIIIIIIIIIIIIIIIIIDDDDD	NH:i:1	HI:i:1	AS:i:30	nM:i:0	NM:i:0	MD:Z:31	jM:B:c,-1	jI:B:i,-1	PS:i:12
r5	16	2L	11073	255	32M	*	0	0	CGTTCCTATTTTGCTGTCCCCGTTCGATTTTT	IHHIIIIIIIIIIIIIIHIIIIHIIIIDDDDD	NH:i:1	HI:i:1	AS:i:31	nM:i:0	NM:i:0	MD:Z:32	jM:B:c,-1	jI:B:i,-1	PS:i:12
r6	16	2L	11077	255	29M	*	0	0	CCTATTTTGCTGTCCCCGTTCGATTTTTA	@CC?FCD<<CG</H@HDHHHIHCF0?@@D	NH:i:1	HI:i:1	AS:i:28	nM:i:0	NM:i:0	MD:Z:29	jM:B:c,-1	jI:B:i,-1	PS:i:12
r7	16	2L	11132	255	31M	*	0	0	AAATTACATCAGGACTAGTACTCGTTTGCGT	IIIHIIHIIIIIHIHIIIHHIGIIIIDDDBD	NH:i:1	HI:i:1	AS:i:30	nM:i:0	NM:i:0	MD:Z:31	jM:B:c,-1	jI:B:i,-1	PS:i:11
r8	16	2L	11138	255	32M	*	0	0	CATCAGGACTAGTACTCGTTTGCGTCGTATTT	1HHHIIHIHGIHIIIIIIHHIIHGGHHDDD@@	NH:i:1	HI:i:1	AS:i:31	nM:i:0	NM:i:0	MD:Z:32	jM:B:c,-1	jI:B:i,-1	PS:i:12
r9	16	2L	11138	255	29M	*	0	0	CATCAGGACTAGTACTCGTTTGCGTCGTA	CFEGFHFCIHIHIIIHDIGIIHCHD<BB?	NH:i:1	HI:i:1	AS:i:28	nM:i:0	NM:i:0	MD:Z:29	jM:B:c,-1	jI:B:i,-1	PS:i:10
r10	16	2L	11140	255	32M	*	0	0	TCAGGACTAGTACTCGTTTGCGTCGTATTTCT	FCCHHCHIHIHHE0HHHDECHIH?IHG@0@D@	NH:i:1	HI:i:1	AS:i:31	nM:i:0	NM:i:0	MD:Z:32	jM:B:c,-1	jI:B:i,-1	PS:i:13
```

It also possible to output alignments with P-site locations only, which can be used for downstream applicaitons such as translated ORF prediction with [RibORF](https://github.com/zhejilab/RibORF).

```bash
psite pbam -f sam -p2 genome.fa sample_RPF.Aligned.sortedByCoord.out.bam output_prefix.gbt.pickle sample_RPF.genome.psite.sam
```

Here are a few lines from an example output:

```
r1      16      1       531180  255     1M      *       0       0       G       J       NH:i:1  HI:i:1  AS:i:30 nM:i:0  NM:i:0  MD:Z:31
r2      16      1       531180  255     1M      *       0       0       G       J       NH:i:1  HI:i:1  AS:i:30 nM:i:0  NM:i:0  MD:Z:31
r3      0       1       629921  255     1M      *       0       0       A       J       NH:i:1  HI:i:1  AS:i:31 nM:i:1  NM:i:1  MD:Z:0C33
r4      0       1       629921  255     1M      *       0       0       A       J       NH:i:1  HI:i:1  AS:i:31 nM:i:1  NM:i:1  MD:Z:0C33
r5      0       1       629922  255     1M      *       0       0       T       J       NH:i:1  HI:i:1  AS:i:32 nM:i:0  NM:i:0  MD:Z:33
r6      0       1       629922  255     1M      *       0       0       T       J       NH:i:1  HI:i:1  AS:i:29 nM:i:1  NM:i:1  MD:Z:0C31
r7      0       1       629922  255     1M      *       0       0       T       J       NH:i:1  HI:i:1  AS:i:29 nM:i:1  NM:i:1  MD:Z:0C31
r8      0       1       629922  255     1M      *       0       0       T       J       NH:i:1  HI:i:1  AS:i:29 nM:i:1  NM:i:1  MD:Z:0C31
r9      0       1       629922  255     1M      *       0       0       T       J       NH:i:1  HI:i:1  AS:i:32 nM:i:0  NM:i:0  MD:Z:33
r10     0       1       629922  255     1M      *       0       0       T       J       NH:i:1  HI:i:1  AS:i:30 nM:i:1  NM:i:1  MD:Z:0C32
```

PSite also has a module for fast calculation of genome or transcriptome P-site coverage of robosome footprints. The alignments should be sort be coordinates before coverage calculation.
```bash
# sort bam
samtools sort -@ 8 -O bam -o sample_RPF.genome.tag.sorted.bam sample_RPF.genome.tag.bam

# calculate coverage
psite coverage -q0 sample_RPF.genome.tag.sorted.bam sample_RPF.psite_cov
```

---------------------------------------

#### Other information
Please use the [issues](https://github.com/gxelab/psite/issues) panel for questions related to PSite, bug reports or feauture fequests.

If you find PSite helpful, please cite it with the zenodo doi:
> Hong Zhang, Yue Chang, Tianyu Lei, 2022, [10.5281/zenodo.7046270](https://doi.org/10.5281/zenodo.7046270).
