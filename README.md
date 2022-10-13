### psite
[![DOI](https://zenodo.org/badge/474568909.svg)](https://zenodo.org/badge/latestdoi/474568909)

Model-based inference of P-site offsets for RPFs

#### Dependency
- `numpy` >= 1.21.2
- `pandas` >= 1.3.4
- `biopython` >= 1.79 (Processing fasta)
- `sklearn` >= 1.1.1 (RandomForest model)
- `pysam` >= 0.17.0 (BAM parsing)
- `click` >= 8.1.2 (Commandline arguments)
- `pyBigWig` >= 0.3.18 (Writing bigWig)

#### Installation
This package is still under development and can be installed and uninstalled by
```bash
# installation
pip install -e .

# uninstall
pip uninstall piste
```

#### Usage
##### Prepare input files
1. Map RPF reads to genome with STAR:
```bash
STAR --runThreadN 16 --outFilterType BySJout --outFilterMismatchNmax 2 --genomeDir genome_index --readFilesIn sample_RPF.fq.gz  --outFileNamePrefix sample_RPF --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --outFilterMultimapNmax 1 --outFilterMatchNmin 16 --alignEndsType EndToEnd --outSAMattributes NH HI AS nM NM MD
```

2. Obtain transcript information from GTF:
```bash
Rscript --vanilla scripts/extract_txinfo_ensembl.R genome.gtf txinfo.tsv
```

3. Calculate transcript expression levels to determine the dominant transcript for each gene (This step is not necessary. In that case, the longest transcript of each gene will be used in analysuis):
```bash
salmon quant -p4 --seqBias --gcBias --posBias -l A -i salmon_index -r sample_RNA.fq.gz -o salmon_results
```

##### Run `psite`
Model training:
```bash
psite train -n2 -i -t salmon -e salmon_results/quant.sf \
    cdna.all.fa.gz sample_RPF.Aligned.toTranscriptome.out.bam fitted_model.pkl txinfo.tsv
```

Predict RPF P-site offsets (offset is stored in `PS` tag) 
```bash
# with transcriptome bam
psite predict -n2 -i all_transcripts.fa  sample_RPF.Aligned.toTranscriptome.out.bam fitted_model.pkl sample_RPF.transcriptome.psite.bam

# with genome bam
psite predict -n2 -i genome.fa  sample_RPF.Aligned.sortedByCoord.out.bam fitted_model.pkl sample_RPF.genome.psite.bam
```

Calculate genome-wide P-site coverage
```bash
# sort bam
samtools sort -@ 8 -O bam -o sample_RPF.genome.psite.sorted.bam sample_RPF.genome.psite.bam

# calculate coverage
psite coverage -q0 sample_RPF.genome.psite.sorted.bam sample_RPF.psite_cov
```

#### Modules
