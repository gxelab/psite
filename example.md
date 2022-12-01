### Analysis of example datasets from Drosophila melanogaster

#### acquire datasets
```bash
# ln -s /nfs_data/leity/poj/fly/data/SRR19387520.fastq.gz embryo_mid_rpf.fq.gz
# ln -s /nfs_data/leity/poj/fly/data/SRR19387518.fastq.gz embryo_late_rna.fq.gz
# ln -s /nfs_data/leity/poj/fly/data/SRR19387519.fastq.gz embryo_late_rpf.fq.gz
# ln -s /nfs_data/leity/poj/fly/data/SRR19387512.fastq.gz embryo_early_rpf.fq.gz
# zcat /nfs_data/leity/poj/fly/patraquim_psite/RNA_seq/SRR11432001.fastq.gz /nfs_data/leity/poj/fly/patraquim_psite/RNA_seq/SRR11432002.fastq.gz | gzip -c > embryo_mid_rna.fq.gz
# ln -s /nfs_data/database/fly_riboseq/public/SRR3031135.fastq.gz s2_normal_rpf.fq.gz
# ln -s /nfs_data/database/fly_riboseq/public/SRR3031124.fastq.gz embryo_2h_rpf.fq.gz
# zcat /nfs_data/database/fly_riboseq/dmel/em_0_2h_mrna_lane1.fq.gz /nfs_data/database/fly_riboseq/dmel/em_0_2h_mrna_lane2.fq.gz | gzip -c >embryo_2h_rna.fq.gz
# zcat /nfs_data/database/fly_riboseq/dmel/s2data/s2_wt_mrna_lane1.fq.gz /nfs_data/database/fly_riboseq/dmel/s2data/s2_wt_mrna_lane2.fq.gz | gzip -c >s2_normal_rna.fq.gz

curl -o embryo_early_rpf.fq.gz https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR19387512
curl -o embryo_early_rna.fq.gz https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR19387513

curl -o embryo_late_rpf.fq.gz https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR19387519
curl -o embryo_late_rna.fq.gz https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR19387518

curl -o embryo_mid_rpf.fq.gz https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR19387520
curl -o SRR11432001.fq.gz https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR11432001
curl -o SRR11432002.fq.gz https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR11432002
zcat SRR11432001.fq.gz SRR11432002.fq.gz gzip -c > embryo_mid_rna.fq.gz


curl -o s2_normal_rna.fq.gz https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR3031122
curl -o s2_normal_rpf.fq.gz https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR3031135
curl -o embryo_2h_rna.fq.gz https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR5075630
curl -o embryo_2h_rpf.fq.gz https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR3031124
```

#### trim adaptors
```bash
cutadapt -a AGATCGGAAGAGCACACGTC -j8 --trim-n -m 18 -o embryo_mid_rna.trim.fq.gz embryo_mid_rna.fq.gz >embryo_mid_rna.trim.log
cutadapt -a TGGAATTCTCGGGTGCCAAGG -j8 --trim-n -m 18 -o s2_normal_rna.trim.fq.gz s2_normal_rna.fq.gz >s2_normal_rna.trim.log
cutadapt -a TGGAATTCTCGGGTGCCAAGG -j8 --trim-n -m 18 -o s2_normal_rpf.trim.fq.gz s2_normal_rpf.fq.gz >s2_normal_rpf.trim.log
cutadapt -a TGGAATTCTCGGGTGCCAAGG -j8 --trim-n -m 18 -o embryo_2h_rna.trim.fq.gz embryo_2h_rna.fq.gz >embryo_2h_rna.trim.log
cutadapt -a TGGAATTCTCGGGTGCCAAGG -j8 --trim-n -m 18 -o embryo_2h_rpf.trim.fq.gz embryo_2h_rpf.fq.gz >embryo_2h_rpf.trim.log
```

#### filter rRNAs and tRNAs
```r
library(Biostrings)
trna <- readDNAStringSet('dmel-all-tRNA-r6.47.fasta.gz')
rrna <- readDNAStringSet('dmel-all-miscRNA-r6.47.fasta.gz')
rrna <- rrna[grepl(' type=rRNA;', names(rrna))]
misc <- c(trna, rrna)
writeXStringSet(misc, 'rRNA_tRNA_combined.fa')
```

```bash
# for i in *.fq.gz; do echo "bowtie2 -p8 --local --un-gz ${i%%.*}.clean.fq.gz -x bt2_index_rtRNA -U $i > /dev/null 2>${i%%.*}.clean.log"; done
bowtie2 -p16 --local --un-gz embryo_2h_rna.clean.fq.gz -x bt2_index_rtRNA -U embryo_2h_rna.trim.fq.gz > /dev/null 2>embryo_2h_rna.clean.log
bowtie2 -p16 --local --un-gz embryo_2h_rpf.clean.fq.gz -x bt2_index_rtRNA -U embryo_2h_rpf.trim.fq.gz > /dev/null 2>embryo_2h_rpf.clean.log
bowtie2 -p16 --local --un-gz embryo_early_rna.clean.fq.gz -x bt2_index_rtRNA -U embryo_early_rna.fq.gz > /dev/null 2>embryo_early_rna.clean.log
bowtie2 -p16 --local --un-gz embryo_early_rpf.clean.fq.gz -x bt2_index_rtRNA -U embryo_early_rpf.fq.gz > /dev/null 2>embryo_early_rpf.clean.log
bowtie2 -p16 --local --un-gz embryo_late_rna.clean.fq.gz -x bt2_index_rtRNA -U embryo_late_rna.fq.gz > /dev/null 2>embryo_late_rna.clean.log
bowtie2 -p16 --local --un-gz embryo_late_rpf.clean.fq.gz -x bt2_index_rtRNA -U embryo_late_rpf.fq.gz > /dev/null 2>embryo_late_rpf.clean.log
bowtie2 -p16 --local --un-gz embryo_mid_rna.clean.fq.gz -x bt2_index_rtRNA -U embryo_mid_rna.trim.fq.gz > /dev/null 2>embryo_mid_rna.clean.log
bowtie2 -p16 --local --un-gz embryo_mid_rpf.clean.fq.gz -x bt2_index_rtRNA -U embryo_mid_rpf.fq.gz > /dev/null 2>embryo_mid_rpf.clean.log
bowtie2 -p16 --local --un-gz s2_normal_rna.clean.fq.gz -x bt2_index_rtRNA -U s2_normal_rna.trim.fq.gz > /dev/null 2>s2_normal_rna.clean.log
bowtie2 -p16 --local --un-gz s2_normal_rpf.clean.fq.gz -x bt2_index_rtRNA -U s2_normal_rpf.trim.fq.gz > /dev/null 2>s2_normal_rpf.clean.log
```

#### mapping
```bash
# for i in *_rna.clean.fq.gz; do echo "STAR --outFilterType BySJout --runThreadN 8 --outFilterMismatchNmax 2 --genomeDir /nfs_data/database/ref_genomes/Dmel_em52/STAR_genome --readFilesIn $i --readFilesCommand zcat --outFileNamePrefix ${i%%.clean.fq.gz}_ --outSAMattributes All --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --outFilterMultimapNmax 1 --outFilterMatchNmin 16 --alignEndsType EndToEnd"; done
STAR --outFilterType BySJout --runThreadN 8 --outFilterMismatchNmax 2 --genomeDir /nfs_data/database/ref_genomes/Dmel_em52/STAR_genome --readFilesIn embryo_2h_rna.clean.fq.gz --readFilesCommand zcat --outFileNamePrefix embryo_2h_rna_ --outSAMattributes All --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --outFilterMultimapNmax 1 --outFilterMatchNmin 16 --alignEndsType EndToEnd
STAR --outFilterType BySJout --runThreadN 8 --outFilterMismatchNmax 2 --genomeDir /nfs_data/database/ref_genomes/Dmel_em52/STAR_genome --readFilesIn embryo_2h_rpf.clean.fq.gz --readFilesCommand zcat --outFileNamePrefix embryo_2h_rpf_ --outSAMattributes All --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --outFilterMultimapNmax 1 --outFilterMatchNmin 16 --alignEndsType EndToEnd
STAR --outFilterType BySJout --runThreadN 8 --outFilterMismatchNmax 2 --genomeDir /nfs_data/database/ref_genomes/Dmel_em52/STAR_genome --readFilesIn embryo_early_rna.clean.fq.gz --readFilesCommand zcat --outFileNamePrefix embryo_early_rna_ --outSAMattributes All --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --outFilterMultimapNmax 1 --outFilterMatchNmin 16 --alignEndsType EndToEnd
STAR --outFilterType BySJout --runThreadN 8 --outFilterMismatchNmax 2 --genomeDir /nfs_data/database/ref_genomes/Dmel_em52/STAR_genome --readFilesIn embryo_early_rpf.clean.fq.gz --readFilesCommand zcat --outFileNamePrefix embryo_early_rpf_ --outSAMattributes All --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --outFilterMultimapNmax 1 --outFilterMatchNmin 16 --alignEndsType EndToEnd
STAR --outFilterType BySJout --runThreadN 8 --outFilterMismatchNmax 2 --genomeDir /nfs_data/database/ref_genomes/Dmel_em52/STAR_genome --readFilesIn embryo_late_rna.clean.fq.gz --readFilesCommand zcat --outFileNamePrefix embryo_late_rna_ --outSAMattributes All --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --outFilterMultimapNmax 1 --outFilterMatchNmin 16 --alignEndsType EndToEnd
STAR --outFilterType BySJout --runThreadN 8 --outFilterMismatchNmax 2 --genomeDir /nfs_data/database/ref_genomes/Dmel_em52/STAR_genome --readFilesIn embryo_late_rpf.clean.fq.gz --readFilesCommand zcat --outFileNamePrefix embryo_late_rpf_ --outSAMattributes All --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --outFilterMultimapNmax 1 --outFilterMatchNmin 16 --alignEndsType EndToEnd
STAR --outFilterType BySJout --runThreadN 8 --outFilterMismatchNmax 2 --genomeDir /nfs_data/database/ref_genomes/Dmel_em52/STAR_genome --readFilesIn embryo_mid_rna.clean.fq.gz --readFilesCommand zcat --outFileNamePrefix embryo_mid_rna_ --outSAMattributes All --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --outFilterMultimapNmax 1 --outFilterMatchNmin 16 --alignEndsType EndToEnd
STAR --outFilterType BySJout --runThreadN 8 --outFilterMismatchNmax 2 --genomeDir /nfs_data/database/ref_genomes/Dmel_em52/STAR_genome --readFilesIn embryo_mid_rpf.clean.fq.gz --readFilesCommand zcat --outFileNamePrefix embryo_mid_rpf_ --outSAMattributes All --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --outFilterMultimapNmax 1 --outFilterMatchNmin 16 --alignEndsType EndToEnd
STAR --outFilterType BySJout --runThreadN 8 --outFilterMismatchNmax 2 --genomeDir /nfs_data/database/ref_genomes/Dmel_em52/STAR_genome --readFilesIn s2_normal_rna.clean.fq.gz --readFilesCommand zcat --outFileNamePrefix s2_normal_rna_ --outSAMattributes All --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --outFilterMultimapNmax 1 --outFilterMatchNmin 16 --alignEndsType EndToEnd
STAR --outFilterType BySJout --runThreadN 8 --outFilterMismatchNmax 2 --genomeDir /nfs_data/database/ref_genomes/Dmel_em52/STAR_genome --readFilesIn s2_normal_rpf.clean.fq.gz --readFilesCommand zcat --outFileNamePrefix s2_normal_rpf_ --outSAMattributes All --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --outFilterMultimapNmax 1 --outFilterMatchNmin 16 --alignEndsType EndToEnd
```

#### salmon quant
```bash
# salmon quant
# for i in *_rna.clean.fq.gz; do echo "salmon quant -p8 --seqBias --gcBias --posBias -l A -i /nfs_data/database/ref_genomes/Dmel_em52/salmon_cdna_ncrna -r $i -o ${i%%.clean.fq.gz}"; done
salmon quant -p8 --seqBias --gcBias --posBias -l A -i /nfs_data/database/ref_genomes/Dmel_em52/salmon_cdna_ncrna -r embryo_2h_rna.clean.fq.gz -o embryo_2h_rna
salmon quant -p8 --seqBias --gcBias --posBias -l A -i /nfs_data/database/ref_genomes/Dmel_em52/salmon_cdna_ncrna -r embryo_early_rna.clean.fq.gz -o embryo_early_rna
salmon quant -p8 --seqBias --gcBias --posBias -l A -i /nfs_data/database/ref_genomes/Dmel_em52/salmon_cdna_ncrna -r embryo_late_rna.clean.fq.gz -o embryo_late_rna
salmon quant -p8 --seqBias --gcBias --posBias -l A -i /nfs_data/database/ref_genomes/Dmel_em52/salmon_cdna_ncrna -r embryo_mid_rna.clean.fq.gz -o embryo_mid_rna
salmon quant -p8 --seqBias --gcBias --posBias -l A -i /nfs_data/database/ref_genomes/Dmel_em52/salmon_cdna_ncrna -r s2_normal_rna.clean.fq.gz -o s2_normal_rna
```
