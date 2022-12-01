### Analysis of example datasets from Drosophila melanogaster

##### acquire datasets
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

