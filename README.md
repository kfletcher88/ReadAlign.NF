# ReadAlign.NF
A nextflow pipeline to align reads to a reference, merge common sample IDs and call SNPs.\
Toy data has been provided using *E. coli* genome assembly [GCA_000005845.2](https://www.ncbi.nlm.nih.gov/assembly/GCF_000005845.2/). These are available under `./ref` and `./reads`.\
The pipeline requires a CSV file, as below. The file `Reads.csv` can be used to run the toy data.
```
SampleID,FCID,Read1,Read2
Ecoli_1,Sim1,Reads/EcoliSS_1.1.fq,Reads/EcoliSS_1.2.fq
Ecoli_2,Sim1,Reads/EcoliSS_2.1.fq,Reads/EcoliSS_2.2.fq
Ecoli_2,Sim2,Reads/EcoliSS_3.1.fq,Reads/EcoliSS_3.2.fq
```

A nextflow installation is required to run the pipeline.\
The pipeline is containerized using conda.\

To run the piepline please:
```
git clone https://github.com/kfletcher88/ReadAlign.NF.git
nextflow run ReadAlign.nf --reference Ref/GCA_000005845.2_ASM584v2_genomic.fna
```

Note that resume is automatic, specified in the `nextflow.config` file.\

Last tested 05/22/2023
