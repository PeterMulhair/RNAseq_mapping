# RNAseq mapping scripts

## Scripts used in the full pipeline of mapping RNA seq reads to genes

Scripts added include programs to: 

* download sra files from ncbi ftp site 
* convert sra files to fastq 
* get information about the quality of sequence reads using FastQC 
* trim reads for adapters and quality and length of the reads 
* map the reads to the genes using bowtie2
* convert SAM output files to sorted BAM files
* show the mapping coverage on the genes with bedtools

---

1. Download SRA files and convert to fastq

`download_fastq.py`

2. Run FastQC on raw reads, trim adapters and reads lower than phred33 and length 35 nucleotides

`getStats_FastQC_Trim_parallel.py`

3. Create bowtie index, map reads (PE or SE), convert output SAM to sorted BAM, get coverage information using bedtools genomecov

`bowtie_run.py`

4. Parse the genomecov output files to find whether the fusion breakpoint is covered by RNA reads

`bedtool_parse.py`


---
