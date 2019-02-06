# mapping_scripts

## Scripts used in the full pipeline of mapping RNA seq reads to genes.

Scripts added include programs to: 

1. download sra files from ncbi ftp site 
2. convert sra files to fastq 
3. get information about the quality of sequence reads using FastQC 
4. trim reads for adapters and quality and length of the reads 
5. map the reads to the genes using bowtie2
6. convert SAM output files to sorted BAM files
7. show the mapping coverage on the genes with bedtools

---

1. Download SRA files and convert to fastq

`download_fastq.py`

2. Run FastQC on raw reads, trim adapters and reads lower than phred33 and length 35 nucleotides

`getStats_FastQC_Trim_parallel.py`

3. Create bowtie index, map reads (PE or SE), convert output SAM to sorted BAM, get coverage information using bedtools genomecov

`bowtie_run.py`

4. Parse the genomecov output files to find whether the fusion breakpoint is covered by RNA reads

`bedtool_parse_tom.py`


---
