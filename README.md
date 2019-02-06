# mapping_scripts

## Scripts used in the full pipeline of mapping RNA seq reads to genes.

Scripts added include programs to 

(1) download sra files from ncbi ftp site 
(2) convert sra files to fastq 
(3) get information about the quality of sequence reads using FastQC 
(4) trim reads for adapters and quality and length of the reads 
(5) map the reads to the genes using bowtie2
(6) convert SAM output files to sorted BAM files
(7) show the mapping coverage on the genes with bedtools
