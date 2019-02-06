import os
from subprocess import call as unix
import glob
from joblib import Parallel, delayed

#List of RNA seq datasets for each taxa
taxa_map = ['ANOCA', 'CAEBR', 'DANPL', 'DANRE', 'MNELE', 'NASVI', 'ONCVO', 'SCHMA', 'TRISP', 'XIPMA', 'DASNO', 'DROME', 'GASAC', 'PETMA', 'TETNG', 'OTOGA', 'BRAFL', 'MYOLU', 'NEMVE', 'RHOPR', 'TAEGU', 'TAKRU', 'LOXAF', 'TRICA']

#List of taxa based on whether they are paired or unpaired reads
SE_reads = ['MNELE', 'TRISP', 'CAEEL', 'NASVI', 'PETMA', 'TETNG', 'XIPMA', 'ANOCA', 'TAEGU', 'MELGA', 'CHICK', 'ORNAN', 'MONDO', 'MACEU', 'LOXAF', 'MOUSE', 'HUMAN']
PE_reads = ['AMPQE', 'NEMVE', 'SCHMA', 'ONCVO', 'CAEBR', 'STRMM', 'DAPPU', 'ZOONE', 'RHOPR', 'TRICA', 'DANPL', 'DROME', 'AEDAE', 'STRPU', 'BRAFL', 'LEPOC', 'DANRE', 'ASTMX', 'GADMO', 'GASAC', 'TAKRU', 'ORENI', 'ORYLA', 'POEFO', 'LATCH', 'XENTR', 'PELSI', 'ANAPL', 'DASNO', 'MYOLU', 'OTOGA']

new_dir = "/data2/bspm/shortRead_mapping/"
old_dir = "/data0/bspm/shortRead_rerun_mapping/"

#Build databases for each gene to be mapped to (found in fusion_nuc_fams dir)
def bowtie_build(OMA):
    os.chdir(new_dir)
    os.mkdir(OMA)
    os.chdir(OMA)
    os.mkdir("bowtie_index")
    os.chdir(old_dir + OMA)
    os.chdir("fusion_nuc_fams")
    for fasta_file in glob.glob("*.fasta"):
        fam_geneID = fasta_file[:-18]
        unix("bowtie2-build -f " + fasta_file + " " + new_dir + OMA + "/bowtie_index/" + fam_geneID, shell=True)

#If taxa has single end reads or paired end reads, run mapping with specified commands in bowtie
def bowtie_map(OMA):
    os.chdir(new_dir + OMA)
    os.mkdir("bowtie_mapping")
    os.chdir("bowtie_index")
    if OMA in SE_reads:
        single_IDs = []
        for single_file in glob.glob("*.bt2"):
            SRR_ID = single_file.split(".")[0]
            single_IDs.append(SRR_ID)
        for single_ID in set(single_IDs):
            for read_file in glob.glob(old_dir + OMA + "/data_fastq_trimmed/*.gz"):
                read = read_file.split("/")[-1]
                read_ID = read[:-14]
                unix("bowtie2 -x " + single_ID + " -U " + read_file + " | samtools view -S -h -F4 - > " + single_ID + "_" + read_ID + "_mapped.sam", shell=True)

    else:
        if OMA in PE_reads:
            IDs = []
            for paired_file in glob.glob("*.bt2"):
                SRR_ID = paired_file.split(".")[0]
                IDs.append(SRR_ID)
            read_files = []
            for ID in set(IDs):
                for read_file in glob.glob(old_dir + OMA + "/data_fastq_trimmed/*.gz"):
                    read = read_file.split("/")[-1]
                    read_ID = read[:-14]
                    read_SRR = read_ID[:-9]
                    trim_unpair = read_ID[-9:]
                    if trim_unpair == '_unpaired':
                        continue
                    else:
                        read_files.append(read_SRR)
                for reads in set(read_files):
                    unix("bowtie2 -x " + ID + " -1 " + old_dir + OMA + "/data_fastq_trimmed/" + reads + "_1_paired_trim.fastq.gz " + "-2 " + old_dir + OMA + "/data_fastq_trimmed/" + reads + "_2_paired_trim.fastq.gz | samtools view -S -h -F4 - > " + ID + "_" + reads + "_mapped.sam", shell=True)
    #Move output SAM files to bowtie_mapping directory
    unix("mv *.sam ../bowtie_mapping/", shell=True)

    
#Convert SAM outputs from bowtie to sorted BAM files
def sam_to_bam(OMA):
    os.chdir(new_dir + OMA)
    os.mkdir("bedtool_coverage")
    os.chdir("bedtool_coverage")
    os.mkdir("bam_files")
    os.mkdir("genCov_files")
    os.chdir(new_dir + OMA + "/bowtie_mapping/")
    for sam_file in glob.glob("*.sam"):
        map_ID = sam_file[:-11]
        unix("samtools view -b " + sam_file + " | samtools sort -o " + map_ID + "_sorted.bam", shell=True)
    unix("mv *.bam ../bedtool_coverage/bam_files", shell=True)

    
#Get positions on the gene where the reads mapped using bedtools genomecov on the bam files
def bedtools_gcov(OMA):
    os.chdir(new_dir + OMA + "/bedtool_coverage/bam_files")
    for bam_file in glob.glob("*.bam"):
        bamID = bam_file[:-11]
        unix("bedtools genomecov -ibam " + bam_file + " -d > " + "positionsCov_" + bamID + ".txt", shell=True)
    unix("mv *.txt ../genCov_files", shell=True)


#Run each command
for ID in taxa_map:
    bowtie_build(ID)
    bowtie_map(ID)
    sam_to_bam(ID)
    bedtools_gcov(ID)


#Run functions for each dataset in parallel on Tomoko
#Parallel(n_jobs=23)(delayed(bowtie_build)(ID) for ID in taxa_map)
#Parallel(n_jobs=23)(delayed(bowtie_map)(ID) for ID in taxa_map)
#Parallel(n_jobs=23)(delayed(sam_to_bam)(ID) for ID in taxa_map)
#Parallel(n_jobs=23)(delayed(bedtools_gcov)(ID) for ID in taxa_map)
