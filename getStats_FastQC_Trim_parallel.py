import os
from subprocess import call as unix
import glob
from joblib import Parallel, delayed
from multiprocessing import Process


taxa_trim = ['ANOCA', 'CAEBR', 'DANPL', 'DANRE', 'MNELE', 'NASVI', 'ONCVO', 'SCHMA', 'TRISP', 'XIPMA', 'DASNO', 'DROME', 'GASAC', 'PETMA', 'TETNG', 'OTOGA', 'BRAFL', 'MYOLU', 'NEMVE', 'RHOPR', 'TAEGU', 'TAKRU', 'LOXAF', 'TRICA']

SE_reads = ['MNELE', 'TRISP', 'CAEEL', 'NASVI', 'PETMA', 'TETNG', 'XIPMA', 'ANOCA', 'TAEGU', 'MELGA', 'CHICK', 'ORNAN', 'MONDO', 'MACEU', 'LOXAF', 'MOUSE', 'HUMAN', 'TRICA']
PE_reads = ['AMPQE', 'NEMVE', 'SCHMA', 'ONCVO', 'CAEBR', 'STRMM', 'DAPPU', 'ZOONE', 'RHOPR', 'DANPL', 'DROME', 'AEDAE', 'STRPU', 'BRAFL', 'LEPOC', 'DANRE', 'ASTMX', 'GADMO', 'GASAC', 'TAKRU', 'ORENI', 'ORYLA', 'POEFO', 'LATCH', 'XENTR', 'PELSI', 'ANAPL', 'DASNO', 'MYOLU', 'OTOGA']



def trim_data(OMA_ID):
    print(OMA_ID, ' is being processed')

    #Trim fastq files to remove adapter (TruSeq3 file) and based on quality and read length
    os.chdir("/data0/bspm/shortRead_rerun_mapping/" + OMA_ID)
    os.mkdir("data_fastq_trimmed")
    os.chdir("data_fastq")
    if OMA_ID in SE_reads:
        unix("cp ../../TruSeq3-SE.fa .", shell=True)
        single_IDs = []
        for single_file in glob.glob("*.gz"):
            SRR_ID = single_file[:-9]
            single_IDs.append(SRR_ID)
        for single_ID in set(single_IDs):
            unix("java -jar /home/bspm/bin/Trimmomatic-0.38/trimmomatic-0.38.jar SE -phred33 -trimlog " + single_ID + "_Logfile.txt " + single_ID + ".fastq.gz" + " " + single_ID + "_trim.fastq.gz " + "ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:10:30 MINLEN:35", shell=True)
    else:
        if OMA_ID in PE_reads:
            unix("cp ../../TruSeq3-PE.fa .", shell=True)
            IDs = []
            for paired_file in glob.glob("*.gz"):
                SRR_ID = paired_file[:-11]
                IDs.append(SRR_ID)
            for ID in set(IDs):
                unix("java -jar /home/bspm/bin/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 -trimlog " + ID + "_Logfile.txt " + ID + "_1.fastq.gz " + ID + "_2.fastq.gz " + ID + "_1_paired_trim.fastq.gz " + ID + "_1_unpaired_trim.fastq.gz " + ID + "_2_paired_trim.fastq.gz " + ID + "_2_unpaired_trim.fastq.gz " + "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:10:30 MINLEN:35", shell=True)
    unix("mv *trim* ../data_fastq_trimmed/", shell=True)
    unix("mv *Logfile* ../data_fastq_trimmed/", shell=True)
    unix("rm *TruSeq3*", shell=True)

    #Run FastQC on trimmed files
    os.chdir("/data0/bspm/shortRead_rerun_mapping/" + OMA_ID)
    os.mkdir("FastQC_postTrim")
    os.chdir("data_fastq_trimmed")
    unix("perl /home/bspm/bin/FastQC/fastqc -o ../FastQC_postTrim --noextract -q *.gz", shell=True)

    print(OMA_ID, ' is finished processed')



for k in taxa_trim:
    trim_data(k)
#Parallel(n_jobs=24)(delayed(trim_data)(k) for k in taxa_trim)
