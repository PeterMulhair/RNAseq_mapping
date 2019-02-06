import os
from subprocess import call as unix
import glob
from joblib import Parallel, delayed
from multiprocessing import Process


#Dictionary of taxa to the ncbi bioproject ID where the raw reads are found
taxa_SRP_dict = {'AMPQE': ['SRP044247'], 'MNELE': ['SRP014138'], 'NEMVE': ['SRP113508'], 'SCHMA': ['ERP000427'], 'TRISP': ['SRP014138'], 'ONCVO': ['ERP001350'], 'CAEBR': ['SRP011366'], 'STRMM': ['SRP041623'], 'RHOPR': ['SRP057515'], 'NASVI': ['SRP029983'], 'DANPL': ['SRP015992'], 'DROME': ['SRP002072'], 'BRAFL': ['SRP056868'], 'PETMA': ['SRP009480'], 'LEPOC': ['SRP042013'], 'DANRE': ['ERP000263'], 'ASTMX': ['SRP058866'], 'GASAC': ['SRP012923'], 'TETNG': ['ERP011338'], 'TAKRU': ['SRP030658'], 'ORYLA': ['SRP032993'], 'ANOCA': ['SRP009813'], 'PELSI': ['SRP119729'], 'TAEGU': ['SRP063457'], 'DASNO': ['SRP012922'], 'MYOLU': ['SRP055976'], 'OTOGA': ['ERP014610']}



def download_data(OMA_ID, SRA_ID):
        print(OMA_ID, ' is being processed')
        dir = os.mkdir(OMA_ID)
        os.chdir("/data1/bspm/shortRead_mapping/" + OMA_ID)
        #Download the sra directories for each taxa bioproject
        for ID in SRA_ID:
                first_ID = ID[:3]
                second_ID = ID[:6]
                final_ID = ID
                get_SRA = "ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByStudy/sra/" + first_ID + "/" + second_ID + "/" + final_ID + "/"
                unix("wget -r " + get_SRA, shell=True) 

	##Put all the .sra read files into data_fastq dir, and remove ftp folder
        os.chdir("/data1/bspm/shortRead_mapping/" + OMA_ID)
        for ID in SRA_ID:
                first_ID = ID[:3]
                second_ID = ID[:6]
                final_ID = ID
                if first_ID == "SRP":
                        read_ID = "SRR"
                else:
                        read_ID = "ERR"
                SRA_files = "ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByStudy/sra/" + first_ID + "/" + second_ID + "/" + final_ID + "/"
                unix("mv " + SRA_files + "*" + read_ID + "*" + "/*.sra " + "/data1/bspm/shortRead_mapping/" + OMA_ID, shell=True)

        os.chdir("/data1/bspm/shortRead_mapping/" + OMA_ID)
        unix("rm -r ftp-trace.ncbi.nih.gov", shell=True)
        unix("mkdir data_sra", shell = True)
        unix("mv *.sra data_sra", shell = True)
        os.chdir("/data1/bspm/shortRead_mapping/")

	##Get fastq file for the SRA files using fastq-dump
        os.chdir("/data1/bspm/shortRead_mapping/" + OMA_ID)
        os.mkdir("data_fastq")
        os.chdir("data_sra")
        for sra_file in glob.glob("*.sra"):
                unix("fastq-dump --gzip --skip-technical --readids --dumpbase --split-files " + sra_file, shell=True)
        unix("mv *.gz ../data_fastq", shell=True)
        unix("rm *.sra", shell=True)
        os.chdir("../")
        unix("rm -r data_sra", shell=True)
        os.chdir("/data1/bspm/shortRead_mapping/")
        print(OMA_ID, ' is finished processing')



for k, v in taxa_SRP_dict.items():
        download_data(k, v)
        
#Parallel(n_jobs=36)(delayed(download_data)(k, v) for k, v in taxa_SRP_dict.items())
