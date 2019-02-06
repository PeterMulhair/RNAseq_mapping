import os
from subprocess import call as unix
import glob
from joblib import Parallel, delayed
import queue

#taxa_list = ['DANPL', 'CAEBR', 'DROME', 'LOXAF', 'MNELE', 'NASVI', 'TRISP', 'XIPMA', 'TETNG', 'TAEGU', 'SCHMA', 'RHOPR', 'OTOGA', 'ONCVO', 'DANRE', 'BRAFL']
taxa_list = ['ANOCA', 'BRAFL', 'DANPL', 'DASNO', 'GASAC', 'MNELE', 'NASVI', 'ONCVO', 'PETMA', 'SCHMA', 'TAKRU', 'TETNG', 'TRISP', 'CAEBR', 'DANRE', 'DROME', 'LOXAF', 'MYOLU', 'NEMVE', 'OTOGA', 'RHOPR', 'TAEGU', 'XIPMA', 'DANRE']

#Create dictionary from dico file where key is gene ID and value is species ID
id2Speices=dict()
with open('/data1/bspm/compSearch_output/blast_allGenomes_Ray.out.cleanNetwork.dico', 'r') as f2:
    for line in f2:
        sp=line.strip().split('\t')[0]
        id=line.strip().split('\t')[1]
        id2Speices[id]=sp

#Create dictionary of geneID to breakpoint gap length                                                                                       
fusionGene_bpGap_len = dict()
with open("/data1/bspm/shortRead_mapping/raw/unambig_dicts/Fam2DomainRegions_allcomps_unambig.csv") as csvfile:
    next(csvfile)
    for line in csvfile:
        domain_regions = []
        lines = line.split('"')
        IDs = lines[0]
        famID = IDs.split(',')[0]
        geneID = IDs.split(',')[1]
        for item in lines:
            if '(' in item:
                domain_regions.append(item)
        domain1 = domain_regions[0].split(',')
        bp_start = domain1[-1].strip(')').strip(' ')
        nuc_bp_start = int(bp_start)*3

        end_domain = domain_regions[-1].split(',')
        bp_end = end_domain[0].strip('(').strip(' ')
        nuc_bp_end = int(bp_end)*3

        bp_start_end = (nuc_bp_start, nuc_bp_end)

        bp_gap_len = int(nuc_bp_end) - int(nuc_bp_start)
        fusionGene_bpGap_len[geneID] = bp_start_end

#Create dict of taxaID to breakpoint gap length                                                                                             
fusionGene_bpGap_len_sp = dict()
for geneID, bp_start_end in fusionGene_bpGap_len.items():
    if geneID[1:] in id2Speices.keys():
        spID = id2Speices[geneID[1:]]
        fusionGene_bpGap_len_sp[spID] = bp_start_end


for taxa in taxa_list:
    print(taxa)
    os.chdir(taxa + "/bedtool_coverage/genCov_files")
    with open("Confirmed_CompFams_01.txt", "w") as outF_01, open("Confirmed_CompFams.txt", "w") as outF:
        for my_file in glob.glob("*positionsCov*"):
            geneID = my_file.split('_')[2]
            famID = my_file.split('_')[1]
            readID = my_file.split('_')[3]
            bp = fusionGene_bpGap_len_sp[geneID]
            bp_start = bp[0]
            bp_start_extended = bp_start - 10 
            bp_end = bp[1]
            bp_end_extended = bp_end + 10
            covered = 'FALSE'
            #Create queue array to save to memory
            q = queue.Queue()
            bp_ext_len = 0
            with open(famID + "_" + geneID + "_breakpoint_cov.txt", "w") as outF1:
                with open(my_file) as f:
                    for line in f:
                        lines = line.split('\t')
                        cov = lines[-1].strip()
                        nuc_pos = lines[1]
                        #If the nucleotide position is between the breakpoints
                        if (int(nuc_pos) >= bp_start) and (int(nuc_pos) <= bp_end):
                            if int(cov) == 0:
                                covered = 'FALSE'
                                break
                            else:
                                q.put(line)
                                covered = 'TRUE'

                    print(famID, geneID, readID, covered)
                    if q.qsize() == (bp_end - bp_start + 1):
                        outF.write(famID + "\t" + geneID + "\t" + readID + "\n")
                        while q.empty() == False:
                            outF1.write(q.get())


                with open(my_file) as f:
                    for line in f:
                        lines = line.split('\t')
                        cov = lines[-1].strip()
                        nuc_pos = lines[1]
                        if (int(nuc_pos) >= bp_start_extended) and (int(nuc_pos) <= bp_end_extended):
                            if int(cov) == 0:
                                covered = 'FALSE'
                                break
                            else:
                                bp_ext_len += 1
                                
                    if bp_ext_len == (bp_end_extended - bp_start_extended + 1):
                        outF_01.write(famID + "\t" + geneID + "\t" + readID + "\n")

    unix("find -iname '*breakpoint_cov.txt' -type f -empty -delete", shell=True)
    os.chdir("../../../")
    
#                if covered == 'TRUE':
#                    outF.write(famID + "\t" + geneID + "\t" + readID + "\n")
#    os.chdir("../../../")
