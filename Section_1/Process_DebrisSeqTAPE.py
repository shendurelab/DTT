# 2024.03.09
# This script takes in:
#   1) input fastq.gz file (R1 or assembled with PEAR, if paired-end)
#   2) output-file-name
#   3) threshold for calling real reads (a fraction of the total reads from the sequencing)

# In addition, for mGASv3 and mGASv5, we read-in the TargetBC previously found from cell-line selection

#   for now, let's run it like:

#python Step1_DebrisSeqTAPE_collect_v1.py \
#    --input_file 'mGASv4sci-P4-C10*fastq.gz' \
#    --output_file 'Step1File_DebrisSeqTAPE_mGASv5.csv' \
#    --Read_threshold 0.002

# But run this through a loop for detecting all the FASTQ files, like:


# for i in `ls -1 mGASv4sci-P*_R1_001.fastq.gz`
# do
#     Name=$(echo $i | cut -d'.' -f 1)
#     echo $i
#     echo ${Name}
#     python Step1_DebrisSeqTAPE_collect_v1.py \
#         --input_file $i \
#         --output_file 'Step1File_DebrisSeqTAPE_mGASv4.csv' \
#         --Read_threshold 0.002
# done

# 

import argparse
import gzip
import subprocess
import os,sys,csv,re
import itertools
import numpy as np
from collections import Counter
from collections import OrderedDict
from operator import itemgetter
from datetime import datetime
from optparse import OptionParser,OptionGroup

#python TAPE_gDNA_process_v1.py \
#    --input_file 'mGASv5-Debris-P1-A1.assembled.fastq.gz' \
#    --output_file 'mGASv5-Debris-P1-A1.csv' \
#    --Read_threshold 200


# Additional function
def is_within_edit_distance(sequence1, sequence2):
    if len(sequence1) != len(sequence2):
        return False

    num_mismatches = 0
    for i in range(len(sequence1)):
        if sequence1[i] != sequence2[i]:
            num_mismatches += 1
            if num_mismatches > 1:
                return False

    return True


def find_sequences_within_edit_distance(given_sequence, sequence_set):
    within_distance = []
    for sequence in sequence_set:
        if is_within_edit_distance(given_sequence, sequence):
            within_distance.append(sequence)

    return within_distance



# For mGASv3 and mGASv5, we have a set of TargetBC we are using (previously found)
# c = 0
# TargetBC_C5 = []
# with open('mESC_Clone5_TapeBC_recovery.csv') as f:
#     reader = csv.reader(f)
#     for row in reader:
#         if c == 0:
#             c += 1
#             continue
#         if c > 0:
#             TargetBC = row[0]
#             DNA = int(row[1])
#             if DNA > 50000:
#                 TargetBC_C5.append(TargetBC)
                
#print(TargetBC_C5)
# print('TargetBC has been loaded.')



def FASTQ_TAPE_reading(read_file):

    # ===========================================================
    ## Reading in the forward PM16 TAPE reads to make one table
    # ===========================================================
    #print ("Reading fwd_reads")

    # Overall reads
    Total_reads = 0
    Seq_table = {}
        
    print ("===== Reading reads for "+ read_file +" =====  \n")
        
        
    cmd = "zcat "+read_file
    
    with os.popen(cmd) as handle: 
        handle = iter(handle)
        for line in handle:
            if line[0] == '@':
                Total_reads +=1
                
                try:
                    seq = next(handle).rstrip('\n')
                except StopIteration:
                    print(seq)
                    
                corrected_barcode=seq[29:41]+','+seq[58:65]+','+seq[78:85]+','+seq[98:105]+','+seq[118:125]+','+seq[138:145]+','+seq[158:165]
                try:
                    Seq_table[corrected_barcode] += 1
                except:
                    Seq_table[corrected_barcode] = 1
                    
    handle.close()
    print(Total_reads)
    
    return Seq_table


def TAPE_correction(plate_loc, Seq_table):

    # ===========================================================
    ## Correcting TargetBC and processing sequential insertion
    # ===========================================================

    # Overall reads
    Total_reads = 0
    Seq_table2 = {}
    Seq_table3 = {}
    TAPE_table = {}
    
    
    for k,nRead in Seq_table.items():
        row = k.split(',')
        TargetBC_correction = [row[0]]
        #TargetBC_correction = find_sequences_within_edit_distance(TargetBC, TargetBC_C5)
        if len(TargetBC_correction) == 1:
            if row[1][3:] == 'GGAT':
                if row[2][3:] == 'GGAT':
                    if row[3][3:] == 'GGAT':
                        if row[4][3:] == 'GGAT':
                            if row[5][3:] == 'GGAT':
                                if row[6][3:] == 'GGAT':
                                    TAPE = ','.join(row[1:7])
                                else:
                                    TAPE = ','.join(row[1:6])+',None'
                            else:
                                TAPE = ','.join(row[1:5])+',None,None'
                        else:
                            TAPE = ','.join(row[1:4])+',None,None,None'
                    else:
                        TAPE = ','.join(row[1:3])+',None,None,None,None'
                else:
                    TAPE = ','.join(row[1:2])+',None,None,None,None,None'

                TAPE = plate_loc+','+TargetBC_correction[0]+'-'+row[1][:3]+','+TAPE
                Total_reads += nRead
                row = TAPE.split(',')


            else:
                TAPE = plate_loc+','+TargetBC_correction[0]+'-Non'+',None,None,None,None,None,None'
                Total_reads += nRead
                row = TAPE.split(',')

            if len(row) == 8:
                try:
                    TAPE_table[TAPE] += nRead
                except:
                    TAPE_table[TAPE] = nRead

    TAPE_table = dict(sorted(TAPE_table.items(), key=lambda item: item[1], reverse=True))
    
    return TAPE_table, Total_reads



def TAPE_MostFreq(TAPE_table):
    TAPE_table2 = {}
    read_table = {}
    
    for k,v in TAPE_table.items():
        row = k.split(',')
        cell_target = ','.join(row[0:3]) # Plate location + TargetBC + Site1
        #rest = ','.join(row[2:])+','+str(v)
        try:
            if read_table[cell_target] < v:
                read_table[cell_target] = v
                TAPE_table2[k] = v
        except KeyError:
            read_table[cell_target] = v
            TAPE_table2[k] = v

    TAPE_table2 = dict(sorted(TAPE_table2.items(), key=lambda item: item[1], reverse=True))
    return TAPE_table2


if __name__ == '__main__':

    parser = argparse.ArgumentParser('Script to generate a file with counts of all TargetBC-TAPE barcodes found within each cell in a 10X experiment.')
    parser.add_argument('--input_file', '-i', help='Fastq.gz file with 6xTAPE reads.')
    parser.add_argument('--output_file', '-o', help='Tab delimited file with cell, mutation barcode, read count, umi count. All observed barcodes correctable to a whitelist are reported.')
    parser.add_argument('--Read_threshold', type=float, default=2, help='Threshold for the read count.')

    args = parser.parse_args()

    read_file=args.input_file
    Read_threshold=args.Read_threshold
    output_file_name=args.output_file

    plate_loc = read_file[10:16]+'...........' # Specific for mGASv4
    Seq_table = FASTQ_TAPE_reading(read_file)
    TAPE_table, Total_reads = TAPE_correction(plate_loc, Seq_table) #, TargetBC_C5) # Specific for mGASv4
    print(Total_reads, Read_threshold)
    read_cutoff = Total_reads * Read_threshold
    print("Read cut-off was: " + str(read_cutoff))
    TAPE_table2 = TAPE_MostFreq(TAPE_table)


    output_file=open(output_file_name, 'a')
    original_stdout = sys.stdout
    #output_file.write(','.join(['name','TargetBC','Site1','Site2','Site3','Site4','Site5','Site6','n_reads','\n']))
    
    # count = 0
    # for k,v in TAPE_table.items():
    #     if count < int(Read_threshold):
    #         if v > read_cutoff:
    #             output_file.write(k+','+str(v)+'\n')
    #             count += 1


    count = 0
    for k,v in TAPE_table2.items():
        #row = k.split(',')
        if count < 1000: #Arbitrary cutoff to limit the script running forever?
            if v > read_cutoff:
                output_file.write(k+','+str(v)+'\n')
                count += 1
            

    
    print("Final file has been saved.")



