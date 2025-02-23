### Here I am going to try a different way, first generating a TapeBC x cell matrix, with values as the read count.
### Then, performing PCA followed by clustering on it. Finally, we map cell clusters (rather than individual cells) to each gastruloids.


### Step-1: we extracted the TapeBC + edit information for individual cells from sci-RNA-seq3 outputs under filtered_sam folder

import sys, os
import re

data_path_list = {}
data_path_list["plate_1"] = "./2023-10-03-sciA-samchoi/nobackup/output-enrich/filtered_sam"
data_path_list["plate_2"] = "./2023-10-03-sciB-samchoi/nobackup/output-enrich/filtered_sam"
data_path_list["plate_3"] = "./2023-10-27-sciA-samchoi/nobackup/output-enrich/filtered_sam"
data_path_list["plate_4"] = "./2023-10-27-sciB-samchoi/nobackup/output-enrich/filtered_sam"
data_path_list["plate_5"] = "./2023-11-01-sciA-samchoi/nobackup/output-enrich/filtered_sam"
data_path_list["plate_6"] = "./2023-11-01-sciB-samchoi/nobackup/output-enrich/filtered_sam"

work_path = "Your_work_path"

plate_id = "plate_" + str(int(sys.argv[1]))
data_path = data_path_list[plate_id]

file_list = os.listdir(data_path)

file_id = file_list[int(sys.argv[2])-1]
name = file_id.replace('.sam', '')
print(file_id)

seed_seq = 'TGATGGTGAGCACG'
poss_seq = set()
for index in range(len(seed_seq)):
    for rep_index in ['A', 'T', 'C', 'G', 'N']:
        tmp_seq = seed_seq[:index] + rep_index + seed_seq[index + 1:]
        poss_seq.add(tmp_seq)

def extract_edit(x):
    y = []
    if x[5:7] == "AA" and x[12:16] == "GGAT" and x[15:29] in poss_seq:
        y.append(x[:5] + x[7:12])
        tmp = x[12:]
        if tmp[0:4] == "GGAT" and tmp[20:24] == "GGAT":
            y.append(tmp[17:20]); tmp = tmp[20:]
            if tmp[0:4] == "GGAT" and tmp[20:24] == "GGAT":
                y.append(tmp[17:20]); tmp = tmp[20:]
                if tmp[0:4] == "GGAT" and tmp[20:24] == "GGAT":
                    y.append(tmp[17:20]); tmp = tmp[20:]
                    if tmp[0:4] == "GGAT" and tmp[20:24] == "GGAT":
                        y.append(tmp[17:20]); tmp = tmp[20:]
                        if tmp[0:4] == "GGAT" and tmp[20:24] == "GGAT":
                            y.append(tmp[17:20]); tmp = tmp[20:]
                            if tmp[0:4] == "GGAT" and tmp[20:24] == "GGAT":
                                y.append(tmp[17:20]); tmp = tmp[20:]
                            else:
                                y.append("None")
                        else:
                            y.append("None"); y.append("None")
                    else:
                        y.append("None"); y.append("None"); y.append("None")
                else:
                    y.append("None"); y.append("None"); y.append("None"); y.append("None")
            else:
                y.append("None"); y.append("None"); y.append("None"); y.append("None"); y.append("None")
        else:
            y.append("None"); y.append("None"); y.append("None"); y.append("None"); y.append("None"); y.append("None")
    return(y)

dat = {}
file = open(os.path.join(data_path, file_id))
line = file.readline()
while line:
    l = line.rstrip().split('\t')
    if len(l) >= 10:
        if l[2] == "SamChoiTAPE":
            tmp = l[0].split(',')[0]
            dat[tmp] = dat.get(tmp, [])
            dat[tmp].append(l[9])          ### save all the possible TapeBC+edit in a list instead of set, since we want to count the read number
    line = file.readline()
file.close()

out = open(os.path.join(work_path, "%s.wo_correction.txt"%name), 'w') ### only output result w/o barcode correction
for i in dat:     ### i refers to a single cell
    tmp = {}
    for j in dat[i]:
        y = extract_edit(j)
        if len(y) > 0:
            x = (y[0], ','.join(y[1:]))    ### tmp[(BC, edit)] = count
            tmp[x] = tmp.get(x, 0)
            tmp[x] += 1
    for k in tmp:
        out.write(name + '.' + i + '\t' + k[0] + '\t' + k[1] + '\t' + str(tmp[k]) + '\n')
out.close()


### Step-2: we created a TapeBC x Cell matrix, in which the values are the read count.

### The TapeBC list was obtained from the Debris-seq. 
### When we mapped the TapeBC from Debris-seq and RNA-seq, we also performed a barcode correction 
### based on the 1-bp mismatch on the TapeBC if their first edit is the same.

import os, sys

work_path = "Your_work_path"

def mismatch(seed_seq, edit_4sites):
    bc_seq = set()
    for index in range(len(seed_seq)):
        for rep_index in ['A', 'T', 'C', 'G', 'N']:
            tmp_seq = seed_seq[:index] + rep_index + seed_seq[index + 1:]
            bc_seq.add(tmp_seq + '-' + edit_4sites)
    return(bc_seq)

### read TapeBC list from Debris-seq output (TapeBC + edit_4sites)
debris_dict = {}
file = open(os.path.join(work_path, "tape_barcode_2", "debris_sample_barcode_uniq.txt"))
for line in file:
    l = line.rstrip().split('\t')
    x = l[1].split('-')
    debris_dict['-'.join([x[0], x[1], x[2], x[3], x[4]])] = l[1]
    debris_dict['-'.join([x[0], x[1], x[2], x[3], 'Non'])] = l[1]
    debris_dict['-'.join([x[0], x[1], x[2], 'Non', 'Non'])] = l[1]
    debris_dict['-'.join([x[0], x[1], 'Non', 'Non', 'Non'])] = l[1]
file.close()

### read white list of cells for individual PCR samples
file = open(os.path.join(work_path, "cell_id.txt"))
PCR_dict = {}
for line in file:
    l = line.rstrip().split('.')
    PCR_dict[l[0]] = PCR_dict.get(l[0], set())
    PCR_dict[l[0]].add(line.rstrip())
file.close()

out = open(os.path.join(work_path, "tape_barcode_3", "Tape_cell_uniq.txt"), 'w')
iter_x = 1
for PCR_id in PCR_dict:
    print(str(iter_x) + ' : ' + PCR_id)
    iter_x += 1

    ### read the TapeBC + edit_4sites from the scRNA-seq
    file = open(os.path.join(work_path, "tape_barcode_1", "%s.wo_correction.txt"%PCR_id))
    cell_debris_map = {}
    for line in file:
        l = line.rstrip().split('\t')
        if l[0] in PCR_dict[PCR_id]:
            edit_4sites = '-'.join([n[0:3] for n in l[2].split(',')][0:4])
            seq_pool = mismatch(l[1], edit_4sites).intersection(debris_dict)
            if len(seq_pool) != 0:
                for cnt in seq_pool:
                    cell_debris_map[(debris_dict[cnt], l[0])] = cell_debris_map.get((debris_dict[cnt], l[0]), 0)
                    cell_debris_map[(debris_dict[cnt], l[0])] += int(l[3])
    file.close()

    ### output
    for cell in cell_debris_map:
        out.write(cell[0] + '\t' + cell[1]  + '\t' + str(cell_debris_map[cell]) + '\n')
out.close()




