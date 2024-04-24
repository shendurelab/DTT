
#############################################################
### In this script, we extracted the TapeBC + edit information for individual cells 
### from sci-RNA-seq3 outputs under filtered_sam folder, followed by performing TapeBC correction within the same cell

import sys, os
import re

data_path_list = {}
data_path_list["plate_1"] = "/net/shendure/vol10/projects/gastruloid-TAPE-sci/2023-10-03-sciA-samchoi/nobackup/output-enrich/filtered_sam"
data_path_list["plate_2"] = "/net/shendure/vol10/projects/gastruloid-TAPE-sci/2023-10-03-sciB-samchoi/nobackup/output-enrich/filtered_sam"
data_path_list["plate_3"] = "/net/shendure/vol10/projects/gastruloid-TAPE-sci/2023-10-27-sciA-samchoi/nobackup/output-enrich/filtered_sam"
data_path_list["plate_4"] = "/net/shendure/vol10/projects/gastruloid-TAPE-sci/2023-10-27-sciB-samchoi/nobackup/output-enrich/filtered_sam"
data_path_list["plate_5"] = "/net/shendure/vol10/projects/gastruloid-TAPE-sci/2023-11-01-sciA-samchoi/nobackup/output-enrich/filtered_sam"
data_path_list["plate_6"] = "/net/shendure/vol10/projects/gastruloid-TAPE-sci/2023-11-01-sciB-samchoi/nobackup/output-enrich/filtered_sam"

work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape/data_sci/tape_barcode_1"

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


### We sought to performing TapeBC correction based on the following rule: 
### if two items (TapeBC + edit) within the same cell, sharing the exact same first edit and only 1-bp mismatch at their TapeBC, 
### then we only retain the item with higher read count.

def mismatch(TapeBC_seq, first_edit):
    bc_seq = set()
    for index in range(len(TapeBC_seq)):
        for rep_index in ['A', 'T', 'C', 'G', 'N']:
            tmp_seq = TapeBC_seq[:index] + rep_index + TapeBC_seq[index + 1:]
            bc_seq.add((tmp_seq, first_edit))
    return(bc_seq)


out_1 = open(os.path.join(work_path, "%s.txt"%name), 'w')
out_2 = open(os.path.join(work_path, "%s.wo_correction.txt"%name), 'w') ### also output result w/o barcode correction
for i in dat:     ### i refers to a single cell
    tmp = {}
    for j in dat[i]:
        y = extract_edit(j)
        if len(y) > 0:
            x = (y[0], ','.join(y[1:]))    ### tmp[(BC, edit)] = count
            tmp[x] = tmp.get(x, 0)
            tmp[x] += 1

    for k in tmp:
        out_2.write(name + '.' + i + '\t' + k[0] + '\t' + k[1] + '\t' + str(tmp[k]) + '\n')

    dat_r = {}
    for k in tmp:
        dat_r[(k[0], k[1].split(',')[0])] = dat_r.get((k[0], k[1].split(',')[0]), [])
        dat_r[(k[0], k[1].split(',')[0])].append((k[0], k[1], tmp[k]))   ### dat_r[(BC, first_edit)] = [(BC, edit, count),(BC, edit, count),...]

    include_item = set()
    tmp_uniq = set()

    for m in dat_r: ### m is (BC, first_edit)
        if m not in include_item:
            m_pool = mismatch(m[0], m[1]).intersection(dat_r)
            max_num = 0
            max_item = ""
            for n in m_pool:
                include_item.add(n)
                for x in dat_r[n]:
                    if x[2] > max_num:
                        max_num = x[2]
                        max_item = x

            tmp_uniq.add(max_item)

    for k in tmp_uniq:
        out_1.write(name + '.' + i + '\t' + k[0] + '\t' + k[1] + '\t' + str(k[2]) + '\n')
out_1.close()
out_2.close()




