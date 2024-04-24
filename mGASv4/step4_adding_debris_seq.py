### step-1: extracting cells which have been profiled with scRNA-seq (after QC filtering)

import os, sys

work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape/data_sci"

file = open(os.path.join(work_path, "cell_id.txt"))
cell_id = {}
for line in file:
    l = line.rstrip().split('.')
    cell_id[l[0]] = cell_id.get(l[0], set())
    cell_id[l[0]].add(line.rstrip())
file.close()

out = open(os.path.join(work_path, "tape_barcode_2", "tape_barcode.txt"), 'w')
cnt = 0
for pcr_id in cell_id:
    cnt += 1
    cell_bc = {}
    file = open(os.path.join(work_path, "tape_barcode_1", "%s.txt"%pcr_id))
    for line in file:
        l = line.rstrip().split('\t')
        if l[0] in cell_id[pcr_id]:
            out.write(line)
out.close()
print(cnt)


out = open(os.path.join(work_path, "tape_barcode_2", "tape_barcode.wo_correction.txt"), 'w')
cnt = 0
for pcr_id in cell_id:
    cnt += 1
    cell_bc = {}
    file = open(os.path.join(work_path, "tape_barcode_1", "%s.wo_correction.txt"%pcr_id))
    for line in file:
        l = line.rstrip().split('\t')
        if l[0] in cell_id[pcr_id]:
            out.write(line)
out.close()
print(cnt)



### step-2: aligning cells from RNA-seq to debris-seq, based on TapeBC + first_edit
### for TapeBC, we retain either perfect match or 1-bp mismatch

import os, sys

work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape/data_sci"

debris = {}
debris_list = []
file = open(os.path.join(work_path, "tape_barcode_2", "debris_sample_barcode.txt"))
for line in file:
    l = line.rstrip().split('\t')
    if len(l[1]) == 16:
        bc = l[1][:5] + l[1][7:16]
        debris[l[0]] = debris.get(l[0], set())
        debris[l[0]].add(bc)
        if l[0] not in debris_list:
            debris_list.append(l[0])
file.close()

print(len(debris_list))

rna = {}
file = open(os.path.join(work_path, "tape_barcode_2", "tape_barcode.txt"))
for line in file:
    l = line.rstrip().split('\t')
    rna[l[0]] = rna.get(l[0], set())
    rna[l[0]].add(l[1] + '-' + l[2][:3])
file.close()

print(len(rna))

out_1 = open(os.path.join(work_path, "tape_barcode_2", "cell_debris_well_perfect.txt"), 'w')
out_1.write('cell' + '\t' + '\t'.join(debris_list) + '\n')
for cell in rna:
    x_per = []
    for well in debris_list:
        x_per.append(str(len(rna[cell].intersection(debris[well]))))
    out_1.write(cell + '\t' + '\t'.join(x_per) + '\n')
out_1.close()

def mismatch(seed_seq):
    bc_seq = set()
    for index in range(len(seed_seq)):
        for rep_index in ['A', 'T', 'C', 'G', 'N']:
            tmp_seq = seed_seq[:index] + rep_index + seed_seq[index + 1:]
            bc_seq.add(tmp_seq)
    return(bc_seq)

out_2 = open(os.path.join(work_path, "tape_barcode_2", "cell_debris_well_1_mismatch.txt"), 'w')
out_2.write('cell' + '\t' + '\t'.join(debris_list) + '\n')
for cell in rna:
    cell_ex = set()
    for i in rna[cell]:
        bc = i.split('-')[0]
        edit = '-' + i.split('-')[1]
        tmp = mismatch(bc)
        cell_ex = cell_ex.union(set([x + edit for x in tmp]))
    x_mis = []
    for well in debris_list:
        x_mis.append(str(len(cell_ex.intersection(debris[well]))))
    out_2.write(cell + '\t' + '\t'.join(x_mis) + '\n')
out_2.close()







