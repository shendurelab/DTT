### Here I am going to create a TapeBC x Cell matrix, with UMIs as values.

### The TapeBC list was obtained from the Debris-seq. 
### When we mapped the TapeBC from Debris-seq and RNA-seq, we also performed a barcode correction 
### based on the 1-bp mismatch on the TapeBC if their first edit is the same.

import os, sys

def mismatch(seed_seq, edit_4sites):
    bc_seq = set()
    for index in range(len(seed_seq)):
        for rep_index in ['A', 'T', 'C', 'G', 'N']:
            tmp_seq = seed_seq[:index] + rep_index + seed_seq[index + 1:]
            bc_seq.add(tmp_seq + '-' + edit_4sites)
    return(bc_seq)


### read TapeBC list from Debris-seq output (TapeBC + edit_4sites)
debris_dict = {}
file = open(os.path.join(work_path, "cell_to_gastruloid", "debris_sample_barcode.txt"))
for line in file:
    l = line.rstrip().split('\t')
    x = l[1].split('-')
    debris_dict['-'.join([x[0], x[1], x[2], x[3], x[4]])] = l[1]
    debris_dict['-'.join([x[0], x[1], x[2], x[3], 'Non'])] = l[1]
    debris_dict['-'.join([x[0], x[1], x[2], 'Non', 'Non'])] = l[1]
    debris_dict['-'.join([x[0], x[1], 'Non', 'Non', 'Non'])] = l[1]
file.close()

### read white list of cells for individual lanes, after filtering based on their transcriptome
file = open(os.path.join(work_path, "cell_id.txt"))
PCR_dict = {}
for line in file:
    l = line.rstrip().split('\t')
    x = l[1].split('_')[1]
    PCR_dict[x] = PCR_dict.get(x, set())
    PCR_dict[x].add(l[0])
file.close()

out = open(os.path.join(work_path, "cell_to_gastruloid", "Tape_cell.txt"), 'w')
iter_x = 1
for PCR_id in PCR_dict:
    print(str(iter_x) + ' : ' + PCR_id)

    if iter_x == 1:
    ### read the TapeBC + edit_4sites from the scRNA-seq
        file = open(os.path.join(work_path, "TapeBC_data", "mGASv3_Lane1_CellByTape_10X_bamExtractV3_t100_Site1collapsed2_v240330.csv"))
    else:
        file = open(os.path.join(work_path, "TapeBC_data", "mGASv3_Lane2_CellByTape_10X_bamExtractV3_t100_Site1collapsed2_v240331.csv"))

    line = file.readline()
    line = file.readline()
    cell_debris_map = {}
    while line:
        l = line.rstrip().split(',')
        if l[0] in PCR_dict[PCR_id]:
            edit_4sites = '-'.join([n[0:3] for n in l[2:6]])
            x = l[1][0:5] + l[1][7:12]
            seq_pool = mismatch(x, edit_4sites).intersection(debris_dict)
            if len(seq_pool) != 0:
                if len(seq_pool) > 1:
                    print(seq_pool)

                for cnt in seq_pool:
                    cell_debris_map[(debris_dict[cnt], l[0])] = cell_debris_map.get((debris_dict[cnt], l[0]), 0)
                    cell_debris_map[(debris_dict[cnt], l[0])] += int(l[9])
        line = file.readline()
    file.close()

    ### output
    for cell in cell_debris_map:
        out.write(cell[0] + '\t' + cell[1]  + '\t' + str(cell_debris_map[cell]) + '\n')

    iter_x += 1
out.close()




