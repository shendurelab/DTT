
### step-3: here we ignore the debris-seq data, calculating how many TapeBC + first_edit are overlapped between each two cells
### we are going to perform the barcode correction, that allows 1bp mismatch on the TapeBC if their first edit are exactly same.

import os, sys

work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape/data_sci"

rna = {}
file = open(os.path.join(work_path, "tape_barcode_2", "tape_barcode.txt"))
for line in file:
    l = line.rstrip().split('\t')
    rna[l[0]] = rna.get(l[0], set())
    rna[l[0]].add(l[1] + '-' + l[2][:3])
file.close()

cell_list = list(rna.keys())

out_1 = open(os.path.join(work_path, "tape_barcode_2", "cell_cell_list.txt"), "w")
for cell_i in range(0, len(cell_list)):
    out_1.write(cell_list[cell_i] + "\n")
out_1.close()

print(len(rna))

def mismatch(seed_seq):
    bc_seq = set()
    for index in range(len(seed_seq)):
        for rep_index in ['A', 'T', 'C', 'G', 'N']:
            tmp_seq = seed_seq[:index] + rep_index + seed_seq[index + 1:]
            bc_seq.add(tmp_seq)
    return(bc_seq)

out_2 = open(os.path.join(work_path, "tape_barcode_2", "cell_cell_1_mismatch.txt"), 'w')
for cell_i in range(0, (len(cell_list)-1)):
    
    cell_ex = set()
    for i in rna[cell_list[cell_i]]:
        bc = i.split('-')[0]
        edit = '-' + i.split('-')[1]
        tmp = mismatch(bc)
        cell_ex = cell_ex.union(set([x + edit for x in tmp]))

    for cell_j in range((cell_i+1), len(cell_list)):
        if len(cell_ex.intersection(rna[cell_list[cell_j]])) >= 2:   ### we only output >=2, otherwise the result is too big
            out_2.write(str(cell_i+1) + '\t' + str(cell_j+1) + '\t' + str(len(cell_ex.intersection(rna[cell_list[cell_j]]))) + '\n')

out_2.close()


### step-4: performing semi-supervised UMAP on the distance matrix generated from cell_cell_1_mismatch.txt
### of note, it works not well, do give up with this analysis

### first, adding two lines as header:
### %%MatrixMarket matrix coordinate integer general
### 229507 229507 332055259

### > cat tmp cell_cell_1_mismatch.txt > cell_pairs.mtx

from scipy.io import mmread
import sys, os
from scipy.sparse import triu, tril
from scipy.sparse import csr_matrix
from umap import UMAP
import numpy as np
import pandas as pd

work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape/data_sci"

cell_list_label = pd.read_csv(os.path.join(work_path, "tape_barcode_2", "cell_cell_list_label.csv"), index_col = 0)

x = mmread(os.path.join(work_path, "tape_barcode_2", "cell_pairs.mtx"))
x.shape

x_l = x.transpose()
x_f = triu(x) + tril(x_l)
x_f.shape

x_f.data = 1.0/(x_f.data + 1) - 1
x_f.setdiag(-1)

cutoff_list = [10, 30, 50]
cutoff = cutoff_list[int(sys.argv[1]) - 1]
print(cutoff)

nonzero_counts = np.diff(x_f.indptr)
rows_to_keep = np.where(nonzero_counts >= cutoff)[0]
x_sub = x_f[rows_to_keep][:, rows_to_keep]
x_sub.shape

umap = UMAP(n_components = 2,
            n_neighbors = cutoff,
            min_dist = 0.1)
embedding = umap.fit_transform(x_sub, y = np.array(cell_list_label['sample_id'])[rows_to_keep])

np.savetxt(os.path.join(work_path, "tape_barcode_2", "umap_embedding.csv.cutoff_%s"%str(cutoff)), embedding, delimiter=",", fmt='%1.3f')











