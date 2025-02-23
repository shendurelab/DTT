
# Importing in tools

import gzip
import subprocess
import os,sys,csv,re
import itertools
from collections import Counter
from collections import OrderedDict
from operator import itemgetter

# Usefule modules
import numpy as np
import edlib
import math
import time

def calculate_distance(a, b):
    cnt = 0
    if a[0] != 'None' and b[0] != 'None':
        if a[0] != 'ETY' or b[0] != 'ETY':
            for i in range(6):
                if (a[i] != 'None' or b[i] != 'None') and a[i] != b[i]:
                    for j in range(i,6):
                        if a[j] != 'None' and a[j] != 'ETY':
                            cnt += 1
                        if b[j] != 'None' and b[j] != 'ETY':
                            cnt += 1
                    break
    return(cnt)

work_path = "./"

start_time = time.time()

file = open(os.path.join(work_path, "cell_to_gastruloid", "cell_assign.txt"))
cell_include = [line.rstrip().split('\t')[0] for line in file]
file.close()

clone_list = ["clone32", "clone05", "clone25"]

clone_i = clone_list[int(sys.argv[1])-1]
print(clone_i)

c = 0
cell_edit_table = {}
cell_list = []
with open(os.path.join(work_path, "tree_sanjay", "loci_%s.csv"%clone_i)) as f:
    reader = csv.reader(f)
    for row in reader:
        if c == 0:
            c += 1
            continue
        else:
            if row[0] in cell_include:
                cell_edit_table[row[0]] = row[1:]
                cell_list.append(row[0])

self_counts = np.zeros((len(cell_list),1))
distance_matrix = np.zeros((len(cell_list), len(cell_list)))
comp_num_matrix = np.zeros((len(cell_list), len(cell_list)))
orig_dis_matrix = np.zeros((len(cell_list), len(cell_list)))
TargetBC_size = len(row[1:])

for ii in range(len(cell_list)):
    match_counts = 0
    for kk in range(0,TargetBC_size,6):
        if (cell_edit_table[cell_list[ii]][kk] != 'None') & (cell_edit_table[cell_list[ii]][kk] != 'ETY'):
            match_counts += 1
            if (cell_edit_table[cell_list[ii]][kk+1] != 'None'):
                match_counts += 1
                if (cell_edit_table[cell_list[ii]][kk+2] != 'None'):
                    match_counts += 1
                    if (cell_edit_table[cell_list[ii]][kk+3] != 'None'):
                        match_counts += 1
                        if (cell_edit_table[cell_list[ii]][kk+4] != 'None'):
                            match_counts += 1
                            if (cell_edit_table[cell_list[ii]][kk+5] != 'None'):
                                match_counts += 1
    self_counts[ii,0] = match_counts ### For individual cell, how many edits

end_time = time.time()
elapsed_time = end_time - start_time

print("Self_counts calculated")
print("Elapsed time so far: " + str(elapsed_time))

for ii in range(len(cell_list)):        
    for jj in range(ii+1,len(cell_list)):
        dist_ij = 0
        comp_ij = 0
        for kk in range(0,TargetBC_size,6):
            a = cell_edit_table[cell_list[ii]][kk:(kk+6)]
            b = cell_edit_table[cell_list[jj]][kk:(kk+6)]
            dist_ij += calculate_distance(a, b)
            
            if a[0] != 'None' and b[0] != 'None':
                comp_ij += 1

        distance_matrix[ii,jj] = distance_matrix[jj,ii] = dist_ij / math.sqrt(self_counts[ii] * self_counts[jj] + 1)
        comp_num_matrix[ii,jj] = comp_num_matrix[jj,ii] = comp_ij
        orig_dis_matrix[ii,jj] = orig_dis_matrix[jj,ii] = dist_ij

print('DM calculated')

np.savetxt(os.path.join(work_path, 'tree_sanjay','cell_cell_comp_%s.csv'%clone_i), comp_num_matrix, delimiter=',', fmt='%d')

np.savetxt(os.path.join(work_path, 'tree_sanjay','cell_cell_orig_distance_%s.csv'%clone_i), orig_dis_matrix, delimiter=',', fmt='%d')

out = open(os.path.join(work_path, 'tree_sanjay','cell_cell_%s.cell_list.csv'%clone_i), 'w')
for i in cell_list:
    out.write(i + '\n')
out.close()

