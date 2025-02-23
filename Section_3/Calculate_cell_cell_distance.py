
#############################################################
### STEP-1:
### For each single DTT, we calculated: 
### 1) cell-cell distance; 2) cell-cell if they are both None

import numpy as np
import csv
import os, sys
import time

work_path = "PATH_1"
save_path = "PATH_2"

DTT_i = int(sys.argv[1]) ### 1-231
step = 6

c = 0
cell_edit_table = {}
cell_list = []
with open(os.path.join(work_path, "tree_sanjay", "loci.csv")) as f:
    reader = csv.reader(f)
    for row in reader:
        if c == 0:
            c += 1
            continue
        else:
            cell_edit_table[row[0]] = row[((DTT_i - 1) * step + 1):(DTT_i * step + 1)]
            cell_list.append(row[0])

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

start_time = time.time()

orig_dis_matrix = np.zeros((len(cell_list), len(cell_list)), dtype=np.int8)
comp_num_matrix = np.zeros((len(cell_list), len(cell_list)), dtype=np.int8)

for ii in range(len(cell_list)):
    print(ii)
    a = cell_edit_table[cell_list[ii]]
    for jj in range(ii+1, len(cell_list)):
        b = cell_edit_table[cell_list[jj]]
        dist_ij = calculate_distance(a, b)
        if a[0] != 'None' and b[0] != 'None':
            comp_ij = 1
        else:
            comp_ij = 0
        orig_dis_matrix[ii,jj] = orig_dis_matrix[jj,ii] = dist_ij
        comp_num_matrix[ii,jj] = comp_num_matrix[jj,ii] = comp_ij

end_time = time.time()
print(end_time - start_time)

np.savez_compressed(os.path.join(save_path, 'tree_sanjay', 'each_DTT', 'cell_cell_distance_%s.npz'%str(DTT_i)), arr1=orig_dis_matrix, arr2=comp_num_matrix)




########################################################
### STEP-2:
### Aggregating 231 DTTs together, as cell-cell distance

import numpy as np
import os, sys

work_path = "PATH_1"
save_path = "PATH_2"

cell_num = 58283
orig_dis_matrix = np.zeros((cell_num, cell_num), dtype=np.int16)
comp_num_matrix = np.zeros((cell_num, cell_num), dtype=np.int16)

N = 231

for DTT_i in range(N):
    print(DTT_i + 1)
    dat = np.load(os.path.join(save_path, 'tree_sanjay', 'each_DTT', 'cell_cell_distance_%s.npz'%str(DTT_i + 1)))
    orig_dis_matrix += dat['arr1']
    comp_num_matrix += dat['arr2']

np.savetxt(os.path.join(work_path, 'tree_sanjay', 'cell_cell_dis.csv'), orig_dis_matrix, delimiter=',', fmt='%d')
np.savetxt(os.path.join(work_path, 'tree_sanjay', 'cell_cell_num.csv'), comp_num_matrix, delimiter=',', fmt='%d')




###########################
### STEP-3:
### Bootstrapping 100 times

import numpy as np
import os, sys
import random

work_path = "PATH_1"
save_path = "PATH_2"

boot_i = int(sys.argv[1]) ### 1-100

cell_num = 58283
orig_dis_matrix = np.zeros((cell_num, cell_num), dtype=np.int16)
comp_num_matrix = np.zeros((cell_num, cell_num), dtype=np.int16)

N = 231

bootstrapped_indices = random.choices(range(N), k=N)

for DTT_i in bootstrapped_indices:
    print(DTT_i + 1)
    dat = np.load(os.path.join(save_path, 'tree_sanjay', 'each_DTT', 'cell_cell_distance_%s.npz'%str(DTT_i + 1)))
    orig_dis_matrix += dat['arr1']
    comp_num_matrix += dat['arr2']

np.savetxt(os.path.join(save_path, 'tree_sanjay', 'bootstrap', 'cell_cell_dis_%s.csv'%str(boot_i)), orig_dis_matrix, delimiter=',', fmt='%d')
np.savetxt(os.path.join(save_path, 'tree_sanjay', 'bootstrap', 'cell_cell_num_%s.csv'%str(boot_i)), comp_num_matrix, delimiter=',', fmt='%d')














