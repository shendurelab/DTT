
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
import random

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

def bootstrap_cell_edit_table(cell_edit_table):
    # Select a reference cell (arbitrary first cell in the dictionary)
    first_cell = next(iter(cell_edit_table))
    values = cell_edit_table[first_cell]
    
    # Determine N (the number of groups of 6)
    N = len(values) // 6
    
    # Perform bootstrap resampling on the reference cell to get N group indices
    bootstrapped_indices = random.choices(range(N), k=N)

    # Initialize a new dictionary for the bootstrapped data
    new_cell_edit_table = {}
    
    for cell, values in cell_edit_table.items():
        # Divide each cell's values into groups of 6 elements
        grouped_values = [values[i * 6:(i + 1) * 6] for i in range(N)]
        
        # Use the same bootstrapped indices for all cells
        bootstrapped_groups = [grouped_values[i] for i in bootstrapped_indices]
        
        # Flatten the list of bootstrapped groups back into a single list
        new_values = [element for group in bootstrapped_groups for element in group]
        
        # Assign the bootstrapped list to the new dictionary
        new_cell_edit_table[cell] = new_values
    
    return new_cell_edit_table

work_path = "./"

boot_i = str(int(sys.argv[1]))

start_time = time.time()

file = open(os.path.join(work_path, "cell_to_gastruloid", "cell_assign.txt"))
cell_include = [line.rstrip().split('\t')[0] for line in file]
file.close()

clone_list = ["clone05", "clone32", "clone25"]

for clone_i in clone_list:

    print(clone_i)

    c = 0
    orig_cell_edit_table = {}
    cell_list = []
    with open(os.path.join(work_path, "tree_sanjay", "loci_%s.csv"%clone_i)) as f:
        reader = csv.reader(f)
        for row in reader:
            if c == 0:
                c += 1
                continue
            else:
                if row[0] in cell_include:
                    orig_cell_edit_table[row[0]] = row[1:]
                    cell_list.append(row[0])

    cell_edit_table = bootstrap_cell_edit_table(orig_cell_edit_table)

    self_counts = np.zeros((len(cell_list),1))
    #distance_matrix = np.zeros((len(cell_list), len(cell_list)))
    comp_num_matrix = np.zeros((len(cell_list), len(cell_list)))
    orig_dis_matrix = np.zeros((len(cell_list), len(cell_list)))

    first_cell = next(iter(cell_edit_table))
    values = cell_edit_table[first_cell]
    TargetBC_size = len(values)

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

            #distance_matrix[ii,jj] = dist_ij / math.sqrt(self_counts[ii] * self_counts[jj] + 1)
            comp_num_matrix[ii,jj] = comp_num_matrix[jj,ii] = comp_ij
            orig_dis_matrix[ii,jj] = orig_dis_matrix[jj,ii] = dist_ij

    print(ii)
    print('DM calculated')

    np.savetxt(os.path.join(work_path, 'tree_sanjay', 'bootstrapping', 'cell_cell_distance', 'cell_cell_orig_distance_%s_%s.csv'%(clone_i, boot_i)), orig_dis_matrix, delimiter=',', fmt='%d')
    np.savetxt(os.path.join(work_path, 'tree_sanjay', 'bootstrapping', 'cell_cell_distance', 'cell_cell_comp_%s_%s.csv'%(clone_i, boot_i)), comp_num_matrix, delimiter=',', fmt='%d')

    

