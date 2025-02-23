
##################################################
### Focusing on ancestor nodes with top 1000 depth

import sys, os
import pandas as pd
import numpy as np
import time
from scipy.io import mmread
from scipy.sparse import eye, hstack, coo_matrix, csr_matrix, lil_matrix

work_path = "PATH_1"
save_path = "PATH_2"

dat_orig = mmread(os.path.join(work_path, "tree_sanjay", "ancestor_matrix_sub.sparseMatrix.mtx")).tocsr().T

kk = int(sys.argv[1]) ### 1-100

dat_anc = mmread(os.path.join(save_path, "tree_sanjay", "bootstrap", "ancestor_matrix_%s.mtx" % str(kk))).tocsr()
dat_tip = eye(dat_anc.shape[0], format='csr')
dat = hstack([dat_tip, dat_anc], format='csr')

chunk_size = 50
N_row = dat_orig.shape[0]
N_col = dat.shape[1]

N_row_chunk = len(range(0, N_row, chunk_size)); print(N_row_chunk)
N_col_chunk = len(range(0, N_col, chunk_size)); print(N_col_chunk)

result_matrix = lil_matrix((N_row, N_col_chunk), dtype=int)

start_time = time.time()
for i in range(0, N_row, chunk_size):
    dat_orig_sub = dat_orig[i : min(i + chunk_size, N_row), :]
    dat_orig_rev = csr_matrix(np.ones(dat_orig_sub.shape)) - dat_orig_sub
    cnt = 0
    for j in range(0, N_col, chunk_size):
        print(str(i)+'/'+str(cnt))
        dat_sub = dat[:, j : min(j + chunk_size, N_col)]
        dat_rev = csr_matrix(np.ones(dat_sub.shape)) - dat_sub
        dat_1 = dat_orig_sub.dot(dat_sub) + dat_orig_rev.dot(dat_rev)
        dat_2 = dat_orig_rev.dot(dat_sub) + dat_orig_sub.dot(dat_rev)
        dat_min = hstack([dat_1, dat_2], format='csr').min(axis=1).toarray().flatten().astype(int)
        result_matrix[i : min(i + chunk_size, N_row), cnt] = dat_min[:, None]
        cnt += 1
        end_time = time.time()
        print(f"Running time: {((end_time - start_time)/60):.3f} minutes")

row_min_values = result_matrix.toarray().min(axis=1)
df = pd.DataFrame({
    'node': np.arange(1, len(row_min_values) + 1),
    'score': row_min_values
})
output_csv = os.path.join(save_path, "tree_sanjay", "bootstrap", "TBE_%s.csv" % str(kk))
df.to_csv(output_csv, index=False)


#########################################################################
### Focusing on all the ancestor nodes, to calculate TBE (it's very slow)

import sys, os
import pandas as pd
import numpy as np
from scipy.io import mmread
from scipy.sparse import eye, hstack, coo_matrix, csr_matrix, lil_matrix

work_path = "PATH_1"
save_path = "PATH_2"

dat_orig = mmread(os.path.join(work_path, "tree_sanjay", "ancestor_matrix.mtx")).tocsr().T
#dat_orig = mmread(os.path.join("/net/gs/vol1/home/cxqiu/share", "ancestor_matrix.mtx")).tocsr().T

kk = int(sys.argv[1]) ### 1-100

dat_anc = mmread(os.path.join(save_path, "tree_sanjay", "bootstrap", "ancestor_matrix_%s.mtx" % str(kk))).tocsr()
#dat_anc = mmread(os.path.join("/net/gs/vol1/home/cxqiu/share", "ancestor_matrix_1.mtx")).tocsr()
dat_tip = eye(dat_anc.shape[0], format='csr')
dat = hstack([dat_tip, dat_anc], format='csr')

chunk_size = 1000
N_row = dat_orig.shape[0]
N_col = dat.shape[1]

N_row_chunk = len(range(0, N_row, chunk_size)); print(N_row_chunk)
N_col_chunk = len(range(0, N_col, chunk_size)); print(N_col_chunk)

result_matrix = lil_matrix((N_row, N_col_chunk), dtype=int)

for i in range(0, N_row, chunk_size):
    print(str(i)+'/'+str(N_row_chunk))
    dat_orig_sub = dat_orig[i : min(i + chunk_size, N_row), :]
    dat_orig_rev = csr_matrix(np.ones(dat_orig_sub.shape)) - dat_orig_sub
    cnt = 0
    for j in range(0, N_col, chunk_size):
        print(str(i)+'/'+str(j))
        dat_sub = dat[:, j : min(j + chunk_size, N_col)]
        dat_rev = csr_matrix(np.ones(dat_sub.shape)) - dat_sub
        dat_1 = dat_orig_sub.dot(dat_sub) + dat_orig_rev.dot(dat_rev)
        dat_2 = dat_orig_rev.dot(dat_sub) + dat_orig_sub.dot(dat_rev)
        dat_min = hstack([dat_1, dat_2], format='csr').min(axis=1).toarray().flatten().astype(int)
        result_matrix[i : min(i + chunk_size, N_row), cnt] = dat_min[:, None]
        cnt += 1

row_min_values = result_matrix.toarray().min(axis=1)
df = pd.DataFrame({
    'node': np.arange(1, len(row_min_values) + 1),
    'score': row_min_values
})
output_csv = os.path.join(save_path, "tree_sanjay", "bootstrap", "TBE_%s.csv" % str(kk))
df.to_csv(output_csv, index=False)

