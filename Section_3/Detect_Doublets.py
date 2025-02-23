
###########################################
### running scrublet to detect doublets ###
###########################################

import scanpy as sc
import pandas as pd
import numpy as np
import scrublet as scr
import os, sys

experiment_id_1 = "data_mGASv5"

experiment_id_2 = "mGASv5_8"

WORK_PATH = './'
adata = sc.read_mtx(os.path.join(WORK_PATH, experiment_id_1, 'processing', experiment_id_2, 'gene_count.mtx'))
pdata = pd.read_csv(os.path.join(WORK_PATH, experiment_id_1, 'processing', experiment_id_2, 'df_cell.csv'), index_col = 0)
fdata = pd.read_csv(os.path.join(WORK_PATH, experiment_id_1, 'processing', experiment_id_2, 'df_gene.csv'), index_col = 0)
adata.obs_names = list(pdata['cell_id'])
adata.var_names = list(fdata['gene_ID'])

min_counts = 3
min_cells = 3
vscore_percentile = 85
n_pc = 30
expected_doublet_rate = 0.06 
sim_doublet_ratio = 2
n_neighbors = 30
scaling_method = 'log'
scrublet_results = scr.compute_doublet_scores(
    adata.X, 
    min_counts = min_counts, 
    min_cells = min_cells, 
    vscore_percentile = vscore_percentile, 
    n_prin_comps = n_pc,
    scaling_method = scaling_method,
    expected_doublet_rate = expected_doublet_rate,
    sim_doublet_ratio = sim_doublet_ratio,
    n_neighbors = n_neighbors, 
    use_approx_neighbors = True, 
    get_doublet_neighbor_parents = False
)

pd.DataFrame(scrublet_results['doublet_scores_observed_cells']).to_csv(os.path.join(WORK_PATH, experiment_id_1, 'processing', experiment_id_2, "doublet_scores_observed_cells.csv"), index = False, header = None)
pd.DataFrame(scrublet_results['doublet_scores_simulated_doublets']).to_csv(os.path.join(WORK_PATH, experiment_id_1, 'processing', experiment_id_2, "doublet_scores_simulated_doublets.csv"), index = False, header = None)

