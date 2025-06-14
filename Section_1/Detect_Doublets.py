

### We performed two steps with the goal of exhaustively detecting and removing potential doublets. 


###########################################
### running scrublet to detect doublets ###
###########################################

import scanpy as sc
import pandas as pd
import numpy as np
import scrublet as scr
import os
import sys

batch_id = str(int(sys.argv[1]))

work_path = "/Your_work_path/data_sci/plate_" + batch_id

adata = sc.read_mtx(os.path.join(work_path, 'gene_count.mtx'))
pdata = pd.read_csv(os.path.join(work_path, 'df_cell.csv'), index_col = 0)
fdata = pd.read_csv(os.path.join(work_path, 'df_gene.csv'), index_col = 0)
adata.obs_names = list(pdata['sample'])
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

pd.DataFrame(scrublet_results['doublet_scores_observed_cells']).to_csv(os.path.join(work_path, "doublet_scores_observed_cells.csv"), index = False, header = None)
pd.DataFrame(scrublet_results['doublet_scores_simulated_doublets']).to_csv(os.path.join(work_path, "doublet_scores_simulated_doublets.csv"), index = False, header = None)




###########################################
### running scrublet to detect doublets ###
###########################################

import scanpy as sc
import pandas as pd
import numpy as np
import scrublet as scr
import os, sys

batch_id = str(int(sys.argv[1]))

work_path = "/Your_work_path/data_sci/plate_" + batch_id

os.mkdir(os.path.join(work_path, "doublet_cluster"))

adata = sc.read_mtx(os.path.join(work_path, 'gene_count.mtx'))
pdata = pd.read_csv(os.path.join(work_path, 'df_cell.csv'), index_col = 0)
fdata = pd.read_csv(os.path.join(work_path, 'df_gene.csv'), index_col = 0)
fdata.index = fdata['gene_ID']
adata.obs = pdata
adata.var = fdata

adata_orig = adata

### remove sex genes
adata = adata_orig[:, ~adata_orig.var['chr'].isin(['chrX', 'chrY'])]
### high variable genes
sc.pp.filter_genes(adata, min_cells=1)
sc.pp.normalize_total(adata, target_sum=1e5)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=3000)
filter_genes = list(adata.var.loc[adata.var['highly_variable'] == True, 'gene_ID'])

### 
adata = adata_orig[:, adata_orig.var['gene_ID'].isin(filter_genes)]
sc.pp.normalize_total(adata, target_sum=1e5)
sc.pp.log1p(adata)
sc.pp.scale(adata)
###
sc.tl.pca(adata, svd_solver='arpack', n_comps = 30)
sc.pp.neighbors(adata, n_neighbors=50, n_pcs=30)
sc.tl.louvain(adata)
sc.tl.umap(adata, min_dist=0.1)

adata.obs['umap_1'] = list(adata.obsm['X_umap'][:,0])
adata.obs['umap_2'] = list(adata.obsm['X_umap'][:,1])
name = "global.csv"
adata.obs.to_csv(os.path.join(work_path, 'doublet_cluster', name))

obs_all = adata.obs
obs_all['louvain'].value_counts()
cluster_list = list(set(list(obs_all['louvain'])))

xx = 0
for cnt in range(len(cluster_list)):
    xx += 1
    print('Processing: ' + str(xx) + '/' + str(len(cluster_list)))
    cluster = cluster_list[cnt]
    include_cell = list(obs_all.loc[obs_all['louvain'] == cluster, 'sample'])
    adata = adata_orig[adata_orig.obs['sample'].isin(include_cell)]
    adata = adata[:, ~adata.var['chr'].isin(['chrX', 'chrY'])]
    sc.pp.filter_genes(adata, min_cells=1)
    sc.pp.normalize_total(adata, target_sum=1e5)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=3000)
    filter_genes = list(adata.var.loc[adata.var['highly_variable'] == True, 'gene_ID'])
    adata = adata_orig[adata_orig.obs['sample'].isin(include_cell)]
    adata = adata[:, adata.var['gene_ID'].isin(filter_genes)]
    sc.pp.normalize_total(adata, target_sum=1e5)
    sc.pp.log1p(adata)
    sc.pp.scale(adata)
    sc.tl.pca(adata, svd_solver='arpack', n_comps = 30)
    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30)
    sc.tl.louvain(adata, resolution = 3)
    sc.tl.umap(adata, min_dist=0.1)
    adata.obs['umap_1'] = list(adata.obsm['X_umap'][:,0])
    adata.obs['umap_2'] = list(adata.obsm['X_umap'][:,1])
    name = 'adata.obs.louvain_' + cluster + '.csv'
    adata.obs.to_csv(os.path.join(work_path, 'doublet_cluster', name))

