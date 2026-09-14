
##########################################################################
### Assigning cells from scRNA-seq of monoclonal gastruloids to debris-seq

work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape"
data_path = "/net/shendure/vol10/projects/prime_editing/nobackup/data/2025_11_21_revision"
source("~/work/scripts/utils.R")
library(Signac)

obj = readRDS(paste0(work_path, "/revision/data_10X_multiome/obj_combine_processed.rds"))
emb = Embeddings(obj, reduction = "umap")
pd = data.frame(obj[[]])
pd$UMAP_1 = emb[,1]
pd$UMAP_2 = emb[,2]

well_barcode = read.csv(paste0(work_path, "/revision/debris_seq/well_barcode_set.csv"), header=T)
well_barcode$well_cluster = paste0("cluster_", well_barcode$well_cluster)

well_barcode_content = unique(well_barcode[,c("well_id", "well_cluster")]) 
well_barcode_content_num = well_barcode_content %>% group_by(well_cluster) %>% tally()
well_cluster_list = paste0("cluster_", 1:nrow(well_barcode_content_num))

run_list = paste0("Step1_MT0", 1:4, "_CellByTapeV3_t2.csv")

dat_all = NULL
for(kk in 1:length(run_list)){
    print(paste0(kk, "/", length(run_list)))
    dat_i = read.csv(paste0(data_path, "/data/", run_list[kk]))[,c(1:4)]
    dat_i$cell_id = unlist(lapply(as.vector(dat_i$cBC), function(x) strsplit(x,"[-]")[[1]][1]))
    dat_i$cell_id = paste0(dat_i$cell_id, "-", kk)
    dat_i$barcode = dat_i$mBC
    dat_i = dat_i[dat_i$cell_id %in% pd$cell_id & dat_i$barcode %in% well_barcode$barcode,]
    print(nrow(dat_i))
    dat_all = rbind(dat_all, dat_i[,c("cell_id", "barcode", "n_UMI")])
}

well_barcode = well_barcode[well_barcode$barcode %in% dat_all$barcode,]


##########################################
### Method-2: also consider the read count

dat_1 = well_barcode[,c("well_id", "barcode", "read_count")] %>% dcast(barcode ~ well_id, fill = 0)
rownames(dat_1) = dat_1[,1]; dat_1 = as.matrix(dat_1[,-1])
dat_1 = t(t(dat_1)/apply(dat_1, 2, sum))

dat_2 = dat_all[,c("cell_id", "barcode", "n_UMI")] %>% dcast(barcode ~ cell_id, fill = 0)
rownames(dat_2) = dat_2[,1]; dat_2 = as.matrix(dat_2[,-1])
dat_2 = t(t(dat_2)/apply(dat_2, 2, sum))

### 934 barcodes x 100 wells x 55694 cells
mat = t(dat_2) %*% dat_1

mat_regroup = NULL
for(i in well_cluster_list){
  print(i)
  well_include = well_barcode_content %>% filter(well_cluster == i) %>% pull(well_id)
  mat_regroup = cbind(mat_regroup, apply(mat[,colnames(mat) %in% well_include, drop = F], 1, mean))
}

rownames(mat_regroup) = rownames(mat)
colnames(mat_regroup) = well_cluster_list
mat = mat_regroup

max_idx = max.col(mat, ties.method = "first")
max_val = mat[cbind(seq_len(nrow(mat)), max_idx)]
mat2 = mat
mat2[cbind(seq_len(nrow(mat)), max_idx)] = -Inf
second_idx = max.col(mat2, ties.method = "first")
second_val = mat2[cbind(seq_len(nrow(mat)), second_idx)]

df = data.frame(cell_id = rownames(mat),
    max_idx = max_idx, max_val = max_val, second_idx = second_idx, second_val = second_val)
df$ratio = df$max_val/df$second_val
df$well_cluster = paste0("cluster_", df$max_idx)

df_sub = df[df$ratio >=3,]
saveRDS(df_sub, paste0(work_path, "/revision/data_10X_multiome/cell_well_cluster_assignment.rds"))
### n = 26,876 cells assigned; 42% of cells were assigned (26876/63889 = 42%)

df_sub_x = readRDS(paste0(work_path, "/revision/monoclonal_RNAseq/cell_well_cluster_assignment.rds"))
well_cluster_list_final = df_sub_x %>% group_by(well_cluster) %>% tally() %>% 
 filter(n >= 200) %>% arrange(-n) %>% pull(well_cluster)
### n = 26 well clusters

df_sub_200 = df_sub[df_sub$well_cluster %in% well_cluster_list_final,]
### 25,366 cells assigned to 26 well clusters (25366/63889 = 40%)
x = df_sub_200 %>% group_by(well_cluster) %>% tally() 
### 99 - 5526, mean 976, sd 1048


### making a heatmap
mat_aggr = NULL
for(i in well_cluster_list_final){
    cell_include = df_sub %>% filter(well_cluster == i) %>% pull(cell_id)
    mat_aggr = rbind(mat_aggr, apply(mat[rownames(mat) %in% cell_include,], 2, mean))
}
rownames(mat_aggr) = well_cluster_list_final
mat_aggr = mat_aggr[,well_cluster_list_final]

library(gplots)
library(RColorBrewer)
library(viridis)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
pdf("~/share/assignment_heatmap.pdf", 8, 5)
heatmap.2(as.matrix(mat_aggr), 
          col=Colors, 
          #col=viridis(100), 
          scale="row", 
          Rowv = F, 
          Colv = F, 
          key=T, 
          density.info="none", 
          trace="none", 
          cexRow = 0.5, 
          cexCol = 0.5,
          margins = c(10,5))
dev.off()



############################
### color the ATAC UMAP using assigned clones

work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape"
data_path = "/net/shendure/vol10/projects/prime_editing/nobackup/data/2025_11_21_revision"
source("~/work/scripts/utils.R")
library(Signac)

obj = readRDS(paste0(work_path, "/revision/data_10X_multiome/obj_combine_processed.rds"))
obj = RunUMAP(object = obj, reduction = 'lsi', dims = 2:50)

pd = data.frame(obj[[]])
pd$UMAP_1_ATAC = Embeddings(obj, reduction = "umap")[,1]
pd$UMAP_2_ATAC = Embeddings(obj, reduction = "umap")[,2]

cell_assignment = readRDS(paste0(work_path, "/revision/data_10X_multiome/cell_well_cluster_assignment.rds"))

pd = pd %>% select(UMAP_1_ATAC, UMAP_2_ATAC, cell_id) %>% 
    left_join(cell_assignment[,c("cell_id", "well_cluster")], by = "cell_id")

saveRDS(pd, paste0(work_path, "/revision/data_10X_multiome/pd_combine_processed.ATAC_umap.rds"))

mat_epiblast = paste0("cluster_", c(8,16,18,20))
mat_hindbrain = paste0("cluster_", c(3,4,5,21,22))
mat_somite = paste0("cluster_", c(9,14,15,19))
mat_other = paste0("cluster_", c(10,11,17,2,23,24,25,26,27,28,29,6,7))

pd$clone_group = rep(NA, nrow(pd))
pd$clone_group[!is.na(pd$well_cluster) & pd$well_cluster %in% mat_epiblast] = "Epiblast"
pd$clone_group[!is.na(pd$well_cluster) & pd$well_cluster %in% mat_hindbrain] = "Hindbrain"
pd$clone_group[!is.na(pd$well_cluster) & pd$well_cluster %in% mat_somite] = "Somites"
pd$clone_group[!is.na(pd$well_cluster) & pd$well_cluster %in% mat_other] = "Others"
pd_sub = pd[!is.na(pd$clone_group),]

clone_color_map = c("#52c05a", "#528ac0", "#dc87ba", "grey50")
names(clone_color_map) = c("Epiblast", "Hindbrain", "Somites", "Others")

for(i in names(clone_color_map)){
    p = ggplot() +
    geom_point(data = pd, aes(x = UMAP_1_ATAC, y = UMAP_2_ATAC), color = "grey80",size = 0.1) +
    geom_point(data = pd_sub[pd_sub$clone_group == i,], aes(x = UMAP_1_ATAC, y = UMAP_2_ATAC, color = clone_group), size = 0.15) +
    theme_void() +
    theme(legend.position="none") +
    scale_color_manual(values = clone_color_map)
    ggsave(paste0("~/share/data_mESCs_ATAC_", i, ".png"), p, width = 3, height = 3, dpi = 300)
}





############################
### color the RNA UMAP using assigned clones

work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape"
data_path = "/net/shendure/vol10/projects/prime_editing/nobackup/data/2025_11_21_revision"
source("~/work/scripts/utils.R")
library(Signac)

obj = readRDS(paste0(work_path, "/revision/data_10X_multiome/obj_combine_processed.rds"))

pd = data.frame(obj[[]])
pd$UMAP_1_RNA = Embeddings(obj, reduction = "umap")[,1]
pd$UMAP_2_RNA = Embeddings(obj, reduction = "umap")[,2]

cell_assignment = readRDS(paste0(work_path, "/revision/data_10X_multiome/cell_well_cluster_assignment.rds"))

pd = pd %>% select(UMAP_1_RNA, UMAP_2_RNA, cell_id) %>% 
    left_join(cell_assignment[,c("cell_id", "well_cluster")], by = "cell_id")

saveRDS(pd, paste0(work_path, "/revision/data_10X_multiome/pd_combine_processed.RNA_umap.rds"))

mat_epiblast = paste0("cluster_", c(8,16,18,20))
mat_hindbrain = paste0("cluster_", c(3,4,5,21,22))
mat_somite = paste0("cluster_", c(9,14,15,19))
mat_other = paste0("cluster_", c(10,11,17,2,23,24,25,26,27,28,29,6,7))

pd$clone_group = rep(NA, nrow(pd))
pd$clone_group[!is.na(pd$well_cluster) & pd$well_cluster %in% mat_epiblast] = "Epiblast"
pd$clone_group[!is.na(pd$well_cluster) & pd$well_cluster %in% mat_hindbrain] = "Hindbrain"
pd$clone_group[!is.na(pd$well_cluster) & pd$well_cluster %in% mat_somite] = "Somites"
pd$clone_group[!is.na(pd$well_cluster) & pd$well_cluster %in% mat_other] = "Others"
pd_sub = pd[!is.na(pd$clone_group),]

clone_color_map = c("#52c05a", "#528ac0", "#dc87ba", "grey50")
names(clone_color_map) = c("Epiblast", "Hindbrain", "Somites", "Others")

for(i in names(clone_color_map)){
    p = ggplot() +
    geom_point(data = pd, aes(x = UMAP_1_RNA, y = UMAP_2_RNA), color = "grey80",size = 0.1) +
    geom_point(data = pd_sub[pd_sub$clone_group == i,], aes(x = UMAP_1_RNA, y = UMAP_2_RNA, color = clone_group), size = 0.15) +
    theme_void() +
    theme(legend.position="none") +
    scale_color_manual(values = clone_color_map)
    ggsave(paste0("~/share/data_mESCs_RNA_", i, ".png"), p, width = 3, height = 3, dpi = 300)
}



