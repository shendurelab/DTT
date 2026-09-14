
##########################################################################
### Assigning cells from scRNA-seq of monoclonal gastruloids to debris-seq

work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape"
data_path = "/net/shendure/vol10/projects/prime_editing/nobackup/data/2025_11_21_revision"
source("~/work/scripts/utils.R")

pd = readRDS(paste0(work_path, "/revision/monoclonal_polyclonal_RNAseq/obj_batchelor_pd.rds"))

well_barcode = read.csv(paste0(work_path, "/revision/debris_seq/well_barcode_set.csv"), header=T)
well_barcode$well_cluster = paste0("cluster_", well_barcode$well_cluster)

well_barcode_content = unique(well_barcode[,c("well_id", "well_cluster")]) 
well_barcode_content_num = well_barcode_content %>% group_by(well_cluster) %>% tally()
well_cluster_list = paste0("cluster_", 1:nrow(well_barcode_content_num))

run_list = paste0("Mono_L0", 1:9, "_GRCm39")
run_list = c(run_list, paste0("Mono_L1", 0:2, "_GRCm39"))

run_short_list = paste0("Step1_ML0", 1:9, "_CellByTapeV3_t2.csv")
run_short_list = c(run_short_list, paste0("Step1_ML1", 0:2, "_CellByTapeV3_t2.csv"))

dat_all = NULL
for(kk in 1:length(run_list)){
    print(paste0(kk, "/", length(run_list)))
    dat_i = read.csv(paste0(data_path, "/data/", run_short_list[kk]))[,c(1:4)]
    dat_i$cell_id = paste0(run_list[kk], "_", dat_i$cBC)
    dat_i$barcode = dat_i$mBC
    dat_i = dat_i[dat_i$cell_id %in% pd$cell_id & dat_i$barcode %in% well_barcode$barcode,]
    dat_all = rbind(dat_all, dat_i[,c("cell_id", "barcode", "n_UMI")])
}

well_barcode = well_barcode[well_barcode$barcode %in% dat_all$barcode,]

### 55063/56436 = 98% cells at least one barcode recoved

##########################################
### We generated a cells × barcodes matrix from single-cell profiling and retained the intersecting barcodes shared 
### with the debris-seq table. Each matrix was then normalized by total read counts per cell or per well, respectively. 
### The resulting matrices (100 wells × 893 barcodes × 55,063 cells) were multiplied to obtain assignment scores for each cell-well pair.

dat_1 = well_barcode[,c("well_id", "barcode", "read_count")] %>% dcast(barcode ~ well_id, fill = 0)
rownames(dat_1) = dat_1[,1]; dat_1 = as.matrix(dat_1[,-1])
dat_1 = t(t(dat_1)/apply(dat_1, 2, sum))

dat_2 = dat_all[,c("cell_id", "barcode", "n_UMI")] %>% dcast(barcode ~ cell_id, fill = 0)
rownames(dat_2) = dat_2[,1]; dat_2 = as.matrix(dat_2[,-1])
dat_2 = t(t(dat_2)/apply(dat_2, 2, sum))

### 893 barcodes x 100 wells
### 893 barcodes x 55,063 cells
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
saveRDS(df_sub, paste0(work_path, "/revision/monoclonal_RNAseq/cell_well_cluster_assignment.rds"))
### n = 45,935 cells assigned to 28 clusters; 81% of cells were assigned (45935/56436 = 81%)

well_cluster_list_final = df_sub %>% group_by(well_cluster) %>% tally() %>% 
 filter(n >= 200) %>% arrange(-n) %>% pull(well_cluster)
### n = 26 well clusters

df_sub_200 = df_sub[df_sub$well_cluster %in% well_cluster_list_final,]
### 45,861 cells assigned to 26 well clusters (45861/56436 = 81%)
x = df_sub_200 %>% group_by(well_cluster) %>% tally() 
### 231 - 6938, mean 1764, sd 1175


### making a heatmap
mat_aggr = NULL
for(i in well_cluster_list_final){
    cell_include = df_sub_200 %>% filter(well_cluster == i) %>% pull(cell_id)
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


#########################################################
### Plot the cell type compositions for each well cluster

source("~/work/scripts/utils.R")
work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape"

pd = readRDS(paste0(work_path, "/revision/monoclonal_polyclonal_RNAseq/obj_batchelor_pd.rds"))
cell_assignment = readRDS(paste0(work_path, "/revision/monoclonal_RNAseq/cell_well_cluster_assignment.rds"))

well_cluster_list_final = cell_assignment %>% group_by(well_cluster) %>% tally() %>% 
  filter(n >= 200) %>% arrange(-n) %>% pull(well_cluster)

#### cell type compositions
df = cell_assignment %>% filter(well_cluster %in% well_cluster_list_final) %>% select(well_cluster, cell_id) %>%
    left_join(pd[,c("cell_id", "celltype")], by = "cell_id") %>% group_by(well_cluster, celltype) %>% tally() %>%
    dcast(celltype~well_cluster, fill = 0)
rownames(df) = df[,1]
df = as.matrix(df[,-1])

### some examples:
### hindbrian bias: 3,4,5,21,22
### somite bias: 9,14,15,19
### epiblast bias: 8,16,18,20


### regroup into three categories
group_1 = c("Hindbrain", "Early neurons")
group_2 = c("Epiblast", "Transitional cells")
group_3 = c("Somites")

df_regroup = NULL
df_regroup = rbind(df_regroup, apply(df[rownames(df) %in% group_1,,drop=F], 2, sum))
df_regroup = rbind(df_regroup, apply(df[rownames(df) %in% group_2,,drop=F], 2, sum))
df_regroup = rbind(df_regroup, apply(df[rownames(df) %in% group_3,,drop=F], 2, sum))
df_regroup = rbind(df_regroup, apply(df[!rownames(df) %in% c(group_1, group_2, group_3), ], 2, sum))
rownames(df_regroup) = c("neural", "psc", "meso", "other")

df_regroup = data.frame(t(df_regroup))
df_regroup$well_cluster = gsub("cluster_", "", rownames(df_regroup))

df_regroup <- df_regroup %>%
  mutate(cluster_color = case_when(
    well_cluster %in% c(3,4,5,21,22) ~ "#528ac0",
    well_cluster %in% c(8,16,18,20)  ~ "#52c05a",
    well_cluster %in% c(9,14,15,19)  ~ "#dc87ba",
    TRUE                             ~ "grey50"),
  cluster_group = case_when(
    well_cluster %in% c(3,4,5,21,22) ~ "Hindbrain",
    well_cluster %in% c(8,16,18,20)  ~ "Epiblast",
    well_cluster %in% c(9,14,15,19)  ~ "Somites",
    TRUE                             ~ "Others"))

library(ggtern)
p <- ggtern(df_regroup, aes(x=neural, y=psc, z=meso)) +
  geom_point(aes(color = cluster_color), size=3) +
  geom_text(aes(label = well_cluster), size=2) +
  scale_color_identity()
ggsave("~/share/ggtern_plot.pdf", p, width = 5, height = 5)


### make a boxplot across groups

clone_color_map = c("#52c05a", "#528ac0", "#dc87ba", "grey50")
names(clone_color_map) = c("Epiblast", "Hindbrain", "Somites", "Others")

dat = NULL
for(group_i in c("Epiblast", "Hindbrain", "Somites", "Others")){
    dat_i = df_regroup[df_regroup$cluster_group == group_i, c(1:4)]
    dat_i = 100*dat_i/apply(dat_i, 1, sum)
    dat_i = data.frame(cluster_id = rep(rownames(dat_i), ncol(dat_i)), 
                       frac = c(as.matrix(dat_i)), 
                       celltype_group = rep(colnames(dat_i), each = nrow(dat_i)),
                       categary = group_i)
    dat = rbind(dat, dat_i)
}

dat$celltype_group = factor(dat$celltype_group, levels = c("psc", "neural", "meso", "other"))
dat$categary = factor(dat$categary, levels = c("Epiblast", "Hindbrain", "Somites", "Others"))

p = ggplot(dat, aes(celltype_group, frac, fill = categary)) + geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2) +
    labs(x="", y="% of cells") +
    theme_classic(base_size = 12) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values=clone_color_map) +
    theme(axis.text.x = element_text(color="black", angle = 45), axis.text.y = element_text(color="black")) +
    facet_wrap(~ categary, ncol = 1, scales = "free_y")
ggsave("~/share/box_plot.pdf", p, width = 3, height = 10)


########################################################################
### color the UMAP using cell types or clones for monoclonal gastruloids

source("~/work/scripts/utils.R")
work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape"

pd = readRDS(paste0(work_path, "/revision/monoclonal_polyclonal_RNAseq/obj_batchelor_pd.rds"))
pd = pd[pd$dataset == "monoclonal_RNAseq",]
### monoclonal gastruloids: 56,436 cells

cell_assignment = readRDS(paste0(work_path, "/revision/monoclonal_RNAseq/cell_well_cluster_assignment.rds"))
well_cluster_list_final = cell_assignment %>% group_by(well_cluster) %>% tally() %>% 
  filter(n >= 200) %>% arrange(-n) %>% pull(well_cluster)
cell_assignment = cell_assignment[cell_assignment$well_cluster %in% well_cluster_list_final,]

pd = pd %>% select(UMAP_1, UMAP_2, cell_id, celltype) %>% 
    left_join(cell_assignment[,c("cell_id", "well_cluster")], by = "cell_id")

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
    geom_point(data = pd, aes(x = UMAP_1, y = UMAP_2), color = "grey80",size = 0.1) +
    geom_point(data = pd_sub[pd_sub$clone_group == i,], aes(x = UMAP_1, y = UMAP_2, color = clone_group), size = 0.15) +
    theme_void() +
    theme(legend.position="none") +
    scale_color_manual(values = clone_color_map)
    ggsave(paste0("~/share/data_monoclonal_gastruloids_umap_clones_", i, ".png"), p, width = 3.5, height = 3, dpi = 300)
}

### plot individual clone clusters in umap

well_barcode = read.csv(paste0(work_path, "/revision/debris_seq/well_barcode_set.csv"), header=T)
well_barcode$well_cluster = paste0("cluster_", well_barcode$well_cluster)
well_barcode_content = unique(well_barcode[,c("well_id", "well_cluster")]) 
well_barcode_content_num = well_barcode_content %>% group_by(well_cluster) %>% tally()

plot_list = list()
for (i in mat_other) {
    print(i)
    pd_x = pd_sub %>% filter(well_cluster == i)
    num = well_barcode_content_num$n[well_barcode_content_num$well_cluster == i]
    if (num == 1){num_x = paste0(num, " gastruloid")} else {num_x = paste0(num, " gastruloids")}
    plot_list[[i]] = ggplot() +
        geom_point(data = pd, aes(x = UMAP_1, y = UMAP_2), color = "grey80",size = 0.1) +
        geom_point(data = pd_x, aes(x = UMAP_1, y = UMAP_2, color = clone_group), size = 0.15) +
        theme_void() +
        scale_color_manual(values = clone_color_map) +
        labs(title = paste0(i, ": ", nrow(pd_x), " cells from ", num_x)) +
        theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 9))
}
combined_plot <- wrap_plots(plot_list, nrow = 2)

ggsave(paste0("~/share/data_monoclonal_gastruloids_umap_clones_", "Others", ".split_individual.png"), combined_plot, width = 21, height = 6, dpi = 300)



##################### BACKUP #########################




