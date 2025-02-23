
##################################
### Performing pseudobulk analysis
### Chengxiang Qiu
### Feb-20, 2025


###################################################################################################
### Step-1: Aggregating cells from each individual gastruloids, followed by performing PCA analysis

source("help_script.R")
work_path = "Your_work_path"

pd = readRDS(paste0(work_path, "/obj_processed_pd.rds"))
obj = readRDS(paste0(work_path, "/obj_processed.rds"))
count = GetAssayData(obj, slot = "counts")

### keep gastruloids with at least 100 cells captured (n=121)
cell_assign = readRDS(paste0(work_path, "/tape_barcode_3/cell_assign.rds"))

df = cell_assign %>%
    group_by(well, celltype) %>% tally() %>%
    dcast(celltype ~ well, fill = 0)

rownames(df) = as.vector(df[[1]])
df = df[,-1]


well_list = unique(cell_assign$well)
count_aggr = NULL
for(i in 1:length(well_list)){
    well_i = well_list[i]
    print(well_i)
    
    count_sub = count[,colnames(count) %in% as.vector(cell_assign$cell_id[cell_assign$well == well_i])]
    count_aggr = cbind(count_aggr,
                       Matrix::rowSums(count_sub))
}
colnames(count_aggr) = well_list

### Here, we are only using genes which overlapped with the data from the "tree of trees" experiment
count_mGASv5 = readRDS("../data_mGASv5/count.rds")
count_aggr = count_aggr[rownames(count_aggr) %in% rownames(count_mGASv5),]

obj = CreateSeuratObject(count_aggr)
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2500)
genes_include = VariableFeatures(obj)
saveRDS(obj, paste0(work_path, "/tape_barcode_3/pca/obj_pseudobulk.rds"))

cds = new_cell_data_set(as.matrix(count_aggr))
saveRDS(cds, paste0(work_path, "/tape_barcode_3/pca/cds_pseudobulk.rds"))

set.seed(2016)
FM = monocle3:::normalize_expr_data(cds, 
                                    norm_method = "log",
                                    pseudo_count = 1)

FM = FM[genes_include,]

num_dim = 10
scaling = TRUE
set.seed(2016)
irlba_res = my_sparse_prcomp_irlba(Matrix::t(FM), 
                                   n = min(num_dim, min(dim(FM)) - 1), 
                                   center = scaling, 
                                   scale. = scaling)
preproc_res = irlba_res$x
row.names(preproc_res) = colnames(cds)

saveRDS(irlba_res, paste0(work_path, "/tape_barcode_3/pca/PCA_pseudobulk.rds"))

prop_var_expl = irlba_res$sdev^2/sum(irlba_res$sdev^2)
print(prop_var_expl)
### 35.7%, 16.5%, 10.7%



######################################
### Step-2: Making scatter plot on PCA

df_x = data.frame(t(df))
df_y = data.frame(sample_id = rownames(preproc_res),
                  PC_1 = preproc_res[,1],
                  PC_2 = preproc_res[,2],
                  PC_3 = preproc_res[,3])
df_y = df_y[rownames(df_x),]
print(sum(rownames(df_x) == df_y$sample_id))
df_x$PC_1 = as.vector(df_y$PC_1)
df_x$PC_2 = as.vector(df_y$PC_2)
df_x$PC_3 = as.vector(df_y$PC_3)

saveRDS(df_x, paste0(work_path, "/tape_barcode_3/pca/cds_pseudobulk_pca.rds"))

gastruloid_celltype_color_code = c("Cardiac.mesoderm" = "#dc9436",
                                   "Definitive.endoderm" = "#5d71db",
                                   "Early.neurons" = "#85b937",
                                   "Endothelial.cells" = "#904bb8",
                                   "Epiblast" = "#52c05a",
                                   "Fibroblasts" = "#cc79dc",
                                   "Floor.plate" = "#c2b73d",
                                   "Node.like.cells" = "#7e72b7",
                                   "Hindbrain" = "#8e8927",
                                   "Anterior.mesendoderm" = "#5a9bd5",
                                   "Mesodermal.progenitors" = "#dd5732",
                                   "MHB" = "#3fc1bf",
                                   "Motor.neurons" = "#dd4168",
                                   "NMPs" = "#5fc08c",
                                   "Notochord" = "#a14c78",
                                   "Extraembryonic.endoderm" = "#37835d",
                                   "Somites" = "#dc87ba",
                                   "Spinal.cord" = "#a4b46c",
                                   "Transitional.cells" = "#ac4c55")


p = ggplot() + 
    geom_scatterpie(data = df_x, aes(x=PC_1, y=PC_2), cols=colnames(df_x)[1:16], color=NA, alpha=.9) + 
    coord_equal() +
    theme_classic(base_size = 10) +
    labs(x = "PC_1 (35.7%)", y = "PC_2 (16.5%)") + 
    scale_fill_manual(values=gastruloid_celltype_color_code) +
    # theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
pdf(paste0(work_path, "/tape_barcode_3/pca/PCA.pdf"), 12, 10)
print(p)
dev.off()

p = ggplot() + 
    geom_scatterpie(data = df_x, aes(x=PC_1, y=PC_2), cols=colnames(df_x)[1:16], color=NA, alpha=.9) + 
    coord_equal() +
    theme_classic(base_size = 10) +
    labs(x = "PC_1 (35.7%)", y = "PC_2 (16.5%)") + 
    scale_fill_manual(values=gastruloid_celltype_color_code) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    xlim(-44, 45) +
    ylim(-60, 25) 
pdf(paste0(work_path, "/tape_barcode_3/pca/PCA_nolegend.pdf"), 10, 10)
print(p)
dev.off()

p = ggplot() + 
    #geom_scatterpie(data = df_x, aes(x=PC_1, y=PC_2), cols=colnames(df_x)[1:16], color=NA, alpha=.9) + 
    #coord_equal() +
    geom_text(data = df_x, aes(x = PC_1, y = PC_2, label = rownames(df_x))) +
    theme_classic(base_size = 10) +
    labs(x = "PC_1 (35.7%)", y = "PC_2 (16.5%)") + 
   # scale_fill_manual(values=gastruloid_celltype_color_code) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    xlim(-44, 45) +
    ylim(-60, 25) 
pdf(paste0(work_path, "/tape_barcode_3/pca/PCA_label.pdf"), 10, 10)
print(p)
dev.off()

df_pct = t(df)/apply(df,2,sum)
df_pct = df_pct[rownames(df_x),]
corr = cor(df_pct, df_x[,c("PC_1","PC_2","PC_3")])
print(corr)


##########################################################################
### Step-3: Making line plot to show cell-type composition changes over PC

dat_1 = data.frame(pct = c(df_pct[,9], df_pct[,14], df_pct[,15]),
                   celltype = rep(colnames(df_pct)[c(9,14,15)], each = nrow(df_pct)),
                   PC_1 = rep(df_x[,"PC_1"], 3))
dat_1$pct = 100 * dat_1$pct

p = ggplot(dat_1, aes(x=PC_1, y=pct, color=celltype)) +
    geom_smooth(method = loess, se = FALSE) +
    theme_classic(base_size = 10) +
    labs(x = "PC_1 (35.7%)", y = "% of cells") + 
    scale_color_manual(values=gastruloid_celltype_color_code) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    xlim(-44, 45) 
pdf(paste0(work_path, "/tape_barcode_3/pca/PCA_corr_PC1.pdf"), 10, 2)
print(p)
dev.off()


dat_2 = data.frame(pct = c(df_pct[,6]),
                   celltype = rep(colnames(df_pct)[c(6)], each = nrow(df_pct)),
                   PC_2 = rep(df_x[,"PC_2"], 1))
dat_2$pct = 100 * dat_2$pct

p = ggplot(dat_2, aes(x=PC_2, y=pct, color=celltype)) +
    geom_smooth(method = loess, se = FALSE) +
    theme_classic(base_size = 10) +
    labs(x = "PC_2 (16.5%)", y = "% of cells") + 
    scale_color_manual(values=gastruloid_celltype_color_code) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    xlim(-60, 25) 
pdf(paste0(work_path, "/tape_barcode_3/pca/PCA_corr_PC2.pdf"), 10, 2)
print(p)
dev.off()


###################################################################
### Step-4: Performing k-means clustering based on cell-prop matrix

df_pct = t(t(df)/apply(df,2,sum))
df_pct = df_pct[,rownames(df_x)]
cluster_res = kmeans(t(df_pct), centers = 4, nstart = 25)$cluster
cluster_res = cluster_res[rownames(df_x)]

df_x$kmeans_cluster = factor(as.vector(cluster_res))

saveRDS(df_x, paste0(work_path, "/tape_barcode_3/pca/cds_pseudobulk_pca.rds"))

gastruloid_cluster = rep("ND", nrow(df_x))
gastruloid_cluster[df_x$kmeans_cluster == 2] = "INT"
gastruloid_cluster[df_x$kmeans_cluster == 3] = "MESO"
gastruloid_cluster[df_x$kmeans_cluster == 4] = "NEURO"
df_x$gastruloid_cluster = as.vector(gastruloid_cluster)

highlight_list = c("P3-B10", "P3-G7", "P3-C6", "P3-C11", "P4-D9", "P3-E12",
                   "P3-G5", "P3-F3", "P3-C8", "P4-F10", "P4-D8", "P3-D5")

p = ggplot() +
    geom_point(data = df_x, aes(PC_1, PC_2, color = gastruloid_cluster), size = 3.5) + 
    geom_point(data = df_x[rownames(df_x) %in% highlight_list,], aes(PC_1, PC_2), color = "black", size = 5) + 
    geom_point(data = df_x[rownames(df_x) %in% highlight_list,], aes(PC_1, PC_2, color = gastruloid_cluster), size = 3.5) + 
    theme_void() +
    theme(legend.position="none") +
    scale_color_manual(values=c("ND" = "#52c05a", "INT" = "#c1d088", "MESO" = "#dc87ba", "NEURO" = "#6e681f")) +
    xlim(-44, 45) +
    ylim(-60, 25) 
pdf(paste0(work_path, "/tape_barcode_3/pca/PCA_cluster.pdf"), 5, 5)
print(p)
dev.off()


###################################################
### Step-5: cluster compositions for each cell type

### keep gastruloids with at least 100 cells captured (n=121)
cell_assign = readRDS(paste0(work_path, "/tape_barcode_3/cell_assign.rds"))

df = cell_assign %>%
    group_by(well, celltype) %>% tally() %>%
    dcast(celltype ~ well, fill = 0)

rownames(df) = as.vector(df[[1]])
df = df[,-1]

df_x = readRDS(paste0(work_path, "/tape_barcode_3/pca/cds_pseudobulk_pca.rds"))
sum(colnames(df) == rownames(df_x))
gastruloid_cluster = rep("ND", nrow(df_x))
gastruloid_cluster[df_x$kmeans_cluster == 2] = "INT"
gastruloid_cluster[df_x$kmeans_cluster == 3] = "MESO"
gastruloid_cluster[df_x$kmeans_cluster == 4] = "NEURO"
df_x$gastruloid_cluster = as.vector(gastruloid_cluster)

df_pct = t(df)/apply(df,2,sum)
df_pct = df_pct[rownames(df_x),]

dat = data.frame(frac = 100*c(as.matrix(df_pct)),
                 celltype = rep(colnames(df_pct), each = nrow(df_pct)),
                 gastruloid_cluster = rep(df_x$gastruloid_cluster, times = ncol(df_pct)))
dat = dat %>% group_by(celltype, gastruloid_cluster) %>% summarize(mean_frac = mean(frac))

dat$gastruloid_cluster = factor(dat$gastruloid_cluster, levels = c("ND","MESO","INT","NEURO"))
dat$celltype = factor(dat$celltype, levels = c("Epiblast", "Node-like cells", "Definitive endoderm", "Notochord",
                                               "Transitional cells", "NMPs", "Mesodermal progenitors",
                                               "Spinal cord", "Somites", "Anterior mesendoderm",
                                               "Cardiac mesoderm", "Endothelial cells", "Early neurons",
                                               "Hindbrain", "Floor plate", "Extraembryonic endoderm"))

p = dat %>%
    ggplot(aes(x=gastruloid_cluster, y=mean_frac, fill=celltype)) + geom_bar(stat='identity') + 
    facet_wrap(~celltype, scales = "free_y", ncol = 4) + 
    labs(x='',y='% of cells') +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    scale_fill_manual(values = gastruloid_celltype_color_code) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black", angle = 45, hjust=1), axis.text.y = element_text(color="black")) 
pdf(paste0(work_path, "/tape_barcode_3/pca/Cell_compositions.pdf"), 5, 5)
print(p)
dev.off()




######################################################################
### Step-6: cell number of individual gastruloids across four clusters

### keep gastruloids with at least 100 cells captured (n=121)
cell_assign = readRDS(paste0(work_path, "/tape_barcode_3/cell_assign.rds"))
cell_num = cell_assign %>% group_by(well) %>% tally()

df_x = readRDS(paste0(work_path, "/tape_barcode_3/pca/cds_pseudobulk_pca.rds"))
gastruloid_cluster = rep("ND", nrow(df_x))
gastruloid_cluster[df_x$kmeans_cluster == 2] = "INT"
gastruloid_cluster[df_x$kmeans_cluster == 3] = "MESO"
gastruloid_cluster[df_x$kmeans_cluster == 4] = "NEURO"
df_x$gastruloid_cluster = as.vector(gastruloid_cluster)

df = data.frame(well = rownames(df_x), gastruloid_cluster = as.vector(df_x$gastruloid_cluster)) %>%
    left_join(cell_num, by = "well") %>% mutate(log2_num = log2(n))
df$gastruloid_cluster = factor(df$gastruloid_cluster, levels = c("ND","MESO","INT","NEURO"))

p = ggplot(df, aes(gastruloid_cluster, log2_num, fill = gastruloid_cluster)) + geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2) +
    labs(x="", y="Log2 (cell #) of individual gastruloids") +
    theme_classic(base_size = 12) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values=c("ND" = "#52c05a", "INT" = "#c1d088", "MESO" = "#dc87ba", "NEURO" = "#6e681f")) +
    theme(axis.text.x = element_text(color="black", angle = 45), axis.text.y = element_text(color="black")) 

pdf(paste0(work_path, "/tape_barcode_3/cell_num_gastruloids.pdf"), 5, 4)
print(p)
dev.off()

summary(aov(log2_num ~ factor(gastruloid_cluster), data = df))

df_x$log2_num = as.vector(df$log2_num)
df_x %>% 
    ggplot(aes(PC_2, log2_num)) + geom_point(aes(color=gastruloid_cluster), size=3) + geom_smooth(method = "lm", se = FALSE) +
    labs(x="PC_2", y="log2_cell_num", title="") +
    theme_classic(base_size = 10) +
    scale_color_manual(values=c("ND" = "#52c05a", "INT" = "#c1d088", "MESO" = "#dc87ba", "NEURO" = "#6e681f")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    ggsave("~/share/PC2_log2_num.pdf")







