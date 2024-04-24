
### After assigning cells to the gastruloids, we are going to perform pseudobulk analysis by PCA

work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape/data_sci"
source("~/work/scripts/tome/utils.R")
library(reshape2)

cell_assign = readRDS(paste0(work_path, "/tape_barcode_2/cell_assigning_res.rds"))
pd = readRDS(paste0(work_path, "/obj_processed_pd.rds"))
obj = readRDS(paste0(work_path, "/obj_processed.rds"))
count = GetAssayData(obj, slot = "counts")

df = cell_assign %>% left_join(pd %>% select(cell = sample, celltype), by = "cell") %>%
    group_by(sample, celltype) %>% tally() %>%
    dcast(celltype ~ sample)

rownames(df) = as.vector(df[[1]])
df = df[,-1]
df[is.na(df)] = 0

### keep gastruloids with at least 100 cells captured (n=134)
df_colsum = apply(df,2,sum)
df = df[,df_colsum > 100]

cell_assign_sub = cell_assign[cell_assign$sample %in% colnames(df),]
sample_list = colnames(df)
count_aggr = NULL
for(i in 1:length(sample_list)){
    sample_i = sample_list[i]
    print(sample_i)
    
    count_sub = count[,colnames(count) %in% as.vector(cell_assign_sub$cell[cell_assign_sub$sample == sample_i])]
    count_aggr = cbind(count_aggr,
                       Matrix::rowSums(count_sub))
}
colnames(count_aggr) = sample_list

obj = CreateSeuratObject(count_aggr)
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2500)
genes_include = VariableFeatures(obj)

cds = new_cell_data_set(as.matrix(count_aggr))
saveRDS(cds, paste0(work_path, "/tape_barcode_2/cds_pseudobulk.rds"))

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

prop_var_expl = irlba_res$sdev^2/sum(irlba_res$sdev^2)
print(prop_var_expl)

### making UMAP?
set.seed(2016)
emb_train_umap = uwot::umap(as.matrix(preproc_res), 
                            n_components = 2,
                            n_neighbors = 5,
                            min_dist = 0.3,
                            metric = "cosine",
                            fast_sgd = FALSE,
                            nn_method = "annoy",
                            ret_model = TRUE,
                            n_threads = 1,
                            verbose = TRUE)

set.seed(2016)
umap_coor = uwot::umap_transform(as.matrix(preproc_res),
                                     emb_train_umap)


### making scatter plot
library(ggplot2)
library(scatterpie)
df_x = data.frame(t(df))
df_y = data.frame(sample_id = rownames(preproc_res),
                  PC_1 = preproc_res[,1],
                  PC_2 = preproc_res[,2],
                  PC_3 = preproc_res[,3],
                  UMAP_1 = umap_coor[,1],
                  UMAP_2 = umap_coor[,2])
print(sum(rownames(df_x) == df_y$sample_id))
df_x$PC_1 = as.vector(df_y$PC_1)
df_x$PC_2 = as.vector(df_y$PC_2)
df_x$PC_3 = as.vector(df_y$PC_3)
df_x$UMAP_1 = as.vector(df_y$UMAP_1)
df_x$UMAP_2 = as.vector(df_y$UMAP_2)

saveRDS(list(df = df,
             df_x = df_x), paste0(work_path, "/tape_barcode_2/cds_pseudobulk_pca.rds"))

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
    labs(x = "PC_1 (38.2%)", y = "PC_2 (15.7%)") + 
    scale_fill_manual(values=gastruloid_celltype_color_code) +
    # theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
pdf("~/share/gastruloid_PCA.pdf", 12, 10)
print(p)
dev.off()

p = ggplot() + 
    geom_scatterpie(data = df_x, aes(x=PC_1, y=PC_2), cols=colnames(df_x)[1:16], color=NA, alpha=.9) + 
    coord_equal() +
    theme_classic(base_size = 10) +
    labs(x = "PC_1 (38.2%)", y = "PC_2 (15.7%)") + 
    scale_fill_manual(values=gastruloid_celltype_color_code) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    xlim(-40, 48) + 
    ylim(-62, 32)
pdf("~/share/gastruloid_PCA_nolegend.pdf", 10, 10)
print(p)
dev.off()



df_pct = t(df)/apply(df,2,sum)
corr = cor(df_pct, df_x[,c("PC_1","PC_2","PC_3")])
print(corr)

### making line plot to show cell-type composition changes over PC

dat_1 = data.frame(pct = c(df_pct[,9], df_pct[,14], df_pct[,15]),
                   celltype = rep(colnames(df_pct)[c(9,14,15)], each = nrow(df_pct)),
                   PC_1 = rep(df_x[,"PC_1"], 3))
dat_1$pct = 100 * dat_1$pct

p = ggplot(dat_1, aes(x=PC_1, y=pct, color=celltype)) +
    geom_smooth(method = loess, se = FALSE) +
    theme_classic(base_size = 10) +
    labs(x = "PC_1 (38.2%)", y = "% of cells") + 
    scale_color_manual(values=gastruloid_celltype_color_code) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    xlim(-40, 48) 
pdf("~/share/gastruloid_PCA_scatter_plot_PC1.pdf", 10, 2)
print(p)
dev.off()


dat_2 = data.frame(pct = c(df_pct[,6]),
                   celltype = rep(colnames(df_pct)[c(6)], each = nrow(df_pct)),
                   PC_2 = rep(df_x[,"PC_2"], 1))
dat_2$pct = 100 * dat_2$pct

p = ggplot(dat_2, aes(x=PC_2, y=pct, color=celltype)) +
    geom_smooth(method = loess, se = FALSE) +
    theme_classic(base_size = 10) +
    labs(x = "PC_2 (15.7%)", y = "% of cells") + 
    scale_color_manual(values=gastruloid_celltype_color_code) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    xlim(-62, 32) 
pdf("~/share/gastruloid_PCA_scatter_plot_PC2.pdf", 10, 2)
print(p)
dev.off()



### can we do k-means clustering based on cell-prop matrix


work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape/data_sci"
source("~/work/scripts/utils.R")


df_pct = t(t(df)/apply(df,2,sum))
cluster_res = kmeans(t(df_pct), centers = 4, nstart = 25)$cluster
cluster_res = cluster_res[rownames(df_x)]

df_x$kmeans_cluster = factor(as.vector(cluster_res))

saveRDS(list(df = df,
             df_x = df_x), paste0(work_path, "/tape_barcode_2/cds_pseudobulk_pca.rds"))


p = ggplot() +
    geom_point(data = df_x, aes(PC_1, PC_2, color = kmeans_cluster), size = 3.5) + 
    theme_void() +
    theme(legend.position="none") +
    scale_color_manual(values=c("1" = "#b98d3e", "2" = "#9970c1", "3" = "#64a860", "4" = "#cc545e")) +
    xlim(-40, 48) + 
    ylim(-62, 32)
pdf("gastruloid_PCA_kmeans_clustering.pdf", 5, 5)
print(p)
dev.off()

df_pct = t(t(df)/apply(df,2,sum))
df_pct = df_pct[,row.names(df_x)]
df_pct = 100*df_pct

dat_x = data.frame(pct = c(df_pct),
                   celltype = rep(rownames(df_pct), times = ncol(df_pct)),
                   kmeans_cluster = rep(as.vector(df_x$kmeans_cluster), each = nrow(df_pct)))
dat_x$celltype = factor(dat_x$celltype, levels = names(gastruloid_celltype_color_code)[names(gastruloid_celltype_color_code) %in% dat_x$celltype])


p = ggplot(dat_x, aes(factor(celltype), pct, fill = celltype)) + 
    geom_boxplot() +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    scale_fill_manual(values=gastruloid_celltype_color_code) +
    facet_grid(rows = vars(kmeans_cluster)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
pdf("gastruloid_PCA_kmeans_clustering_pct.pdf", 8, 3)
print(p)
dev.off()

p = ggplot(dat_x, aes(factor(celltype), log2(pct), fill = celltype)) + 
    geom_boxplot() +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    scale_fill_manual(values=gastruloid_celltype_color_code) +
    facet_grid(rows = vars(kmeans_cluster)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
pdf("~/share/gastruloid_PCA_kmeans_clustering_pct_3.pdf", 3, 8)
print(p)
dev.off()

### cluster-compositions for each cell type

dat = melt(as.matrix(df), value.name = "cell_num"); colnames(dat) = c("celltype", "gastruloid", "cell_num")
dat_x = data.frame(gastruloid = rownames(df_x), kmeans_cluster = df_x$kmeans_cluster)
dat = dat %>% left_join(dat_x, by = "gastruloid") 
dat_1 = dat %>% group_by(celltype, kmeans_cluster) %>% summarize(cell_num_x = mean(cell_num)) 
dat_2 = dat_1 %>%
    left_join(dat_1 %>% group_by(celltype) %>% summarize(total_cell_num = sum(cell_num_x)), by = "celltype") %>% 
    mutate(percentage = cell_num_x/total_cell_num)

dat_3 = dat_2 %>% filter(kmeans_cluster == 1) %>% arrange(percentage)

dat_2$celltype = factor(dat_2$celltype, levels = as.vector(dat_3$celltype))

p = ggplot(dat_2, aes(x = celltype, y = percentage*100, fill = kmeans_cluster)) +
    geom_bar(stat="identity", width = 1) +
    labs(x = "", y = "") +
    theme_classic(base_size = 12) +
    scale_fill_manual(values=c("1" = "#b98d3e", "2" = "#9970c1", "3" = "#64a860", "4" = "#cc545e")) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black", angle = 90, hjust = 1), axis.text.y = element_text(color="black")) 
pdf("~/share/gastruloid_PCA_kmeans_clustering_pct_2.pdf")
print(p)
dev.off()



##################### BACKUP #######################################################


### After assiging cells to indivudal gastruloids, we are looking for the heterogenity of cell-type compositions

work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape/data_sci"
source("~/work/scripts/tome/utils.R")
library(reshape2)

cell_assign = readRDS(paste0(work_path, "/tape_barcode_2/cell_assigning_res.rds"))
pd = readRDS(paste0(work_path, "/obj_processed_pd.rds"))
#obj = readRDS(paste0(work_path, "/obj_processed.rds"))

df = cell_assign %>% left_join(pd %>% select(cell = sample, celltype), by = "cell") %>%
    group_by(sample, celltype) %>% tally() %>%
    dcast(celltype ~ sample)

rownames(df) = as.vector(df[[1]])
df = df[,-1]
df[is.na(df)] = 0

### keep gastruloids with at least 100 cells captured (n=134)
df_colsum = apply(df,2,sum)
df = df[,df_colsum > 100]

df_pct = t(t(df)/apply(df,2,sum))

cds = new_cell_data_set(as.matrix(df))

set.seed(2016)
FM = monocle3:::normalize_expr_data(cds, 
                                    norm_method = "log")

num_dim = 10
scaling = TRUE
set.seed(2016)
irlba_res = my_sparse_prcomp_irlba(Matrix::t(FM), 
                                   n = min(num_dim, min(dim(FM)) - 1), 
                                   center = scaling, 
                                   scale. = scaling)
preproc_res = irlba_res$x
row.names(preproc_res) = colnames(cds)

prop_var_expl = irlba_res$sdev^2/sum(irlba_res$sdev^2)
print(prop_var_expl)

### making scatter plot
library(ggplot2)
library(scatterpie)
df_x = data.frame(t(df))
df_y = data.frame(sample_id = rownames(preproc_res),
                  PC_1 = preproc_res[,1],
                  PC_2 = preproc_res[,2],
                  PC_3 = preproc_res[,3])
print(sum(rownames(df_x) == df_y$sample_id))
df_x$PC_1 = as.vector(df_y$PC_1)
df_x$PC_2 = as.vector(df_y$PC_2)
df_x$PC_3 = as.vector(df_y$PC_3)

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
    labs(x = "PC_1 (37.0%)", y = "PC_2 (20.5%)") + 
    scale_fill_manual(values=gastruloid_celltype_color_code) +
    # theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
pdf("~/share/gastruloid_PCA_2.pdf", 12, 10)
print(p)
dev.off()





dat_1 = data.frame(pct = c(df_pct),
                   celltype = rep(colnames(df_pct), each = nrow(df_pct)),
                   PC_1 = rep(df_x[,"PC_1"], 16))
dat_1$pct = 100 * dat_1$pct

p = ggplot(dat_1, aes(x=PC_1, y=pct, color=celltype)) +
    geom_smooth(method = loess, se = FALSE) +
    theme_classic(base_size = 10) +
    labs(x = "PC_1 (38.2%)", y = "% of cells") + 
    scale_color_manual(values=gastruloid_celltype_color_code) +
    #theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    xlim(-40, 48) 
pdf("~/share/gastruloid_PCA_scatter_plot_PC1_all.pdf", 10, 5)
print(p)
dev.off()


dat_2 = data.frame(pct = c(df_pct),
                   celltype = rep(colnames(df_pct), each = nrow(df_pct)),
                   PC_2 = rep(df_x[,"PC_2"], 1))
dat_2$pct = 100 * dat_2$pct

p = ggplot(dat_2, aes(x=PC_2, y=pct, color=celltype)) +
    geom_smooth(method = loess, se = FALSE) +
    theme_classic(base_size = 10) +
    labs(x = "PC_2 (15.7%)", y = "% of cells") + 
    scale_color_manual(values=gastruloid_celltype_color_code) +
    #theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    xlim(-62, 32) 
pdf("~/share/gastruloid_PCA_scatter_plot_PC2_all.pdf", 10, 5)
print(p)
dev.off()



