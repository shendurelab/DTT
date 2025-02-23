
#######################################################
### Integration with other in vivo or in vitro datasets
### Chengxiang Qiu
### Feb-20, 2025

######################################
### Step-1: Integrating seven datasets

source("help_script.R")
mouse_gene = read.table("mouse39-samchoiTAPE.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

library(future)
library(future.apply)
plan("multicore", workers = 4)
options(future.globals.maxSize = 90000 * 1024^2)

### 1st dataset - Our 144 gastruloids dataset, generated using sci-RNA-seq3
dat = readRDS(paste0(work_path, "/obj_processed.rds"))
pd = readRDS(paste0(work_path, "/obj_processed_pd.rds"))
count_1 = GetAssayData(dat, slot = "counts")
pd_1 = data.frame(pd)[,c("plate","celltype")]
pd_1$group = "gastruloid"

### downsampling to 100K cells
keep = sample(1:nrow(pd_1), 100000)
count_1 = count_1[,keep]
pd_1 = pd_1[keep,]

fd_1 = mouse_gene[rownames(count_1),]
fd_1$rowSum = Matrix::rowSums(count_1)
fd_1_x = fd_1 %>% group_by(gene_short_name) %>% slice_max(order_by = rowSum, n = 1, with_ties = F)
count_1 = count_1[as.vector(fd_1_x$gene_ID),]
rownames(count_1) = as.vector(fd_1_x$gene_short_name)


### 2nd dataset - Suppinger's dataset
### reference: https://pubmed.ncbi.nlm.nih.gov/37209681/
count_2 = Matrix::readMM(paste0(data_path, "/gastruloid_Liberali/GSE229513_UMI_counts.mtx"))
rownames(count_2) = as.vector(read.table(paste0(data_path, "/gastruloid_Liberali/GSE229513_genes.tsv"), header=T)$gene)
colnames(count_2) = as.vector(read.table(paste0(data_path, "/gastruloid_Liberali/GSE229513_barcodes.tsv"), sep="\t")$V1)
pd_2 = readRDS(paste0(data_path, "/gastruloid_Liberali/df_cell.rds"))
pd_2 = pd_2[,c(1:4, 9, 13)]
pd_2$celltype = pd_2$celltypeannotation
pd_2$plate = pd_2$batch
pd_2 = pd_2[,c("plate","timepoints","celltype")]
pd_2$group = "Liberali"


### 3rd dataset - TLS data
### reference: https://pubmed.ncbi.nlm.nih.gov/33303587/
count_TLS_96h = Read10X(paste0(data_path, "/gastruloid_TLS/TLS_96h"), gene.column = 1, strip.suffix = T)
count_TLS_108h = Read10X(paste0(data_path, "/gastruloid_TLS/TLS_108h"), gene.column = 1, strip.suffix = T)
count_TLS_120h = Read10X(paste0(data_path, "/gastruloid_TLS/TLS_120h"), gene.column = 1, strip.suffix = T)

meta_data = read.table(paste0(data_path, "/gastruloid_TLS/TLS_meta_data.tsv"), as.is=T, header=T)
meta_data$sample = paste0(meta_data$TP, "_", meta_data$BC)
rownames(meta_data) = as.vector(meta_data$sample)

colnames(count_TLS_96h) = paste0("TLS_96h_", colnames(count_TLS_96h))
count_TLS_96h = count_TLS_96h[,colnames(count_TLS_96h) %in% meta_data$sample]
colnames(count_TLS_108h) = paste0("TLS_108h_", colnames(count_TLS_108h))
count_TLS_108h = count_TLS_108h[,colnames(count_TLS_108h) %in% meta_data$sample]
colnames(count_TLS_120h) = paste0("TLS_120h_", colnames(count_TLS_120h))
count_TLS_120h = count_TLS_120h[,colnames(count_TLS_120h) %in% meta_data$sample]

count_3 = cbind(count_TLS_96h, count_TLS_108h, count_TLS_120h)
pd_3 = meta_data[colnames(count_3),]
pd_3$celltype = pd_3$cell_state
pd_3$plate = pd_3$TP
pd_3$group = "TLS"
pd_3 = pd_3[,c("plate","celltype","group")]

fd_3 = mouse_gene[rownames(count_3),]
fd_3$rowSum = Matrix::rowSums(count_3)
fd_3_x = fd_3 %>% filter(!is.na(gene_short_name)) %>% group_by(gene_short_name) %>% slice_max(order_by = rowSum, n = 1, with_ties = F)
count_3 = count_3[as.vector(fd_3_x$gene_ID),]
rownames(count_3) = as.vector(fd_3_x$gene_short_name)


### 4th dataset - natural embryos during organogenesis
### reference: https://pubmed.ncbi.nlm.nih.gov/38355799/
pd_all = readRDS(paste0(data_path_3, "/mtx/adata_scale.obs.rds"))

pd_x = readRDS("./obj_integrated_embryo_pd_old.rds")
pd = pd_all[pd_all$cell_id %in% rownames(pd_x),]
rownames(pd) = as.vector(pd$cell_id)

embryo_list = as.vector(names(table(pd$embryo_id)))

count_4 = NULL
for(j in embryo_list){
    
    print(paste0(j,"/",length(embryo_list)))
    count_j = readRDS(paste0(data_path_3, "/embryo/", j, "_gene_count.rds"))
    keep = colnames(count_j) %in% rownames(pd)
    count_4 = cbind(count_4, count_j[,keep,drop=FALSE])
    rm(count_j)
}
pd_4 = pd[colnames(count_4),]
pd_4$group = "jax"
pd_4$plate = pd_4$sequencing_batch
pd_4$celltype = pd_4$major_trajectory
pd_4$timepoints = pd_4$day
pd_4 = pd_4[,c("plate","celltype","timepoints","group","celltype_update")]

fd_4 = mouse_gene[rownames(count_4),]
fd_4$rowSum = Matrix::rowSums(count_4)
fd_4_x = fd_4 %>% filter(!is.na(gene_short_name)) %>% group_by(gene_short_name) %>% slice_max(order_by = rowSum, n = 1, with_ties = F)
count_4 = count_4[as.vector(fd_4_x$gene_ID),]
rownames(count_4) = as.vector(fd_4_x$gene_short_name)

### downsampling to 100K cells
keep = sample(1:nrow(pd_1), 100000)
count_4 = count_4[,keep]
pd_4 = pd_4[keep,]


### 5th dataset - natural embryos during gastrulation
### reference: https://pubmed.ncbi.nlm.nih.gov/30787436/
obj_5 = readRDS("./Pijuan.rds")
pd_5 = data.frame(obj_5[[]])
pd_5$celltype = pd_5$pre_celltype
pd_5$plate = paste0("batch_",pd_5$group)
pd_5$group = "pijuan"
pd_5$timepoints = pd_5$day
pd_5 = pd_5[,c("plate","celltype","timepoints","group")]
count_5 = GetAssayData(obj_5, slot = "counts")
rm(obj_5)

fd_5 = mouse_gene[rownames(count_5),]
fd_5$rowSum = Matrix::rowSums(count_5)
fd_5_x = fd_5 %>% filter(!is.na(gene_short_name)) %>% group_by(gene_short_name) %>% slice_max(order_by = rowSum, n = 1, with_ties = F)
count_5 = count_5[as.vector(fd_5_x$gene_ID),]
rownames(count_5) = as.vector(fd_5_x$gene_short_name)


### 6th dataset - Rosen's gastruloid
### reference: https://www.biorxiv.org/content/10.1101/2022.09.27.509783v1
obj_6 = readRDS(paste0(data_path, "/gastruloid_Rosen/GSE212050_seurat_final.rds"))
fd_6 = read.table(paste0(data_path, "/gastruloid_Rosen/GSE212050_feature_metadata_final.txt"), header=T, as.is=T, sep=",")[,c(1,2)]
names(fd_6) = c("gene_ID","gene_short_name")
pd_6 = data.frame(obj_6[[]])[,c("experiment", "cluster","timepoint")]
pd_6$group = "Rosen"
names(pd_6) = c("plate","celltype","timepoints","group")

count_6 = GetAssayData(obj_6, slot = "counts")
rm(obj_6)

sum(fd_6$gene_ID == rownames(count_6))
fd_6$rowSum = Matrix::rowSums(count_6)
fd_6_x = fd_6 %>% filter(!is.na(gene_short_name)) %>% group_by(gene_short_name) %>% slice_max(order_by = rowSum, n = 1, with_ties = F)
count_6 = count_6[as.vector(fd_6_x$gene_ID),]
rownames(count_6) = as.vector(fd_6_x$gene_short_name)



### 7th dataset - Van's gastruloid
### reference: https://pubmed.ncbi.nlm.nih.gov/32076263/
dat = readRDS(paste0(data_path, "/gastruloid_van/dat.rds"))
count_7 = dat[['count']]
fd_7 = dat[['fd']]
pd_7 = dat[['pd']]
names(fd_7) = c("gene_ID","gene_short_name")
pd_7 = pd_7[,c("batch", "celltye")]
pd_7$group = "Van"
names(pd_7) = c("plate","celltype","group")

sum(fd_7$gene_ID == rownames(count_7))
fd_7$rowSum = Matrix::rowSums(count_7)
fd_7_x = fd_7 %>% filter(!is.na(gene_short_name)) %>% group_by(gene_short_name) %>% slice_max(order_by = rowSum, n = 1, with_ties = F)
count_7 = count_7[as.vector(fd_7_x$gene_ID),]
rownames(count_7) = as.vector(fd_7_x$gene_short_name)


### merge seven datasets
gene_overlap = intersect(rownames(count_1), intersect(rownames(count_2), rownames(count_3)))
gene_overlap = intersect(rownames(count_4), gene_overlap)
gene_overlap = intersect(rownames(count_5), gene_overlap)
gene_overlap = intersect(rownames(count_6), gene_overlap)
gene_overlap = intersect(rownames(count_7), gene_overlap)
print(length(gene_overlap))

obj_1 = CreateSeuratObject(count_1[gene_overlap,], meta.data = pd_1)
obj_2 = CreateSeuratObject(count_2[gene_overlap,], meta.data = pd_2)
obj_3 = CreateSeuratObject(count_3[gene_overlap,], meta.data = pd_3)
obj_4 = CreateSeuratObject(count_4[gene_overlap,], meta.data = pd_4)
obj_5 = CreateSeuratObject(count_5[gene_overlap,], meta.data = pd_5)
obj_6 = CreateSeuratObject(count_6[gene_overlap,], meta.data = pd_6)
obj_7 = CreateSeuratObject(count_7[gene_overlap,], meta.data = pd_7)

obj = merge(x = obj_1, y = c(obj_2, obj_3, obj_4, obj_5, obj_6, obj_7))

print(table(obj$group))
print(dim(obj))


obj.list <- SplitObject(obj, split.by = "group")
obj.list <- future_lapply(X = obj.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = obj.list)
obj.list <- future_lapply(X = obj.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = obj.list, reduction = "rpca", 
                                  dims = 1:50)
obj.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)

obj.integrated <- ScaleData(obj.integrated, verbose = FALSE)
obj.integrated <- RunPCA(obj.integrated, npcs = 30, verbose = FALSE)
obj.integrated <- RunUMAP(obj.integrated, dims = 1:30, min.dist = 0.3, n.components = 3)

pd = data.frame(obj.integrated[[]])
pd$UMAP_3d_1 = Embeddings(obj.integrated, reduction = "umap")[,1]
pd$UMAP_3d_2 = Embeddings(obj.integrated, reduction = "umap")[,2]
pd$UMAP_3d_3 = Embeddings(obj.integrated, reduction = "umap")[,3]

saveRDS(obj.integrated, paste0(work_path, "/integration/obj_integrated_seven_datasets.rds"))

pca_coor = Embeddings(obj.integrated, reduction = "pca")
saveRDS(pca_coor, paste0(work_path, "/integration/obj_integrated_seven_datasets_pca.rds"))

obj.integrated <- RunUMAP(obj.integrated, dims = 1:30, min.dist = 0.3, n.components = 2)
pd$UMAP_2d_1 = Embeddings(obj.integrated, reduction = "umap")[,1]
pd$UMAP_2d_2 = Embeddings(obj.integrated, reduction = "umap")[,2]

saveRDS(pd, paste0(work_path, "/integration/obj_integrated_seven_datasets_pd.rds"))


#############################################
### making 2D UMAP on the integrated datasets

group_i = "gastruloid"
pd_x = pd %>% filter(group == group_i) %>% group_by(celltype) %>%
    summarize(UMAP_2d_1_mean = mean(UMAP_2d_1), UMAP_2d_2_mean = mean(UMAP_2d_2)) 
pd_y = data.frame(celltype = names(gastruloid_celltype_color_code)[names(gastruloid_celltype_color_code) %in% pd_x$celltype],
                  celltype_id = 1:sum(names(gastruloid_celltype_color_code) %in% pd_x$celltype))
pd_x = pd_x %>%
    left_join(pd_y, by = "celltype")
try(ggplot(pd) +
        geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.1, color = "grey80") +
        geom_point(data = subset(pd, group == group_i),
                   aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype), size=0.1) +
        ggrepel::geom_text_repel(data = pd_x, 
                                 aes(x = UMAP_2d_1_mean, y = UMAP_2d_2_mean, label = celltype_id), color = "black", size = 4, family = "Arial") +
        theme_void() +
        #        labs(title = group_i) +
        theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=gastruloid_celltype_color_code) +
        ggsave(paste0(work_path, "/plot/Integration_seven_datasets_color_by_celltype_1_nolabel.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)


group_i = "Liberali"
pd_x = pd %>% filter(group == group_i) %>% group_by(celltype) %>%
    summarize(UMAP_2d_1_mean = mean(UMAP_2d_1), UMAP_2d_2_mean = mean(UMAP_2d_2))
pd_y = data.frame(celltype = names(Liberali_celltype_color_code)[names(Liberali_celltype_color_code) %in% pd_x$celltype],
                  celltype_id = 1:sum(names(Liberali_celltype_color_code) %in% pd_x$celltype))
pd_x = pd_x %>%
    left_join(pd_y, by = "celltype")
try(ggplot(pd) +
        geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.1, color = "grey80") +
        geom_point(data = subset(pd, group == group_i),
                   aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype), size=0.1) +
        ggrepel::geom_text_repel(data = pd_x, 
                                 aes(x = UMAP_2d_1_mean, y = UMAP_2d_2_mean, label = celltype_id), color = "black", size = 4, family = "Arial") +
        theme_void() +
        #        labs(title = group_i) +
        theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=Liberali_celltype_color_code) +
        ggsave(paste0(work_path, "/plot/Integration_seven_datasets_color_by_celltype_2_nolabel.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)


group_i = "TLS"
pd_x = pd %>% filter(group == group_i) %>% group_by(celltype) %>%
    summarize(UMAP_2d_1_mean = mean(UMAP_2d_1), UMAP_2d_2_mean = mean(UMAP_2d_2))
pd_y = data.frame(celltype = names(TLS_celltype_color_code)[names(TLS_celltype_color_code) %in% pd_x$celltype],
                  celltype_id = 1:sum(names(TLS_celltype_color_code) %in% pd_x$celltype))
pd_x = pd_x %>%
    left_join(pd_y, by = "celltype")
try(ggplot(pd) +
        geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.1, color = "grey80") +
        geom_point(data = subset(pd, group == group_i),
                   aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype), size=0.1) +
        ggrepel::geom_text_repel(data = pd_x, 
                                 aes(x = UMAP_2d_1_mean, y = UMAP_2d_2_mean, label = celltype_id), color = "black", size = 4, family = "Arial") +
        theme_void() +
        #        labs(title = group_i) +
        theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=TLS_celltype_color_code) +
        ggsave(paste0(work_path, "/plot/Integration_seven_datasets_color_by_celltype_3_nolabel.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)


group_i = "pijuan"
pd_x = pd %>% filter(group == group_i) %>% group_by(celltype) %>%
    summarize(UMAP_2d_1_mean = mean(UMAP_2d_1), UMAP_2d_2_mean = mean(UMAP_2d_2))
pd_y = data.frame(celltype = names(pijuan_celltype_color_code)[names(pijuan_celltype_color_code) %in% pd_x$celltype],
                  celltype_id = 1:sum(names(pijuan_celltype_color_code) %in% pd_x$celltype))
pd_x = pd_x %>%
    left_join(pd_y, by = "celltype")
try(ggplot(pd) +
        geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.1, color = "grey80") +
        geom_point(data = subset(pd, group == group_i),
                   aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype), size=0.1) +
        ggrepel::geom_text_repel(data = pd_x, 
                                 aes(x = UMAP_2d_1_mean, y = UMAP_2d_2_mean, label = celltype_id), color = "black", size = 4, family = "Arial") +
        theme_void() +
        #        labs(title = group_i) +
        theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=pijuan_celltype_color_code) +
        ggsave(paste0(work_path, "/plot/Integration_seven_datasets_color_by_celltype_4_nolabel.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)

group_i = "jax"
pd_x = pd %>% filter(group == group_i) %>% group_by(celltype) %>%
    summarize(UMAP_2d_1_mean = mean(UMAP_2d_1), UMAP_2d_2_mean = mean(UMAP_2d_2))
pd_y = data.frame(celltype = names(jax_celltype_color_code)[names(jax_celltype_color_code) %in% pd_x$celltype],
                  celltype_id = 1:sum(names(jax_celltype_color_code) %in% pd_x$celltype))
pd_x = pd_x %>%
    left_join(pd_y, by = "celltype")
try(ggplot(pd) +
        geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.1, color = "grey80") +
        geom_point(data = subset(pd, group == group_i),
                   aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype), size=0.1) +
        ggrepel::geom_text_repel(data = pd_x, 
                                 aes(x = UMAP_2d_1_mean, y = UMAP_2d_2_mean, label = celltype_id), color = "black", size = 4, family = "Arial") +
        theme_void() +
        #        labs(title = group_i) +
        theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=jax_celltype_color_code) +
        ggsave(paste0(work_path, "/plot/Integration_seven_datasets_color_by_celltype_5_nolabel.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)


group_i = "Rosen"
pd_x = pd %>% filter(group == group_i) %>% group_by(celltype) %>%
    summarize(UMAP_2d_1_mean = mean(UMAP_2d_1), UMAP_2d_2_mean = mean(UMAP_2d_2))
pd_y = data.frame(celltype = names(Rosen_celltype_color_code)[names(Rosen_celltype_color_code) %in% pd_x$celltype],
                  celltype_id = 1:sum(names(Rosen_celltype_color_code) %in% pd_x$celltype))
pd_x = pd_x %>%
    left_join(pd_y, by = "celltype")
try(ggplot(pd) +
        geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.1, color = "grey80") +
        geom_point(data = subset(pd, group == group_i),
                   aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype), size=0.1) +
        ggrepel::geom_text_repel(data = pd_x, 
                                 aes(x = UMAP_2d_1_mean, y = UMAP_2d_2_mean, label = celltype_id), color = "black", size = 4, family = "Arial") +
        theme_void() +
        #        labs(title = group_i) +
        theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=Rosen_celltype_color_code) +
        ggsave(paste0(work_path, "/plot/Integration_seven_datasets_color_by_celltype_6_nolabel.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)


group_i = "Van"
pd_x = pd %>% filter(group == group_i) %>% group_by(celltype) %>%
    summarize(UMAP_2d_1_mean = mean(UMAP_2d_1), UMAP_2d_2_mean = mean(UMAP_2d_2))
pd_y = data.frame(celltype = names(Van_celltype_color_code)[names(Van_celltype_color_code) %in% pd_x$celltype],
                  celltype_id = 1:sum(names(Van_celltype_color_code) %in% pd_x$celltype))
pd_x = pd_x %>%
    left_join(pd_y, by = "celltype")
try(ggplot(pd) +
        geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.1, color = "grey80") +
        geom_point(data = subset(pd, group == group_i),
                   aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype), size=0.1) +
        ggrepel::geom_text_repel(data = pd_x, 
                                 aes(x = UMAP_2d_1_mean, y = UMAP_2d_2_mean, label = celltype_id), color = "black", size = 4, family = "Arial") +
        theme_void() +
        #        labs(title = group_i) +
        theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=Van_celltype_color_code) +
        ggsave(paste0(work_path, "/plot/Integration_seven_datasets_color_by_celltype_7_nolabel.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)


#############################################
### Step-2: Integration subset of brain cells

### 1st dataset - our gastruloid data
dat = readRDS(paste0(work_path, "/obj_processed.rds"))
pd = readRDS(paste0(work_path, "/obj_processed_pd.rds"))
count_1 = GetAssayData(dat, slot = "counts")
pd_1 = data.frame(pd)[,c("plate","celltype")]
pd_1$group = "gastruloid"


keep = pd_1$celltype %in% c("Early neurons","Floor plate","Hindbrain","NMPs","Spinal cord")
count_1 = count_1[,keep]
pd_1 = pd_1[keep,]


### 2nd dataset - natural embryo during organogenesis
pd_all = readRDS(paste0(data_path_3, "/mtx/adata_scale.obs.rds"))

pd_x = readRDS("./Neurons_plus_big_adata_scale.obs.rds")
pd = pd_all[pd_all$cell_id %in% rownames(pd_x),]
rownames(pd) = as.vector(pd$cell_id)

pd = subset(pd, day %in% c("E8.5","E8.75","E9.0","E9.25","E9.5","E9.75","E10.0") & !major_trajectory %in% c("Ependymal_cells", "Intermediate_neuronal_progenitors"))
pd = pd[sample(1:nrow(pd), 200000),]

embryo_list = as.vector(names(table(pd$embryo_id)))

count_4 = NULL
for(j in embryo_list){
    
    print(paste0(j,"/",length(embryo_list)))
    count_j = readRDS(paste0(data_path_3, "/embryo/", j, "_gene_count.rds"))
    keep = colnames(count_j) %in% rownames(pd)
    count_4 = cbind(count_4, count_j[,keep,drop=FALSE])
    rm(count_j)
}
pd_4 = pd[colnames(count_4),]
pd_4$group = "jax"
pd_4$plate = pd_4$sequencing_batch
pd_4$celltype = pd_4$major_trajectory
pd_4$timepoints = pd_4$day
pd_4 = pd_4[,c("plate","celltype","timepoints","group","celltype_update")]

gene_overlap = intersect(rownames(count_1), rownames(count_4))
print(length(gene_overlap))

obj_1 = CreateSeuratObject(count_1[gene_overlap,], meta.data = pd_1)
obj_4 = CreateSeuratObject(count_4[gene_overlap,], meta.data = pd_4)

obj = merge(obj_1, obj_4)

print(table(obj$group))
print(dim(obj))

obj.list <- SplitObject(obj, split.by = "group")
obj.list <- future_lapply(X = obj.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = obj.list)
obj.list <- future_lapply(X = obj.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = obj.list, reduction = "rpca", 
                                  dims = 1:50)
obj.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)

obj.integrated <- ScaleData(obj.integrated, verbose = FALSE)
obj.integrated <- RunPCA(obj.integrated, npcs = 30, verbose = FALSE)
obj.integrated <- RunUMAP(obj.integrated, dims = 1:30, min.dist = 0.3, n.components = 3)

pd = data.frame(obj.integrated[[]])
pd$UMAP_3d_1 = Embeddings(obj.integrated, reduction = "umap")[,1]
pd$UMAP_3d_2 = Embeddings(obj.integrated, reduction = "umap")[,2]
pd$UMAP_3d_3 = Embeddings(obj.integrated, reduction = "umap")[,3]

saveRDS(obj.integrated, paste0(work_path, "/integration/obj_integrated_brain.rds"))

pca_coor = Embeddings(obj.integrated, reduction = "pca")
saveRDS(pca_coor, paste0(work_path, "/integration/obj_integrated_brain_pca.rds"))

obj.integrated <- RunUMAP(obj.integrated, dims = 1:30, min.dist = 0.3, n.components = 2)
pd$UMAP_2d_1 = Embeddings(obj.integrated, reduction = "umap")[,1]
pd$UMAP_2d_2 = Embeddings(obj.integrated, reduction = "umap")[,2]

saveRDS(pd, paste0(work_path, "/integration/obj_integrated_brain_pd.rds"))

########################
### Making 2D UMAP plots

group_i = "gastruloid"
pd_x = pd %>% filter(group == group_i) %>% group_by(celltype) %>%
    summarize(UMAP_2d_1_mean = mean(UMAP_2d_1), UMAP_2d_2_mean = mean(UMAP_2d_2)) 
pd_y = data.frame(celltype = names(gastruloid_celltype_color_code)[names(gastruloid_celltype_color_code) %in% pd_x$celltype],
                  celltype_id = 1:sum(names(gastruloid_celltype_color_code) %in% pd_x$celltype))
pd_x = pd_x %>%
    left_join(pd_y, by = "celltype")
try(ggplot(pd) +
        geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.1, color = "grey80") +
        geom_point(data = subset(pd, group == group_i),
                   aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype), size=0.1) +
        ggrepel::geom_text_repel(data = pd_x, 
                                 aes(x = UMAP_2d_1_mean, y = UMAP_2d_2_mean, label = celltype_id), color = "black", size = 4, family = "Arial") +
        theme_void() +
        #        labs(title = group_i) +
        theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=gastruloid_celltype_color_code) +
        ggsave(paste0(work_path, "/plot/Integration_brain_color_by_celltype_1_nolabel.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)

x = unique(pd$celltype_update[pd$group == "jax"])
names(Rosen_celltype_color_code) = x

group_i = "jax"
pd_x = pd %>% filter(group == group_i) %>% group_by(celltype_update) %>%
    summarize(UMAP_2d_1_mean = mean(UMAP_2d_1), UMAP_2d_2_mean = mean(UMAP_2d_2))
pd_y = data.frame(celltype_update = names(Rosen_celltype_color_code)[names(Rosen_celltype_color_code) %in% pd_x$celltype_update],
                  celltype_id = 1:sum(names(Rosen_celltype_color_code) %in% pd_x$celltype_update))
pd_x = pd_x %>%
    left_join(pd_y, by = "celltype_update")
try(ggplot(pd) +
        geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.1, color = "grey80") +
        geom_point(data = subset(pd, group == group_i),
                   aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype_update), size=0.1) +
        ggrepel::geom_text_repel(data = pd_x, 
                                 aes(x = UMAP_2d_1_mean, y = UMAP_2d_2_mean, label = celltype_id), color = "black", size = 4, family = "Arial") +
        theme_void() +
        #        labs(title = group_i) +
        theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=Rosen_celltype_color_code) +
        ggsave(paste0(work_path, "/plot/Integration_brain_color_by_celltype_2_nolabel.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)

pd_all = readRDS(paste0(data_path_3, "/mtx/adata_scale.obs.rds"))
group_i = "jax"
pd_x = pd %>% mutate(cell_id = rownames(pd)) %>%
    filter(group == group_i) %>% left_join(pd_all[,c("cell_id", "somite_count")], by = "cell_id")
pd_x$somite_count = factor(pd_x$somite_count, levels = somite_list[somite_list %in% pd_x$somite_count])
pd_x = pd_x %>% group_by(somite_count) %>% slice_sample(n = 3000)
try(ggplot(pd) +
        geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.1, color = "grey80") +
        geom_point(data = pd_x,
                   aes(x = UMAP_2d_1, y = UMAP_2d_2, color = somite_count), size=0.1) +
        theme_void() +
        theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=somite_color_plate) +
        ggsave(paste0(work_path, "/plot/Integration_brain_color_by_celltype_3_nolabel.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)



################################################
### Step-3: Integration subset of neuronal cells

### 1st dataset - the gastruloid data
dat = readRDS(paste0(work_path, "/obj_processed.rds"))
pd = readRDS(paste0(work_path, "/obj_processed_pd.rds"))
count_1 = GetAssayData(dat, slot = "counts")
pd_1 = data.frame(pd)[,c("plate","celltype")]
pd_1$group = "gastruloid"

### downsampling to 100K cells
keep = pd_1$celltype %in% c("Early neurons")
count_1 = count_1[,keep]
pd_1 = pd_1[keep,]

### 4th dataset - natural embryo during organogenesis
pd_all = readRDS(paste0(data_path_3, "/mtx/adata_scale.obs.rds"))

pd_x = readRDS("./Neurons_plus_big_adata_scale.obs.rds")
pd = pd_all[pd_all$cell_id %in% rownames(pd_x),]
rownames(pd) = as.vector(pd$cell_id)

pd = subset(pd, day %in% c("E8.5","E8.75","E9.0","E9.25","E9.5","E9.75","E10.0") & major_trajectory %in% c("CNS_neurons"))

embryo_list = as.vector(names(table(pd$embryo_id)))

count_4 = NULL
for(j in embryo_list){
    
    print(paste0(j,"/",length(embryo_list)))
    count_j = readRDS(paste0(data_path_3, "/embryo/", j, "_gene_count.rds"))
    keep = colnames(count_j) %in% rownames(pd)
    count_4 = cbind(count_4, count_j[,keep,drop=FALSE])
    rm(count_j)
}
pd_4 = pd[colnames(count_4),]
pd_4$group = "jax"
pd_4$plate = pd_4$sequencing_batch
pd_4$celltype = pd_4$major_trajectory
pd_4$timepoints = pd_4$day
pd_4 = pd_4[,c("plate","celltype","timepoints","group","celltype_update")]

gene_overlap = intersect(rownames(count_1), rownames(count_4))
print(length(gene_overlap))

obj_1 = CreateSeuratObject(count_1[gene_overlap,], meta.data = pd_1)
obj_4 = CreateSeuratObject(count_4[gene_overlap,], meta.data = pd_4)

obj = merge(obj_1, obj_4)

print(table(obj$group))
print(dim(obj))

obj.list <- SplitObject(obj, split.by = "group")
obj.list <- future_lapply(X = obj.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = obj.list)
obj.list <- future_lapply(X = obj.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = obj.list, reduction = "rpca", 
                                  dims = 1:50)
obj.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)

obj.integrated <- ScaleData(obj.integrated, verbose = FALSE)
obj.integrated <- RunPCA(obj.integrated, npcs = 30, verbose = FALSE)
obj.integrated <- RunUMAP(obj.integrated, dims = 1:30, min.dist = 0.3, n.components = 3)

pd = data.frame(obj.integrated[[]])
pd$UMAP_3d_1 = Embeddings(obj.integrated, reduction = "umap")[,1]
pd$UMAP_3d_2 = Embeddings(obj.integrated, reduction = "umap")[,2]
pd$UMAP_3d_3 = Embeddings(obj.integrated, reduction = "umap")[,3]

saveRDS(obj.integrated, paste0(work_path, "/integration/obj_integrated_neurons.rds"))

pca_coor = Embeddings(obj.integrated, reduction = "pca")
saveRDS(pca_coor, paste0(work_path, "/integration/obj_integrated_neurons_pca.rds"))

obj.integrated <- RunUMAP(obj.integrated, dims = 1:30, min.dist = 0.3, n.components = 2)
pd$UMAP_2d_1 = Embeddings(obj.integrated, reduction = "umap")[,1]
pd$UMAP_2d_2 = Embeddings(obj.integrated, reduction = "umap")[,2]

saveRDS(pd, paste0(work_path, "/integration/obj_integrated_neurons_pd.rds"))

########################
### Making 2D UMAP plots

group_i = "gastruloid"
pd_x = pd %>% filter(group == group_i) %>% group_by(celltype) %>%
    summarize(UMAP_2d_1_mean = mean(UMAP_2d_1), UMAP_2d_2_mean = mean(UMAP_2d_2)) 
pd_y = data.frame(celltype = names(gastruloid_celltype_color_code)[names(gastruloid_celltype_color_code) %in% pd_x$celltype],
                  celltype_id = 1:sum(names(gastruloid_celltype_color_code) %in% pd_x$celltype))
pd_x = pd_x %>%
    left_join(pd_y, by = "celltype")
try(ggplot(pd) +
        geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.1, color = "grey80") +
        geom_point(data = subset(pd, group == group_i),
                   aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype), size=0.1) +
        #ggrepel::geom_text_repel(data = pd_x, 
        #                         aes(x = UMAP_2d_1_mean, y = UMAP_2d_2_mean, label = celltype_id), color = "black", size = 4, family = "Arial") +
        theme_void() +
        #        labs(title = group_i) +
        theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=gastruloid_celltype_color_code) +
        ggsave(paste0(work_path, "/plot/Integration_neurons_color_by_celltype_1_nolabel.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)

pd_all = readRDS(paste0(data_path_3, "/mtx/adata_scale.obs.rds"))
pd_z = pd_all[rownames(pd_all) %in% rownames(pd),]
pd_z_1 = pd_z[!is.na(pd_z$neurons_sub_clustering),]
pd_z_2 = pd_z[is.na(pd_z$neurons_sub_clustering),]; pd_z_2$neurons_sub_clustering = as.vector(pd_z_2$celltype_update)
pd_z = rbind(pd_z_1, pd_z_2)
x = table(pd_z$neurons_sub_clustering)

tmp_celltype_color_code = Rosen_celltype_color_code[1:length(x)]
names(tmp_celltype_color_code) = names(x)

group_i = "jax"
pd_x = pd %>% mutate(cell_id = rownames(pd)) %>%
    filter(group == group_i) %>% 
    left_join(pd_z[,c("cell_id","neurons_sub_clustering")], by = "cell_id") %>% 
    group_by(neurons_sub_clustering) %>%
    summarize(UMAP_2d_1_mean = mean(UMAP_2d_1), UMAP_2d_2_mean = mean(UMAP_2d_2))
pd_y = data.frame(neurons_sub_clustering = names(tmp_celltype_color_code)[names(tmp_celltype_color_code) %in% pd_x$neurons_sub_clustering],
                  celltype_id = 1:sum(names(tmp_celltype_color_code) %in% pd_x$neurons_sub_clustering))
pd_x = pd_x %>%
    left_join(pd_y, by = "neurons_sub_clustering")
try(ggplot(pd) +
        geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.1, color = "grey80") +
        geom_point(data = pd %>% mutate(cell_id = rownames(pd)) %>%
                       filter(group == group_i) %>% 
                       left_join(pd_z[,c("cell_id","neurons_sub_clustering")], by = "cell_id"),
                   aes(x = UMAP_2d_1, y = UMAP_2d_2, color = neurons_sub_clustering), size=0.1) +
        ggrepel::geom_text_repel(data = pd_x, 
                                 aes(x = UMAP_2d_1_mean, y = UMAP_2d_2_mean, label = celltype_id), color = "black", size = 4, family = "Arial") +
        theme_void() +
        #        labs(title = group_i) +
        theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=tmp_celltype_color_code) +
        ggsave(paste0(work_path, "/plot/Integration_neurons_color_by_celltype_2_nolabel.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)

pd_all = readRDS(paste0(data_path_3, "/mtx/adata_scale.obs.rds"))

group_i = "jax"
pd_x = pd %>% mutate(cell_id = rownames(pd)) %>%
    filter(group == group_i) %>% left_join(pd_all[,c("cell_id", "somite_count")], by = "cell_id")
pd_x$somite_count = factor(pd_x$somite_count, levels = somite_list[somite_list %in% pd_x$somite_count])
pd_x = pd_x %>% group_by(somite_count) %>% slice_sample(n = 2000)
try(ggplot(pd) +
        geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.1, color = "grey80") +
        geom_point(data = pd_x,
                   aes(x = UMAP_2d_1, y = UMAP_2d_2, color = somite_count), size=0.1) +
        theme_void() +
        theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=somite_color_plate) +
        ggsave(paste0(work_path, "/plot/Integration_neurons_color_by_celltype_3_nolabel.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)


###########################################################
### Step-4: Integration subset of data from specific stages

### 1st dataset - Sam's gastruloid data
dat = readRDS(paste0(work_path, "/obj_processed.rds"))
pd = readRDS(paste0(work_path, "/obj_processed_pd.rds"))
count_1 = GetAssayData(dat, slot = "counts")
pd_1 = data.frame(pd)[,c("plate","celltype")]
pd_1$group = "gastruloid"

### downsampling to 100K cells
keep = rownames(pd_1) %in% rownames(pd_five)
count_1 = count_1[,keep]
pd_1 = pd_1[keep,]

fd_1 = mouse_gene[rownames(count_1),]
fd_1$rowSum = Matrix::rowSums(count_1)
fd_1_x = fd_1 %>% group_by(gene_short_name) %>% slice_max(order_by = rowSum, n = 1, with_ties = F)
count_1 = count_1[as.vector(fd_1_x$gene_ID),]
rownames(count_1) = as.vector(fd_1_x$gene_short_name)


### 2nd dataset - Liberali's data
count_2 = Matrix::readMM(paste0(data_path, "/gastruloid_Liberali/GSE229513_UMI_counts.mtx"))
rownames(count_2) = as.vector(read.table(paste0(data_path, "/gastruloid_Liberali/GSE229513_genes.tsv"), header=T)$gene)
colnames(count_2) = as.vector(read.table(paste0(data_path, "/gastruloid_Liberali/GSE229513_barcodes.tsv"), sep="\t")$V1)
pd_2 = readRDS(paste0(data_path, "/gastruloid_Liberali/df_cell.rds"))
pd_2 = pd_2[,c(1:4, 9, 13)]
pd_2$celltype = pd_2$celltypeannotation
pd_2$plate = pd_2$batch
pd_2 = pd_2[,c("plate","timepoints","celltype")]
pd_2$group = "Liberali"

### downsampling to 100K cells
keep = pd_2$timepoints == "120h"
count_2 = count_2[,keep]
pd_2 = pd_2[keep,]


### 3rd dataset - TLS data
count_TLS_96h = Read10X(paste0(data_path, "/gastruloid_TLS/TLS_96h"), gene.column = 1, strip.suffix = T)
count_TLS_108h = Read10X(paste0(data_path, "/gastruloid_TLS/TLS_108h"), gene.column = 1, strip.suffix = T)
count_TLS_120h = Read10X(paste0(data_path, "/gastruloid_TLS/TLS_120h"), gene.column = 1, strip.suffix = T)

meta_data = read.table(paste0(data_path, "/gastruloid_TLS/TLS_meta_data.tsv"), as.is=T, header=T)
meta_data$sample = paste0(meta_data$TP, "_", meta_data$BC)
rownames(meta_data) = as.vector(meta_data$sample)

colnames(count_TLS_96h) = paste0("TLS_96h_", colnames(count_TLS_96h))
count_TLS_96h = count_TLS_96h[,colnames(count_TLS_96h) %in% meta_data$sample]
colnames(count_TLS_108h) = paste0("TLS_108h_", colnames(count_TLS_108h))
count_TLS_108h = count_TLS_108h[,colnames(count_TLS_108h) %in% meta_data$sample]
colnames(count_TLS_120h) = paste0("TLS_120h_", colnames(count_TLS_120h))
count_TLS_120h = count_TLS_120h[,colnames(count_TLS_120h) %in% meta_data$sample]

count_3 = cbind(count_TLS_96h, count_TLS_108h, count_TLS_120h)
pd_3 = meta_data[colnames(count_3),]
pd_3$celltype = pd_3$cell_state
pd_3$plate = pd_3$TP
pd_3$group = "TLS"
pd_3 = pd_3[,c("plate","celltype","group")]

fd_3 = mouse_gene[rownames(count_3),]
fd_3$rowSum = Matrix::rowSums(count_3)
fd_3_x = fd_3 %>% filter(!is.na(gene_short_name)) %>% group_by(gene_short_name) %>% slice_max(order_by = rowSum, n = 1, with_ties = F)
count_3 = count_3[as.vector(fd_3_x$gene_ID),]
rownames(count_3) = as.vector(fd_3_x$gene_short_name)


### downsampling to 100K cells
keep = pd_3$plate == "TLS_120h"
count_3 = count_3[,keep]
pd_3 = pd_3[keep,]


### sixth dataset - Rosen's gastruloid
obj_6 = readRDS(paste0(data_path, "/gastruloid_Rosen/GSE212050_seurat_final.rds"))
fd_6 = read.table(paste0(data_path, "/gastruloid_Rosen/GSE212050_feature_metadata_final.txt"), header=T, as.is=T, sep=",")[,c(1,2)]
names(fd_6) = c("gene_ID","gene_short_name")
pd_6 = data.frame(obj_6[[]])[,c("experiment", "cluster","timepoint")]
pd_6$group = "Rosen"
names(pd_6) = c("plate","celltype","timepoints","group")

count_6 = GetAssayData(obj_6, slot = "counts")
rm(obj_6)

sum(fd_6$gene_ID == rownames(count_6))
fd_6$rowSum = Matrix::rowSums(count_6)
fd_6_x = fd_6 %>% filter(!is.na(gene_short_name)) %>% group_by(gene_short_name) %>% slice_max(order_by = rowSum, n = 1, with_ties = F)
count_6 = count_6[as.vector(fd_6_x$gene_ID),]
rownames(count_6) = as.vector(fd_6_x$gene_short_name)


### downsampling to 100K cells
keep = pd_6$timepoints == "d5"
count_6 = count_6[,keep]
pd_6 = pd_6[keep,]


### invitroth dataset - Van's gastruloid
dat = readRDS(paste0(data_path, "/gastruloid_van/dat.rds"))
count_7 = dat[['count']]
fd_7 = dat[['fd']]
pd_7 = dat[['pd']]
names(fd_7) = c("gene_ID","gene_short_name")
pd_7 = pd_7[,c("batch", "celltye")]
pd_7$group = "Van"
names(pd_7) = c("plate","celltype","group")

sum(fd_7$gene_ID == rownames(count_7))
fd_7$rowSum = Matrix::rowSums(count_7)
fd_7_x = fd_7 %>% filter(!is.na(gene_short_name)) %>% group_by(gene_short_name) %>% slice_max(order_by = rowSum, n = 1, with_ties = F)
count_7 = count_7[as.vector(fd_7_x$gene_ID),]
rownames(count_7) = as.vector(fd_7_x$gene_short_name)


### merge invitro datasets
gene_overlap = intersect(rownames(count_1), intersect(rownames(count_2), rownames(count_3)))
gene_overlap = intersect(rownames(count_6), gene_overlap)
gene_overlap = intersect(rownames(count_7), gene_overlap)
print(length(gene_overlap))

obj_1 = CreateSeuratObject(count_1[gene_overlap,], meta.data = pd_1)
obj_2 = CreateSeuratObject(count_2[gene_overlap,], meta.data = pd_2)
obj_3 = CreateSeuratObject(count_3[gene_overlap,], meta.data = pd_3)
obj_6 = CreateSeuratObject(count_6[gene_overlap,], meta.data = pd_6)
obj_7 = CreateSeuratObject(count_7[gene_overlap,], meta.data = pd_7)



obj = merge(x = obj_1, y = c(obj_2, obj_3, obj_6, obj_7))

print(table(obj$group))
print(dim(obj))


obj.list <- SplitObject(obj, split.by = "group")
obj.list <- future_lapply(X = obj.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = obj.list)
obj.list <- future_lapply(X = obj.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = obj.list, reduction = "rpca", 
                                  dims = 1:50)
obj.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)

obj.integrated <- ScaleData(obj.integrated, verbose = FALSE)
obj.integrated <- RunPCA(obj.integrated, npcs = 30, verbose = FALSE)
obj.integrated <- RunUMAP(obj.integrated, dims = 1:30, min.dist = 0.3, n.components = 3)

pd = data.frame(obj.integrated[[]])
pd$UMAP_3d_1 = Embeddings(obj.integrated, reduction = "umap")[,1]
pd$UMAP_3d_2 = Embeddings(obj.integrated, reduction = "umap")[,2]
pd$UMAP_3d_3 = Embeddings(obj.integrated, reduction = "umap")[,3]

saveRDS(obj.integrated, paste0(work_path, "/integration/obj_integrated_invitro_datasets_120h.rds"))

pca_coor = Embeddings(obj.integrated, reduction = "pca")
saveRDS(pca_coor, paste0(work_path, "/integration/obj_integrated_invitro_datasets_120h_pca.rds"))

obj.integrated <- RunUMAP(obj.integrated, dims = 1:30, min.dist = 0.3, n.components = 2)
pd$UMAP_2d_1 = Embeddings(obj.integrated, reduction = "umap")[,1]
pd$UMAP_2d_2 = Embeddings(obj.integrated, reduction = "umap")[,2]

saveRDS(pd, paste0(work_path, "/integration/obj_integrated_invitro_datasets_120h_pd.rds"))




#############################################
### making 2D UMAP on the integrated datasets


group_i = "gastruloid"
pd_x = pd %>% filter(group == group_i) %>% group_by(celltype) %>%
    summarize(UMAP_2d_1_mean = mean(UMAP_2d_1), UMAP_2d_2_mean = mean(UMAP_2d_2)) 
pd_y = data.frame(celltype = names(gastruloid_celltype_color_code)[names(gastruloid_celltype_color_code) %in% pd_x$celltype],
                  celltype_id = 1:sum(names(gastruloid_celltype_color_code) %in% pd_x$celltype))
pd_x = pd_x %>%
    left_join(pd_y, by = "celltype")
try(ggplot(pd) +
        geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.1, color = "grey80") +
        geom_point(data = subset(pd, group == group_i),
                   aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype), size=0.1) +
        ggrepel::geom_text_repel(data = pd_x, 
                                 aes(x = UMAP_2d_1_mean, y = UMAP_2d_2_mean, label = celltype_id), color = "black", size = 4, family = "Arial") +
        theme_void() +
        #        labs(title = group_i) +
        theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=gastruloid_celltype_color_code) +
        ggsave(paste0(work_path, "/plot/Integration_invitro_datasets_120h_color_by_celltype_1_nolabel.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)


group_i = "Liberali"
pd_x = pd %>% filter(group == group_i) %>% group_by(celltype) %>%
    summarize(UMAP_2d_1_mean = mean(UMAP_2d_1), UMAP_2d_2_mean = mean(UMAP_2d_2))
pd_y = data.frame(celltype = names(Liberali_celltype_color_code),
                  celltype_id = 1:length(Liberali_celltype_color_code))
pd_x = pd_x %>%
    left_join(pd_y, by = "celltype")
try(ggplot(pd) +
        geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.1, color = "grey80") +
        geom_point(data = subset(pd, group == group_i),
                   aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype), size=0.1) +
        ggrepel::geom_text_repel(data = pd_x, 
                                 aes(x = UMAP_2d_1_mean, y = UMAP_2d_2_mean, label = celltype_id), color = "black", size = 4, family = "Arial") +
        theme_void() +
        #        labs(title = group_i) +
        theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=Liberali_celltype_color_code) +
        ggsave(paste0(work_path, "/plot/Integration_invitro_datasets_120h_color_by_celltype_2_nolabel.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)


group_i = "TLS"
pd_x = pd %>% filter(group == group_i) %>% group_by(celltype) %>%
    summarize(UMAP_2d_1_mean = mean(UMAP_2d_1), UMAP_2d_2_mean = mean(UMAP_2d_2))
pd_y = data.frame(celltype = names(TLS_celltype_color_code)[names(TLS_celltype_color_code) %in% pd_x$celltype],
                  celltype_id = 1:sum(names(TLS_celltype_color_code) %in% pd_x$celltype))
pd_x = pd_x %>%
    left_join(pd_y, by = "celltype")
try(ggplot(pd) +
        geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.1, color = "grey80") +
        geom_point(data = subset(pd, group == group_i),
                   aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype), size=0.1) +
        ggrepel::geom_text_repel(data = pd_x, 
                                 aes(x = UMAP_2d_1_mean, y = UMAP_2d_2_mean, label = celltype_id), color = "black", size = 4, family = "Arial") +
        theme_void() +
        #        labs(title = group_i) +
        theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=TLS_celltype_color_code) +
        ggsave(paste0(work_path, "/plot/Integration_invitro_datasets_120h_color_by_celltype_3_nolabel.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)

group_i = "Rosen"
pd_x = pd %>% filter(group == group_i) %>% group_by(celltype) %>%
    summarize(UMAP_2d_1_mean = mean(UMAP_2d_1), UMAP_2d_2_mean = mean(UMAP_2d_2))
pd_y = data.frame(celltype = names(Rosen_celltype_color_code),
                  celltype_id = 1:length(Rosen_celltype_color_code))
pd_x = pd_x %>%
    left_join(pd_y, by = "celltype")
try(ggplot(pd) +
        geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.1, color = "grey80") +
        geom_point(data = subset(pd, group == group_i),
                   aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype), size=0.1) +
        ggrepel::geom_text_repel(data = pd_x, 
                                 aes(x = UMAP_2d_1_mean, y = UMAP_2d_2_mean, label = celltype_id), color = "black", size = 4, family = "Arial") +
        theme_void() +
        #        labs(title = group_i) +
        theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=Rosen_celltype_color_code) +
        ggsave(paste0(work_path, "/plot/Integration_invitro_datasets_120h_color_by_celltype_4_nolabel.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)


group_i = "Van"
pd_x = pd %>% filter(group == group_i) %>% group_by(celltype) %>%
    summarize(UMAP_2d_1_mean = mean(UMAP_2d_1), UMAP_2d_2_mean = mean(UMAP_2d_2))
pd_y = data.frame(celltype = names(Van_celltype_color_code)[names(Van_celltype_color_code) %in% pd_x$celltype],
                  celltype_id = 1:sum(names(Van_celltype_color_code) %in% pd_x$celltype))
pd_x = pd_x %>%
    left_join(pd_y, by = "celltype")
try(ggplot(pd) +
        geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.1, color = "grey80") +
        geom_point(data = subset(pd, group == group_i),
                   aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype), size=0.1) +
        ggrepel::geom_text_repel(data = pd_x, 
                                 aes(x = UMAP_2d_1_mean, y = UMAP_2d_2_mean, label = celltype_id), color = "black", size = 4, family = "Arial") +
        theme_void() +
        #        labs(title = group_i) +
        theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=Van_celltype_color_code) +
        ggsave(paste0(work_path, "/plot/Integration_invitro_datasets_120h_color_by_celltype_5_nolabel.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)






