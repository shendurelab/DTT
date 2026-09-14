
########################################################################
### Integrating the new dataset with the old dataset, to help annotation

source("~/work/scripts/utils.R")
mouse_gene = read.table("~/work/tome/code/refdata-gex-GRCm39-2024-A.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape"
data_path = "/net/shendure/vol10/projects/prime_editing/nobackup/data/2025_11_21_revision"

dataset = "monoclonal_polyclonal_RNAseq"
pd = readRDS(paste0(work_path, "/revision/", dataset, "/obj_batchelor_pd.rds"))

dataset = "monoclonal_RNAseq"
count_1 = readRDS(paste0(work_path, "/revision/", dataset, "/count.rds"))

dataset = "polyclonal_RNAseq"
count_2 = readRDS(paste0(work_path, "/revision/", dataset, "/count.rds"))

count = cbind(count_1, count_2)
mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "M")) &
              mouse_gene$gene_type %in% c("protein_coding", "lncRNA"),]
count = count[rownames(count) %in% as.vector(mouse_gene_sub$gene_ID),]

pd$celltype = paste0("cluster_", as.vector(pd$RNA_snn_res.1))
pd = pd[,c("cell_id", "experiment_id", "log2_umi", "MT_pct", "RIBO_pct", "doublet_score", "dataset", "celltype")]
sum(pd$cell_id == colnames(count))

dataset = "mGASv5"
count_3 = readRDS(paste0(work_path, "/data_mGASv5/count.rds"))
pd_3 = readRDS(paste0(work_path, "/data_mGASv5/obj_sctransform_pd.rds"))
pd_3$dataset = dataset
pd_3 = pd_3[,c("cell_id", "experiment_id", "log2_umi", "MT_pct", "RIBO_pct", "doublet_score", "dataset", "celltype")]
sum(pd_3$cell_id == colnames(count_3))

gene_overlap = intersect(rownames(count), rownames(count_3))
count = cbind(count[gene_overlap,], count_3[gene_overlap,])
pd = rbind(pd, pd_3)

mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "M")) &
              mouse_gene$gene_type %in% c("protein_coding", "lncRNA"),]
count = count[rownames(count) %in% as.vector(mouse_gene_sub$gene_ID),]

obj = CreateSeuratObject(count, meta.data = pd)
dataset = "monoclonal_polyclonal_RNAseq"

obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2500)
obj = ScaleData(object = obj, verbose = FALSE)
obj = RunPCA(object = obj, npcs = 30, verbose = FALSE)

pca_coor = as.matrix(Embeddings(obj, reduction = "pca"))
aligned_coor = batchelor::reducedMNN(as.matrix(pca_coor),
                                     batch = data.frame(obj[[]])[,"dataset"],
                                     k=20)
aligned_coor = aligned_coor$corrected

colnames(aligned_coor) = paste0("aligned_", 1:30)
rownames(aligned_coor) = colnames(obj)

obj[['aligned']] = Seurat::CreateDimReducObject(embeddings=as.matrix(aligned_coor), key='aligned_')

obj = RunUMAP(object = obj, reduction = "aligned", dims = 1:30, min.dist = 0.3, verbose = FALSE)
obj = FindNeighbors(object = obj, reduction = "aligned", dims = 1:30, verbose = FALSE)
obj = FindClusters(object = obj, resolution = 1, verbose = FALSE)
obj_processed = obj
saveRDS(obj_processed, paste0(work_path, "/revision/", dataset, "/obj_batchelor_old_dataset.rds"))

pd = data.frame(obj_processed[[]])
pd$UMAP_1 = as.vector(Embeddings(obj_processed, reduction = "umap")[,1])
pd$UMAP_2 = as.vector(Embeddings(obj_processed, reduction = "umap")[,2])
pd$group = NULL
saveRDS(pd, paste0(work_path, "/revision/", dataset, "/obj_batchelor_old_dataset_pd.rds"))

cds = new_cell_data_set(count[rownames(obj_processed), colnames(obj_processed)],
                        cell_metadata = data.frame(obj_processed[[]]),
                        gene_metadata = mouse_gene[rownames(obj_processed),])
cds = preprocess_cds(cds, use_genes = rownames(obj)[1:1000])
reducedDims(cds)$UMAP = as.matrix(Embeddings(obj_processed, reduction = "umap"))
saveRDS(cds, paste0(work_path, "/revision/", dataset, "/cds_batchelor_old_dataset.rds"))




##########################################################
### step-2: annotating cell types based on the integration using monocle/3

source("~/work/scripts/utils.R")
mouse_gene = read.table("~/work/tome/code/refdata-gex-GRCm39-2024-A.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape"
dataset = "monoclonal_polyclonal_RNAseq"

pd = readRDS(paste0(work_path, "/revision/", dataset, "/obj_batchelor_pd.rds"))

anno = rep(NA, nrow(pd))
anno[pd$RNA_snn_res.1 %in% c(1,13)] = "Epiblast"
anno[pd$RNA_snn_res.1 %in% c(22)] = "ESCs (2-cell state)"
anno[pd$RNA_snn_res.1 %in% c(14,15)] = "Transitional cells"

anno[pd$RNA_snn_res.1 %in% c(0,4,9,10,12,16)] = "Somites"
anno[pd$RNA_snn_res.1 %in% c(20)] = "Endothelial cells"

anno[pd$RNA_snn_res.1 %in% c(6,7)] = "NMPs"
anno[pd$RNA_snn_res.1 %in% c(3,5)] = "Spinal cord"
anno[pd$RNA_snn_res.1 %in% c(2,11)] = "Hindbrain"
anno[pd$RNA_snn_res.1 %in% c(17)] = "Early neurons"

anno[pd$RNA_snn_res.1 %in% c(8)] = "Node-like cells"
anno[pd$RNA_snn_res.1 %in% c(18)] = "Notochord"
anno[pd$RNA_snn_res.1 %in% c(21)] = "Definitive endoderm"

anno[pd$RNA_snn_res.1 %in% c(19)] = "Extraembryonic endoderm"

pd$celltype = as.vector(anno)
saveRDS(pd, paste0(work_path, "/revision/", dataset, "/obj_batchelor_pd.rds"))

pd_x = pd %>% group_by(celltype) %>% summarize(UMAP_1_mean = mean(UMAP_1), UMAP_2_mean = mean(UMAP_2))

p = ggplot() +
        geom_point(data = pd, aes(x = UMAP_1, y = UMAP_2), size=0.15, color = "black") +
        geom_point(data = pd,
                   aes(x = UMAP_1, y = UMAP_2, color = celltype), size=0.1) +
        ggrepel::geom_text_repel(data = pd_x, aes(x = UMAP_1_mean, y = UMAP_2_mean, label = celltype), color = "black", size = 3, family = "Arial") +
        theme_void() +
        theme(legend.position="none") +
        labs(title = "New experiment") +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=gastruloid_celltype_color_code)
ggsave(paste0(work_path, "/revision/", dataset, "/UMAP_celltype_title.png"), p,
               dpi = 300,
               height  = 6, 
               width = 8)

### One cluster (RNA_snn_res.1 = 20) of Transitional cells is a little messy, so I put it in a lower layer when plotting
try(ggplot() +
        geom_point(data = pd, aes(x = UMAP_1, y = UMAP_2), size=0.15, color = "black") +
        geom_point(data = pd[pd$RNA_snn_res.1 == 20,],
                   aes(x = UMAP_1, y = UMAP_2, color = celltype), size=0.1) +
        geom_point(data = pd[pd$RNA_snn_res.1 != 20,],
                   aes(x = UMAP_1, y = UMAP_2, color = celltype), size=0.1) +
        theme_void() +
        theme(legend.position="none") +
        # labs(title = "mGASv5") +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=gastruloid_celltype_color_code) +
        ggsave(paste0(work_path, "/revision/", dataset, "/UMAP_celltype.png"),
               dpi = 300,
               height  = 6, 
               width = 7), silent = TRUE)

try(ggplot() +
        geom_point(data = pd, aes(x = UMAP_1, y = UMAP_2), size=0.1, color = "grey80") +
        geom_point(data = pd[pd$dataset == "monoclonal_RNAseq",], aes(x = UMAP_1, y = UMAP_2), size=0.15, color = "red") +
        theme_void() +
        theme(legend.position="none") +
        # labs(title = "mGASv5") +
        theme(plot.title = element_text(hjust = 0.5)) +
        ggsave(paste0(work_path, "/revision/", dataset, "/UMAP_celltype_monoclonal.png"),
               dpi = 300,
               height  = 6, 
               width = 7), silent = TRUE)

try(ggplot() +
        geom_point(data = pd, aes(x = UMAP_1, y = UMAP_2), size=0.1, color = "grey80") +
        geom_point(data = pd[pd$dataset == "polyclonal_RNAseq",], aes(x = UMAP_1, y = UMAP_2), size=0.15, color = "blue") +
        theme_void() +
        theme(legend.position="none") +
        # labs(title = "mGASv5") +
        theme(plot.title = element_text(hjust = 0.5)) +
        ggsave(paste0(work_path, "/revision/", dataset, "/UMAP_celltype_polyclonal.png"),
               dpi = 300,
               height  = 6, 
               width = 7), silent = TRUE)


> summary(pd$UMI_count)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
   1024    2833    4045    4820    5639   95517

df = pd %>% group_by(dataset, celltype) %>% tally() %>%
left_join(pd %>% group_by(dataset) %>% tally() %>% rename(total_n = n), by = "dataset") %>%
mutate(frac = 100*n/total_n)









