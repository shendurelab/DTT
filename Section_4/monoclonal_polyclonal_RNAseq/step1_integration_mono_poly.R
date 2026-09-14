
#############################################################
### Merging monoclonal gastruloids and polyclonal gastruloids

##################################################################################
### Test-1: directly combining the two datasets without performing any integration

source("~/work/scripts/utils.R")
mouse_gene = read.table("~/work/tome/code/refdata-gex-GRCm39-2024-A.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape"
data_path = "/net/shendure/vol10/projects/prime_editing/nobackup/data/2025_11_21_revision"

dataset = "monoclonal_RNAseq"
count_1 = readRDS(paste0(work_path, "/revision/", dataset, "/count.rds"))
pd_1 = readRDS(paste0(work_path, "/revision/", dataset, "/pd.rds"))
pd_1$dataset = dataset

dataset = "polyclonal_RNAseq"
count_2 = readRDS(paste0(work_path, "/revision/", dataset, "/count.rds"))
pd_2 = readRDS(paste0(work_path, "/revision/", dataset, "/pd.rds"))
pd_2$dataset = dataset

count = cbind(count_1, count_2)
pd = rbind(pd_1, pd_2)

mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "M")) &
              mouse_gene$gene_type %in% c("protein_coding", "lncRNA"),]
count = count[rownames(count) %in% as.vector(mouse_gene_sub$gene_ID),]
obj = CreateSeuratObject(count, meta.data = pd)

dataset = "monoclonal_polyclonal_RNAseq"

obj$group = dataset
obj_processed = doClusterSeurat(obj)
saveRDS(obj_processed, paste0(work_path, "/revision/", dataset, "/obj_processed.rds"))

pd = data.frame(obj_processed[[]])
pd$UMAP_1 = as.vector(Embeddings(obj_processed, reduction = "umap")[,1])
pd$UMAP_2 = as.vector(Embeddings(obj_processed, reduction = "umap")[,2])
pd$group = NULL
saveRDS(pd, paste0(work_path, "/revision/", dataset, "/obj_processed_pd.rds"))

cds = new_cell_data_set(GetAssayData(obj_processed, slot = "counts"),
                        cell_metadata = data.frame(obj_processed[[]]),
                        gene_metadata = mouse_gene[rownames(obj_processed),])
cds = preprocess_cds(cds, use_genes = VariableFeatures(obj_processed))
reducedDims(cds)$UMAP = as.matrix(Embeddings(obj_processed, reduction = "umap"))
saveRDS(cds, paste0(work_path, "/revision/", dataset, "/cds_processed.rds"))


pd = readRDS(paste0(work_path, "/revision/", dataset, "/obj_processed_pd.rds"))
p1 = ggplot() +
     geom_point(data = pd, aes(x = UMAP_1, y = UMAP_2), color = "grey80", size = 0.1) +
     geom_point(data = pd[pd$dataset == "monoclonal_RNAseq",], aes(x = UMAP_1, y = UMAP_2), color = "red", size = 0.15) +
     theme_void() +
     theme(legend.position="none") +
     labs(title = "monoclonal_RNAseq") +
     theme(plot.title = element_text(hjust = 0.5))
p2 = ggplot() +
     geom_point(data = pd, aes(x = UMAP_1, y = UMAP_2), color = "grey80", size = 0.1) +
     geom_point(data = pd[pd$dataset == "polyclonal_RNAseq",], aes(x = UMAP_1, y = UMAP_2), color = "blue", size = 0.15) +
     theme_void() +
     theme(legend.position="none") +
     labs(title = "polyclonal_RNAseq") +
     theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("~/share/", "simple_combine", ".png"), p1 / p2, width = 7, height = 10, dpi = 300)



##############################
### Test-2: trying SCTransform

source("~/work/scripts/utils.R")
mouse_gene = read.table("~/work/tome/code/refdata-gex-GRCm39-2024-A.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

options(future.globals.maxSize = 100 * 1024^3)

work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape"
data_path = "/net/shendure/vol10/projects/prime_editing/nobackup/data/2025_11_21_revision"

dataset = "monoclonal_RNAseq"
count_1 = readRDS(paste0(work_path, "/revision/", dataset, "/count.rds"))
pd_1 = readRDS(paste0(work_path, "/revision/", dataset, "/pd.rds"))
pd_1$dataset = dataset

dataset = "polyclonal_RNAseq"
count_2 = readRDS(paste0(work_path, "/revision/", dataset, "/count.rds"))
pd_2 = readRDS(paste0(work_path, "/revision/", dataset, "/pd.rds"))
pd_2$dataset = dataset

count = cbind(count_1, count_2)
pd = rbind(pd_1, pd_2)

mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "M")) &
              mouse_gene$gene_type %in% c("protein_coding", "lncRNA"),]
count = count[rownames(count) %in% as.vector(mouse_gene_sub$gene_ID),]
obj = CreateSeuratObject(count, meta.data = pd)

dataset = "monoclonal_polyclonal_RNAseq"
obj$group = dataset

obj = SCTransform(object = obj, verbose = FALSE)
obj = RunPCA(object = obj, npcs = 30, verbose = FALSE)
obj = RunUMAP(object = obj, dims = 1:30, min.dist = 0.3, verbose = FALSE)
obj = FindNeighbors(object = obj, dims = 1:30, verbose = FALSE, reduction = "pca")
obj = FindClusters(object = obj, resolution = 1, verbose = FALSE)
saveRDS(obj, paste0(work_path, "/revision/", dataset, "/obj_sctransform.rds"))

pd = data.frame(obj[[]])
pd$UMAP_1 = as.vector(Embeddings(obj, reduction = "umap")[,1])
pd$UMAP_2 = as.vector(Embeddings(obj, reduction = "umap")[,2])
pd$group = NULL
saveRDS(pd, paste0(work_path, "/revision/", dataset, "/obj_sctransform_pd.rds"))

cds = new_cell_data_set(count[rownames(obj), colnames(obj)],
                        cell_metadata = data.frame(obj[[]]),
                        gene_metadata = mouse_gene[rownames(obj),])
cds = preprocess_cds(cds, use_genes = rownames(obj)[1:1000])
reducedDims(cds)$UMAP = as.matrix(Embeddings(obj, reduction = "umap"))
saveRDS(cds, paste0(work_path, "/revision/", dataset, "/cds_sctransform.rds"))


pd = readRDS(paste0(work_path, "/revision/", dataset, "/obj_sctransform_pd.rds"))
p1 = ggplot() +
     geom_point(data = pd, aes(x = UMAP_1, y = UMAP_2), color = "grey80", size = 0.1) +
     geom_point(data = pd[pd$dataset == "monoclonal_RNAseq",], aes(x = UMAP_1, y = UMAP_2), color = "red", size = 0.15) +
     theme_void() +
     theme(legend.position="none") +
     labs(title = "monoclonal_RNAseq") +
     theme(plot.title = element_text(hjust = 0.5))
p2 = ggplot() +
     geom_point(data = pd, aes(x = UMAP_1, y = UMAP_2), color = "grey80", size = 0.1) +
     geom_point(data = pd[pd$dataset == "polyclonal_RNAseq",], aes(x = UMAP_1, y = UMAP_2), color = "blue", size = 0.15) +
     theme_void() +
     theme(legend.position="none") +
     labs(title = "polyclonal_RNAseq") +
     theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("~/share/", "sctransform", ".png"), p1 / p2, width = 7, height = 10, dpi = 300)



###########################################
### Test-3: trying integration using Seurat

source("~/work/scripts/utils.R")
mouse_gene = read.table("~/work/tome/code/refdata-gex-GRCm39-2024-A.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape"
data_path = "/net/shendure/vol10/projects/prime_editing/nobackup/data/2025_11_21_revision"

dataset = "monoclonal_RNAseq"
count_1 = readRDS(paste0(work_path, "/revision/", dataset, "/count.rds"))
pd_1 = readRDS(paste0(work_path, "/revision/", dataset, "/pd.rds"))
pd_1$dataset = dataset

dataset = "polyclonal_RNAseq"
count_2 = readRDS(paste0(work_path, "/revision/", dataset, "/count.rds"))
pd_2 = readRDS(paste0(work_path, "/revision/", dataset, "/pd.rds"))
pd_2$dataset = dataset

count = cbind(count_1, count_2)
pd = rbind(pd_1, pd_2)

mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "M")) &
              mouse_gene$gene_type %in% c("protein_coding", "lncRNA"),]
count = count[rownames(count) %in% as.vector(mouse_gene_sub$gene_ID),]
obj = CreateSeuratObject(count, meta.data = pd)

dataset = "monoclonal_polyclonal_RNAseq"

obj$group = as.vector(obj$dataset)
obj_processed = doClusterSeurat(obj)
saveRDS(obj_processed, paste0(work_path, "/revision/", dataset, "/obj_integration.rds"))

pd = data.frame(obj_processed[[]])
pd$UMAP_1 = as.vector(Embeddings(obj_processed, reduction = "umap")[,1])
pd$UMAP_2 = as.vector(Embeddings(obj_processed, reduction = "umap")[,2])
pd$group = NULL
saveRDS(pd, paste0(work_path, "/revision/", dataset, "/obj_integration_pd.rds"))

cds = new_cell_data_set(count[,colnames(obj_processed)],
                        cell_metadata = data.frame(obj_processed[[]]),
                        gene_metadata = mouse_gene[rownames(count),])
cds = preprocess_cds(cds, use_genes = rownames(count)[1:1000])
reducedDims(cds)$UMAP = as.matrix(Embeddings(obj_processed, reduction = "umap"))
saveRDS(cds, paste0(work_path, "/revision/", dataset, "/cds_integration.rds"))


pd = readRDS(paste0(work_path, "/revision/", dataset, "/obj_integration_pd.rds"))
p1 = ggplot() +
     geom_point(data = pd, aes(x = UMAP_1, y = UMAP_2), color = "grey80", size = 0.1) +
     geom_point(data = pd[pd$dataset == "monoclonal_RNAseq",], aes(x = UMAP_1, y = UMAP_2), color = "red", size = 0.15) +
     theme_void() +
     theme(legend.position="none") +
     labs(title = "monoclonal_RNAseq") +
     theme(plot.title = element_text(hjust = 0.5))
p2 = ggplot() +
     geom_point(data = pd, aes(x = UMAP_1, y = UMAP_2), color = "grey80", size = 0.1) +
     geom_point(data = pd[pd$dataset == "polyclonal_RNAseq",], aes(x = UMAP_1, y = UMAP_2), color = "blue", size = 0.15) +
     theme_void() +
     theme(legend.position="none") +
     labs(title = "polyclonal_RNAseq") +
     theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("~/share/", "integration", ".png"), p1 / p2, width = 7, height = 10, dpi = 300)



##############################################
### Test-4: trying integration using Monocle/3

source("~/work/scripts/utils.R")
mouse_gene = read.table("~/work/tome/code/refdata-gex-GRCm39-2024-A.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape"
data_path = "/net/shendure/vol10/projects/prime_editing/nobackup/data/2025_11_21_revision"

dataset = "monoclonal_RNAseq"
count_1 = readRDS(paste0(work_path, "/revision/", dataset, "/count.rds"))
pd_1 = readRDS(paste0(work_path, "/revision/", dataset, "/pd.rds"))
pd_1$dataset = dataset

dataset = "polyclonal_RNAseq"
count_2 = readRDS(paste0(work_path, "/revision/", dataset, "/count.rds"))
pd_2 = readRDS(paste0(work_path, "/revision/", dataset, "/pd.rds"))
pd_2$dataset = dataset

count = cbind(count_1, count_2)
pd = rbind(pd_1, pd_2)

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
saveRDS(obj_processed, paste0(work_path, "/revision/", dataset, "/obj_batchelor.rds"))

pd = data.frame(obj_processed[[]])
pd$UMAP_1 = as.vector(Embeddings(obj_processed, reduction = "umap")[,1])
pd$UMAP_2 = as.vector(Embeddings(obj_processed, reduction = "umap")[,2])
pd$group = NULL
saveRDS(pd, paste0(work_path, "/revision/", dataset, "/obj_batchelor_pd.rds"))

cds = new_cell_data_set(count[rownames(obj_processed), colnames(obj_processed)],
                        cell_metadata = data.frame(obj_processed[[]]),
                        gene_metadata = mouse_gene[rownames(obj_processed),])
cds = preprocess_cds(cds, use_genes = rownames(obj)[1:1000])
reducedDims(cds)$UMAP = as.matrix(Embeddings(obj_processed, reduction = "umap"))
saveRDS(cds, paste0(work_path, "/revision/", dataset, "/cds_batchelor.rds"))


pd = readRDS(paste0(work_path, "/revision/", dataset, "/obj_batchelor_pd.rds"))
p1 = ggplot() +
     geom_point(data = pd, aes(x = UMAP_1, y = UMAP_2), color = "grey80", size = 0.1) +
     geom_point(data = pd[pd$dataset == "monoclonal_RNAseq",], aes(x = UMAP_1, y = UMAP_2), color = "red", size = 0.15) +
     theme_void() +
     theme(legend.position="none") +
     labs(title = "monoclonal_RNAseq") +
     theme(plot.title = element_text(hjust = 0.5))
p2 = ggplot() +
     geom_point(data = pd, aes(x = UMAP_1, y = UMAP_2), color = "grey80", size = 0.1) +
     geom_point(data = pd[pd$dataset == "polyclonal_RNAseq",], aes(x = UMAP_1, y = UMAP_2), color = "blue", size = 0.15) +
     theme_void() +
     theme(legend.position="none") +
     labs(title = "polyclonal_RNAseq") +
     theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("~/share/", "batchelor", ".png"), p1 / p2, width = 7, height = 10, dpi = 300)













