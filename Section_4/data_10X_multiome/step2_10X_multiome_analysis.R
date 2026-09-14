
#################################################
### Step-2: After QC, analyzing the multiome data

library(Signac)
source("~/work/scripts/utils.R")
work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape"
data_path = "/net/shendure/vol10/projects/akeith/nobackup/analyses/SR_revisions/ES_multiome/count/arc_merged/outs"

mouse_gene = read.table("~/work/tome/code/mouse.v23.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

obj = readRDS(paste0(work_path, "/revision/data_10X_multiome/obj_combine.rds"))

### Gene expression data processing

options(future.globals.maxSize = 100 * 1024^3)

DefaultAssay(obj) = "RNA"
obj = SCTransform(obj)
obj = RunPCA(obj)

### DNA accessibility data processing

DefaultAssay(obj) = "ATAC"
obj = FindTopFeatures(obj, min.cutoff = 5)
obj = RunTFIDF(obj)
obj = RunSVD(obj)

# build a joint neighbor graph using both assays
obj = FindMultiModalNeighbors(
    object = obj,
    reduction.list = list("pca", "lsi"), 
    dims.list = list(1:50, 2:50),
    modality.weight.name = "RNA.weight",
    verbose = TRUE
)

# build a joint UMAP visualization
obj = RunUMAP(
    object = obj,
    nn.name = "weighted.nn",
    assay = "RNA",
    verbose = TRUE
)

saveRDS(obj, paste0(work_path, "/revision/data_10X_multiome/obj_combine_processed.rds"))

p = DimPlot(obj, label = TRUE, repel = TRUE, reduction = "umap", group.by = "lane") + NoLegend()
ggsave("~/share/10X_multiome_UMAP_combine.pdf", p, width = 6, height = 6)

#############
### ATAC only

obj = RunUMAP(object = obj, reduction = 'lsi', dims = 2:50)
p = DimPlot(obj, label = TRUE, repel = TRUE, reduction = "umap", group.by = "lane") + NoLegend()
ggsave("~/share/10X_multiome_UMAP_ATAC.pdf", p, width = 6, height = 6)

############
### RNA only

obj = FindNeighbors(object = obj, dims = 1:50, reduction = "pca")
obj = RunUMAP(object = obj, reduction = 'pca', dims = 1:50)
obj = FindClusters(object = obj, graph.name = "SCT_snn", resolution = 1)

obj_sub = subset(obj, subset = seurat_clusters %in% c(0:13))
p = DimPlot(obj_sub, label = TRUE, repel = TRUE, reduction = "umap", group.by = "seurat_clusters") + NoLegend()
ggsave("~/share/10X_multiome_UMAP_RNA_cluster.pdf", p, width = 6, height = 6)

obj_rna = readRDS(paste0(work_path, "/revision/data_10X_multiome/obj_RNA.rds"))
count = GetAssayData(obj_rna, slot = "counts")[,colnames(obj_sub)]
pd = data.frame(obj_rna[[]])[colnames(obj_sub),]
pd$seurat_clusters = as.vector(obj_sub$seurat_clusters)
obj_rna = CreateSeuratObject(count, meta.data = pd)
Idents(obj_rna) = as.vector(obj_rna$seurat_clusters)
obj_rna_sub = subset(obj_rna, downsample = 1000)
obj_rna_sub = NormalizeData(obj_rna_sub, normalization.method = "LogNormalize", scale.factor = 10000)

res_13 = FindMarkers(obj_rna_sub, ident.1 = 13, only.pos=T)
res_13 = res_13 %>% mutate(gene_ID = rownames(res_13)) %>% 
    left_join(mouse_gene[,c("gene_ID","gene_short_name")], by = "gene_ID")
res_13$cluster = "cluster_13"
res_13 = res_13[res_13$p_val_adj < 0.05,]

res_11 = FindMarkers(obj_rna_sub, ident.1 = 11, only.pos=T)
res_11 = res_11 %>% mutate(gene_ID = rownames(res_11)) %>% 
    left_join(mouse_gene[,c("gene_ID","gene_short_name")], by = "gene_ID")
res_11$cluster = "cluster_11"
res_11 = res_11[res_11$p_val_adj < 0.05,]

res = rbind(res_13, res_11)
write.table(res, "~/share/DEG_cluster_11_13.csv", row.names=F, sep=",", quote=F)


cds = new_cell_data_set(count,
                        cell_metadata = pd,
                        gene_metadata = mouse_gene[rownames(count),])
cds = preprocess_cds(cds, use_genes = sample(rownames(count),1000))
reducedDims(cds)$UMAP = as.matrix(Embeddings(obj_sub, reduction = "umap")[,c("umap_1", "umap_2")])
saveRDS(cds, paste0(work_path, "/revision/data_10X_multiome/cds_rna.rds"))






