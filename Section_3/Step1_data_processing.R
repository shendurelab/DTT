
####################################################
### Processing 10X scRNA-seq data of "tree of trees"
### Chengxiang Qiu
### Feb-20, 2025

### processing mGASv5 dataset, including eight lanes, mGASv5_1-8
### we dissociated all gastruloids together and distributed across 8 lanes.
### detecting doublets, filtering cells by UMI count, mito%, and ribo% if necessary

### The script used for processing 10x data
module load cellranger/7.2.0
#!/bin/bash
INPUTFILES=(1 2 3 4 5 6 7 8)
INPUTFILENAME="${INPUTFILES[$SGE_TASK_ID - 1]}"
echo $INPUTFILENAME
cellranger count \
--fastqs ./mkfastq_output/outs/fastq_path/AAC7MYYHV \
--nosecondary \
--localcores 8 \
--include-introns true \
--sample mGASv4_"$INPUTFILENAME" \
--output-dir ./mGASv5_"$INPUTFILENAME" \
--id mGASv5_"$INPUTFILENAME" \
--transcriptome ./refdata-cellranger-mm10-3.0.0

###################################################
### Step-1: Starting to read the data matrix into R

source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

work_path = "Your_work_path"

run_id = "mGASv5_8"

count_matrix = Read10X(paste0(work_path, "/data_mGASv5/data/", run_id, "/outs/raw_feature_bc_matrix"), 
                       gene.column = 1,
                       strip.suffix = T)
colnames(count_matrix) = paste0(run_id, "_", colnames(count_matrix))

count_matrix_nonzero = count_matrix
count_matrix_nonzero@x[count_matrix_nonzero@x >= 1] = 1

pd = data.frame(cell_id = colnames(count_matrix),
                UMI_count = colSums(count_matrix),
                gene_count = colSums(count_matrix_nonzero))
rownames(pd) = as.vector(pd$cell_id)
pd$experiment_id = run_id

pdf(paste0(work_path, "/data_mGASv5/processing/", run_id, "/hist_log2_umi_orig.pdf"))
hist(log2(pd$UMI_count + 1), 100)
dev.off()

mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "X", "Y", "M")),]

count_matrix_filter = count_matrix[rownames(count_matrix) %in% as.vector(mouse_gene_sub$gene_ID),
                                   pd$UMI_count >= 500 & pd$gene_count >= 250]
pd_filter = pd[pd$UMI_count >= 500 & pd$gene_count >= 250,]
pd_filter$log2_umi = log2(pd_filter$UMI_count)

MT_gene = as.vector(mouse_gene[grep("^mt-",mouse_gene$gene_short_name),]$gene_ID)
Rpl_gene = as.vector(mouse_gene[grep("^Rpl",mouse_gene$gene_short_name),]$gene_ID)
Mrpl_gene = as.vector(mouse_gene[grep("^Mrpl",mouse_gene$gene_short_name),]$gene_ID)
Rps_gene = as.vector(mouse_gene[grep("^Rps",mouse_gene$gene_short_name),]$gene_ID)
Mrps_gene = as.vector(mouse_gene[grep("^Mrps",mouse_gene$gene_short_name),]$gene_ID)
RIBO_gene = c(Rpl_gene, Mrpl_gene, Rps_gene, Mrps_gene)

pd_filter$MT_pct = 100 * Matrix::colSums(count_matrix_filter[rownames(count_matrix_filter) %in% MT_gene,])/Matrix::colSums(count_matrix_filter)
pd_filter$RIBO_pct = 100 * Matrix::colSums(count_matrix_filter[rownames(count_matrix_filter) %in% RIBO_gene,])/Matrix::colSums(count_matrix_filter)

Matrix::writeMM(t(count_matrix_filter), paste0(work_path, "/data_mGASv5/processing/", run_id, "/gene_count.mtx"))
write.csv(pd_filter, paste0(work_path, "/data_mGASv5/processing/", run_id, "/df_cell.csv"))
write.csv(mouse_gene[rownames(count_matrix_filter),], paste0(work_path, "/data_mGASv5/processing/", run_id, "/df_gene.csv"))

### running scrublet to identify doublets
### python Detect_Doublets.py

scrublet_score = read.csv(paste0(work_path, "/data_mGASv5/processing/", run_id, "/doublet_scores_observed_cells.csv"), header=F)
pd_filter$doublet_score = as.vector(scrublet_score$V1)
saveRDS(pd_filter, paste0(work_path, "/data_mGASv5/processing/", run_id, "/df_cell.rds"))

simulated_scrublet_score = read.csv(paste0(work_path, "/data_mGASv5/processing/", run_id, "/doublet_scores_simulated_doublets.csv"), header=F)
df = data.frame(score = c(as.vector(scrublet_score$V1), as.vector(simulated_scrublet_score$V1)),
                group = rep(c("observed", "simulated"), times = c(nrow(scrublet_score), nrow(simulated_scrublet_score))))
p = ggplot(df, aes(score, fill = group)) + geom_histogram(alpha=0.6, binwidth = 0.02)
pdf(paste0(work_path, "/data_mGASv5/processing/", run_id, "/hist_doublet_score.pdf"))
print(p)
dev.off()

pdf(paste0(work_path, "/data_mGASv5/processing/", run_id, "/hist_log2_umi.pdf"))
hist(pd_filter$log2_umi, 100)
dev.off()

pdf(paste0(work_path, "/data_mGASv5/processing/", run_id, "/hist_MT_pct.pdf"))
hist(pd_filter$MT_pct, 100)
dev.off()

### filtering cells 

if(run_id %in% paste0("mGASv5_", c(1,3,4,5,7,8))){
    umi_cutoff = 10.5
} else if(run_id %in% paste0("mGASv5_", c(2,6))) {
    umi_cutoff = 11
} else {
    umi_cutoff = NULL
}
print(umi_cutoff)

pd_filter = pd_filter[pd_filter$log2_umi >= umi_cutoff &
                      pd_filter$UMI_count <= quantile(pd_filter$UMI_count, 0.995) &
                      pd_filter$doublet_score <= 0.2 &
                      pd_filter$MT_pct <= 10 &
                      pd_filter$MT_pct >= 1 &
                      pd_filter$RIBO_pct <= 40,]

count_matrix_filter = count_matrix_filter[,as.vector(pd_filter$cell_id)]

obj = CreateSeuratObject(count_matrix_filter, meta.data = pd_filter)

saveRDS(obj, paste0(work_path, "/data_mGASv5/processing/", run_id, "/obj_", run_id, ".rds"))


### rm df_cell.csv df_gene.csv gene_count.mtx
### mkdir doublet_removing; mv *.pdf *.csv ./doublet_removing




#################################################################################################
### Step-2: Merging eight samples, followed by dimension reduction, to check if any quality issue

run_list = paste0("mGASv5_", c(1:8))

run_id = run_list[1]
obj_1 = readRDS(paste0(work_path, "/data_mGASv5/processing/", run_id, "/obj_", run_id, ".rds"))
run_id = run_list[2]
obj_2 = readRDS(paste0(work_path, "/data_mGASv5/processing/", run_id, "/obj_", run_id, ".rds"))
run_id = run_list[3]
obj_3 = readRDS(paste0(work_path, "/data_mGASv5/processing/", run_id, "/obj_", run_id, ".rds"))
run_id = run_list[4]
obj_4 = readRDS(paste0(work_path, "/data_mGASv5/processing/", run_id, "/obj_", run_id, ".rds"))
run_id = run_list[5]
obj_5 = readRDS(paste0(work_path, "/data_mGASv5/processing/", run_id, "/obj_", run_id, ".rds"))
run_id = run_list[6]
obj_6 = readRDS(paste0(work_path, "/data_mGASv5/processing/", run_id, "/obj_", run_id, ".rds"))
run_id = run_list[7]
obj_7 = readRDS(paste0(work_path, "/data_mGASv5/processing/", run_id, "/obj_", run_id, ".rds"))
run_id = run_list[8]
obj_8 = readRDS(paste0(work_path, "/data_mGASv5/processing/", run_id, "/obj_", run_id, ".rds"))

obj = merge(obj_1, y = c(obj_2, obj_3, obj_4, obj_5, obj_6, obj_7, obj_8))
saveRDS(obj, paste0(work_path, "/data_mGASv5/obj.rds"))

count = GetAssayData(obj, slot = "counts")
pd = data.frame(obj[[]])

saveRDS(count, paste0(work_path, "/data_mGASv5/count.rds"))
saveRDS(pd, paste0(work_path, "/data_mGASv5/pd.rds"))

### processing the data
count = GetAssayData(obj, slot = "counts")
pd = data.frame(obj[[]])

mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "M")) &
                                mouse_gene$gene_type %in% c("protein_coding", "lincRNA"),]
count = count[rownames(count) %in% as.vector(mouse_gene_sub$gene_ID),]
obj = CreateSeuratObject(count, meta.data = pd)

obj$group = "mGASv5"
obj_processed = doClusterSeurat(obj)
obj_processed = FindClusters(object = obj_processed, resolution = 2, verbose = FALSE)
obj_processed = FindClusters(object = obj_processed, resolution = 5, verbose = FALSE)
saveRDS(obj_processed, paste0(work_path, "/data_mGASv5/obj_processed.rds"))

pd = data.frame(obj_processed[[]])
pd$UMAP_1 = as.vector(Embeddings(obj_processed, reduction = "umap")[,1])
pd$UMAP_2 = as.vector(Embeddings(obj_processed, reduction = "umap")[,2])
pd$group = NULL
saveRDS(pd, paste0(work_path, "/data_mGASv5/obj_processed_pd.rds"))



#######################################################
### Step-3: Performing SCTransform for normalizing data

count = readRDS(paste0(work_path, "/data_mGASv5/count.rds"))
pd = readRDS(paste0(work_path, "/data_mGASv5/pd.rds"))

mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "M")) &
                                mouse_gene$gene_type %in% c("protein_coding", "lincRNA"),]
count = count[rownames(count) %in% as.vector(mouse_gene_sub$gene_ID),]
obj = CreateSeuratObject(count, meta.data = pd)
print(dim(obj))

obj = SCTransform(object = obj, verbose = FALSE)
obj = RunPCA(object = obj, npcs = 30, verbose = FALSE)
obj = RunUMAP(object = obj, dims = 1:30, min.dist = 0.3, verbose = FALSE)
obj = FindNeighbors(object = obj, dims = 1:30, verbose = FALSE, reduction = "pca")
obj = FindClusters(object = obj, resolution = 1, verbose = FALSE)
obj = FindClusters(object = obj, resolution = 2, verbose = FALSE)
obj = FindClusters(object = obj, resolution = 5, verbose = FALSE)
saveRDS(obj, paste0(work_path, "/data_mGASv5/obj_sctransform.rds"))

pd = data.frame(obj[[]])
pd$UMAP_1 = as.vector(Embeddings(obj, reduction = "umap")[,1])
pd$UMAP_2 = as.vector(Embeddings(obj, reduction = "umap")[,2])
pd$group = NULL
saveRDS(pd, paste0(work_path, "/data_mGASv5/obj_sctransform_pd.rds"))



#################################
### Step-4: Annotating cell types


pd = readRDS(paste0(work_path, "/data_mGASv5/obj_sctransform_pd.rds"))

anno = rep(NA, nrow(pd))
anno[pd$SCT_snn_res.1 %in% c(1,7,8,16,23)] = "PSC-like cells"
anno[pd$SCT_snn_res.1 %in% c(28)] = "ESCs (2-cell state)"
anno[pd$SCT_snn_res.1 %in% c(6,20,22,24)] = "Transitional cells"

anno[pd$SCT_snn_res.1 %in% c(13)] = "Mesodermal progenitors"
anno[pd$SCT_snn_res.1 %in% c(0,10)] = "Somites"
anno[pd$SCT_snn_res.1 %in% c(27)] = "Endothelial cells"

anno[pd$SCT_snn_res.1 %in% c(2,11)] = "NMPs"
anno[pd$SCT_snn_res.1 %in% c(3,5,25)] = "Spinal cord"
anno[pd$SCT_snn_res.1 %in% c(4,9,12,17)] = "Hindbrain"
anno[pd$SCT_snn_res.1 %in% c(18,19)] = "Early neurons"
anno[pd$SCT_snn_res.1 %in% c(26)] = "Motor neurons"

anno[pd$SCT_snn_res.1 %in% c(21)] = "Node-like cells"
anno[pd$SCT_snn_res.1 %in% c(15)] = "Notochord"
anno[pd$SCT_snn_res.1 %in% c(14)] = "Definitive endoderm"

anno[pd$SCT_snn_res.2 %in% c(32)] = "MHB"
anno[pd$SCT_snn_res.5 %in% c(61)] = "Floor plate"
anno[pd$SCT_snn_res.5 %in% c(76)] = "Extraembryonic endoderm"
anno[pd$SCT_snn_res.5 %in% c(68)] = "Fibroblasts"
anno[pd$SCT_snn_res.1 %in% c(6) & !pd$SCT_snn_res.5 %in% c(52,33)] = "Hindbrain"

pd$celltype = as.vector(anno)
saveRDS(pd, paste0(work_path, "/data_mGASv5/obj_sctransform_pd.rds"))

try(ggplot() +
        geom_point(data = pd, aes(x = UMAP_1, y = UMAP_2), size=0.15, color = "black") +
        geom_point(data = pd,
                   aes(x = UMAP_1, y = UMAP_2, color = celltype), size=0.1) +
        ggrepel::geom_text_repel(data = pd %>% group_by(celltype) %>% sample_n(1), aes(x = UMAP_1, y = UMAP_2, label = celltype), color = "black", size = 3, family = "Arial") +
        theme_void() +
        theme(legend.position="none") +
        labs(title = "mGASv5") +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=gastruloid_celltype_color_code) +
        ggsave(paste0(work_path, "/data_mGASv5/UMAP_celltype_title.png"),
               dpi = 300,
               height  = 6, 
               width = 6), silent = TRUE)

### One cluster (SCT_snn_res.1 = 20) of Transitional cells is a little messy, so I put it in a lower layer when plotting
try(ggplot() +
        geom_point(data = pd, aes(x = UMAP_1, y = UMAP_2), size=0.15, color = "black") +
        geom_point(data = pd[pd$SCT_snn_res.1 == 20,],
                   aes(x = UMAP_1, y = UMAP_2, color = celltype), size=0.1) +
        geom_point(data = pd[pd$SCT_snn_res.1 != 20,],
                   aes(x = UMAP_1, y = UMAP_2, color = celltype), size=0.1) +
        theme_void() +
        theme(legend.position="none") +
        # labs(title = "mGASv5") +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=gastruloid_celltype_color_code) +
        ggsave(paste0(work_path, "/data_mGASv5/UMAP_celltype.png"),
               dpi = 300,
               height  = 6, 
               width = 6), silent = TRUE)











