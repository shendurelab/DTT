
######################################################
### Processing 10X scRNA-seq data of eight gastruloids
### Chengxiang Qiu
### Feb-20, 2025

### processing mGASv3 transcriptome dataset, including two samples, mGas1 and mGas2,
### which are different gastruloids (first lane, mGas1, had 3 more-developed gastruloids and second lane, mGas2, had 5 less-developed gastruloids).
### detecting doublets, filtering cells by UMI count, mito%, and ribo% if necessary

### The script used for processing 10x data

module load cellranger/7.2.0
INPUTFILENAME=mGas1
cellranger count \
--fastqs ./mkfastq_output/outs/fastq_path/AAANT3MHV/,./mkfastq_output_2/outs/fastq_path/AACTKHJM5/ \
--nosecondary \
--localcores 8 \
--include-introns true \
--chemistry SC3Pv3HT \
--sample "$INPUTFILENAME" \
--output-dir ./"$INPUTFILENAME" \
--id "$INPUTFILENAME" \
--transcriptome ./refdata-cellranger-mm10-3.0.0

### Starting processing data

source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

work_path = "Your_work_path"

run_id = "mGas1"

count_matrix = Read10X(paste0(work_path, "/data_mGASv3/data/", run_id, "/outs/raw_feature_bc_matrix"), 
                       gene.column = 1,
                       strip.suffix = T)
colnames(count_matrix) = paste0("mGASv3_", run_id, "_", colnames(count_matrix))

count_matrix_nonzero = count_matrix
count_matrix_nonzero@x[count_matrix_nonzero@x >= 1] = 1

pd = data.frame(cell_id = colnames(count_matrix),
                UMI_count = colSums(count_matrix),
                gene_count = colSums(count_matrix_nonzero))
rownames(pd) = as.vector(pd$cell_id)
pd$experiment_id = paste0("mGASv3_", run_id)

pdf(paste0(work_path, "/data_mGASv3/processing/", run_id, "/hist_log2_umi_orig.pdf"))
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

Matrix::writeMM(t(count_matrix_filter), paste0(work_path, "/data_mGASv3/processing/", run_id, "/gene_count.mtx"))
write.csv(pd_filter, paste0(work_path, "/data_mGASv3/processing/", run_id, "/df_cell.csv"))
write.csv(mouse_gene[rownames(count_matrix_filter),], paste0(work_path, "/data_mGASv3/processing/", run_id, "/df_gene.csv"))

### running scrublet to identify doublets
### python Detect_Doublets.py

scrublet_score = read.csv(paste0(work_path, "/data_mGASv3/processing/", run_id, "/doublet_scores_observed_cells.csv"), header=F)
pd_filter$doublet_score = as.vector(scrublet_score$V1)
saveRDS(pd_filter, paste0(work_path, "/data_mGASv3/processing/", run_id, "/df_cell.rds"))

simulated_scrublet_score = read.csv(paste0(work_path, "/data_mGASv3/processing/", run_id, "/doublet_scores_simulated_doublets.csv"), header=F)
df = data.frame(score = c(as.vector(scrublet_score$V1), as.vector(simulated_scrublet_score$V1)),
                group = rep(c("observed", "simulated"), times = c(nrow(scrublet_score), nrow(simulated_scrublet_score))))
p = ggplot(df, aes(score, fill = group)) + geom_histogram(alpha=0.6, binwidth = 0.02)
pdf(paste0(work_path, "/data_mGASv3/processing/", run_id, "/hist_doublet_score.pdf"))
print(p)
dev.off()

pdf(paste0(work_path, "/data_mGASv3/processing/", run_id, "/hist_log2_umi.pdf"))
hist(pd_filter$log2_umi, 100)
dev.off()

pdf(paste0(work_path, "/data_mGASv3/processing/", run_id, "/hist_MT_pct.pdf"))
hist(pd_filter$MT_pct, 100)
dev.off()

### filtering cells 

if(run_id == "mGas1"){
    umi_cutoff = 12.5
} else {
    umi_cutoff = 12
}

pd_filter = pd_filter[pd_filter$log2_umi >= umi_cutoff &
                      pd_filter$UMI_count <= quantile(pd_filter$UMI_count, 0.995) &
                      pd_filter$doublet_score <= 0.2 &
                      pd_filter$MT_pct <= 10 &
                      pd_filter$MT_pct >= 1 &
                      pd_filter$RIBO_pct <= 40,]

count_matrix_filter = count_matrix_filter[,as.vector(pd_filter$cell_id)]

obj = CreateSeuratObject(count_matrix_filter, meta.data = pd_filter)

saveRDS(obj, paste0(work_path, "/data_mGASv3/processing/", run_id, "/obj_", "mGASv3_", run_id, ".rds"))

### rm df_cell.csv df_gene.csv gene_count.mtx
### mkdir doublet_removing; mv *.pdf *.csv ./doublet_removing

#############################
### making some quality plots

run_id = "mGas1"
df_cell = readRDS(paste0(work_path, "/data_mGASv3/data/", run_id, "/df_cell.rds"))

p1 = ggplot(df_cell, aes(log2_umi)) + geom_histogram(bins = 100) +
    geom_vline(xintercept = c(12.5, quantile(df_cell$log2_umi, 0.995)), color = "red") +
    theme_classic(base_size = 12) +
    labs(x = "Log2 (UMI count per cell)", y = "Count", title = "mGas1 (Well 14/17/25)") +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 

p2 = ggplot(subset(df_cell, MT_pct <= 30), aes(MT_pct)) + geom_histogram(bins = 100) +
    geom_vline(xintercept = c(1, 10), color = "red") +
    theme_classic(base_size = 12) +
    labs(x = "Mito %", y = "Count",  title = "mGas1 (Well 14/17/25)") +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 

run_id = "mGas2"
df_cell = readRDS(paste0(work_path, "/data_mGASv3/data/", run_id, "/df_cell.rds"))

p3 = ggplot(df_cell, aes(log2_umi)) + geom_histogram(bins = 100) +
    geom_vline(xintercept = c(12, quantile(df_cell$log2_umi, 0.995)), color = "red") +
    theme_classic(base_size = 12) +
    labs(x = "Log2 (UMI count per cell)", y = "Count", title = "mGas2 (Well 01/03/21/16/28)") +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 

p4 = ggplot(subset(df_cell, MT_pct <= 30), aes(MT_pct)) + geom_histogram(bins = 100) +
    geom_vline(xintercept = c(1, 10), color = "red") +
    theme_classic(base_size = 12) +
    labs(x = "Mito %", y = "Count", title = "mGas2 (Well 01/03/21/16/28)") +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 

pdf(paste0(work_path, "/data_mGASv3/quality_summary.pdf"))
print(p1 + p2 + p3 + p4)
dev.off()


###################################################################################################
### Step-2: Merging mGas1 and mGas2, followed by dimension reduction, to check if any quality issue

run_list = c("mGas1", "mGas2")

run_id = run_list[1]
obj_1 = readRDS(paste0(work_path, "/data_mGASv3/processing/", run_id, "/obj_", "mGASv3_", run_id, ".rds"))

run_id = run_list[2]
obj_2 = readRDS(paste0(work_path, "/data_mGASv3/processing/", run_id, "/obj_", "mGASv3_", run_id, ".rds"))

obj = merge(obj_1, obj_2)
saveRDS(obj, paste0(work_path, "/data_mGASv3/obj.rds"))

### processing the data
count = GetAssayData(obj, slot = "counts")
pd = data.frame(obj[[]])

mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "M")) &
                                mouse_gene$gene_type %in% c("protein_coding", "lincRNA"),]
count = count[rownames(count) %in% as.vector(mouse_gene_sub$gene_ID),]
obj = CreateSeuratObject(count, meta.data = pd)

obj$group = obj$experiment_id
obj_processed = doClusterSeurat(obj, min.dist = 0.5)
#obj_processed = RunUMAP(object = obj_processed, reduction = "pca", dims = 1:30, min.dist = 0.5, n.components = 2)
obj_processed = FindClusters(object = obj_processed, resolution = 2, verbose = FALSE)
obj_processed = FindClusters(object = obj_processed, resolution = 5, verbose = FALSE)
saveRDS(obj_processed, paste0(work_path, "/data_mGASv3/obj_processed.rds"))

pd = data.frame(obj_processed[[]])
pd$UMAP_1 = as.vector(Embeddings(obj_processed, reduction = "umap")[,1])
pd$UMAP_2 = as.vector(Embeddings(obj_processed, reduction = "umap")[,2])
pd$group = NULL
saveRDS(pd, paste0(work_path, "/data_mGASv3/obj_processed_pd.rds"))



################################################################
### Step-3: manually annotating cell types in the mGASv3 dataset

pd$celltype_old = as.vector(cds$celltype_old)

anno = rep(NA, nrow(pd))
anno[pd$integrated_snn_res.1 %in% c(4,11)] = "PSC-like cells"
anno[pd$integrated_snn_res.1 %in% c(10)] = "Transitional cells"

anno[pd$integrated_snn_res.1 %in% c(2,7)] = "Mesodermal progenitors"
anno[pd$integrated_snn_res.1 %in% c(0,5,13)] = "Somites"

anno[pd$integrated_snn_res.1 %in% c(3,6,8,9)] = "Spinal cord"
anno[pd$integrated_snn_res.1 %in% c(1,12)] = "NMPs"

anno[pd$integrated_snn_res.1 %in% c(14)] = "Definitive endoderm"
anno[pd$integrated_snn_res.1 %in% c(15)] = "Extraembryonic endoderm"
anno[pd$integrated_snn_res.1 %in% c(16)] = "Endothelial cells"

anno[pd$integrated_snn_res.1 %in% c(6) & pd$celltype_old == "NMPs"] = "NMPs"
anno[pd$integrated_snn_res.1 %in% c(13) & pd$integrated_snn_res.5 %in% c(49)] = "Fibroblasts"
anno[pd$integrated_snn_res.1 %in% c(14) & pd$integrated_snn_res.20 %in% c(113)] = "Notochord"

pd$celltype = as.vector(anno)
saveRDS(pd, paste0(work_path, "/data_mGASv3/obj_processed_pd.rds"))


try(ggplot() +
        geom_point(data = pd, aes(x = UMAP_1, y = UMAP_2), size=0.3, color = "black") +
        geom_point(data = pd,
                   aes(x = UMAP_1, y = UMAP_2, color = celltype), size=0.2) +
        ggrepel::geom_text_repel(data = pd %>% group_by(celltype) %>% sample_n(1), aes(x = UMAP_1, y = UMAP_2, label = celltype), color = "black", size = 3, family = "Arial") +
        theme_void() +
        theme(legend.position="none") +
        labs(title = "mGASv3") +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=gastruloid_celltype_color_code) +
        ggsave(paste0(work_path, "/data_mGASv3/UMAP_celltype_title.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)

try(ggplot() +
        geom_point(data = pd, aes(x = UMAP_1, y = UMAP_2), size=0.3, color = "black") +
        geom_point(data = pd,
                   aes(x = UMAP_1, y = UMAP_2, color = celltype), size=0.2) +
       # ggrepel::geom_text_repel(data = pd %>% group_by(celltype) %>% sample_n(1), aes(x = UMAP_1, y = UMAP_2, label = celltype), color = "black", size = 3, family = "Arial") +
        theme_void() +
        theme(legend.position="none") +
       # labs(title = "mGASv3") +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=gastruloid_celltype_color_code) +
        ggsave(paste0(work_path, "/data_mGASv3/UMAP_celltype.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)





