
##########################################################
### Processing the transcriptome of monoclonal-gastruloids
### Jan-17-2026

source("~/work/scripts/utils.R")
mouse_gene = read.table("~/work/tome/code/refdata-gex-GRCm39-2024-A.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape"
data_path = "/net/shendure/vol10/projects/prime_editing/nobackup/data/2025_11_21_revision"

dataset = "monoclonal_RNAseq"

run_list = paste0("Mono_L0", 1:9, "_GRCm39")
run_list = c(run_list, paste0("Mono_L1", 0:2, "_GRCm39"))

for(kk in 1:length(run_list)){
    run_id = run_list[kk]
    print(run_id)

    count_matrix = Matrix::readMM(paste0(data_path, "/", run_id, "/outs/raw_feature_bc_matrix/matrix.mtx.gz"))
    df_gene = read.table(paste0(data_path, "/", run_id, "/outs/raw_feature_bc_matrix/features.tsv.gz"), sep="\t", as.is=T)
    df_cell = read.table(paste0(data_path, "/", run_id, "/outs/raw_feature_bc_matrix/barcodes.tsv.gz"), sep="\t", as.is=T)
    colnames(count_matrix) = paste0(run_id, "_", as.vector(df_cell$V1))
    rownames(count_matrix) = as.vector(df_gene$V1)
    count_matrix = count_matrix[df_gene$V3 != "Custom",]

    count_matrix_nonzero = count_matrix
    count_matrix_nonzero@x[count_matrix_nonzero@x >= 1] = 1

    pd = data.frame(cell_id = colnames(count_matrix),
                    UMI_count = colSums(count_matrix),
                    gene_count = colSums(count_matrix_nonzero))
    rownames(pd) = as.vector(pd$cell_id)
    pd$experiment_id = run_id

    pdf(paste0(work_path, "/revision/", dataset, "/processing/", run_id, "/hist_log2_umi_orig.pdf"))
    hist(log2(pd$UMI_count[pd$UMI_count >= 100] + 1), 100)
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

    Matrix::writeMM(t(count_matrix_filter), paste0(work_path, "/revision/", dataset, "/processing/", run_id, "/gene_count.mtx"))
    write.csv(pd_filter, paste0(work_path, "/revision/", dataset, "/processing/", run_id, "/df_cell.csv"))
    write.csv(mouse_gene[rownames(count_matrix_filter),], paste0(work_path, "/revision/", dataset, "/processing/", run_id, "/df_gene.csv"))
}

### running scrublet to identify doublets
### python step1_scrublet.py

for(kk in 1:length(run_list)){
    run_id = run_list[kk]
    print(run_id)

    pd_filter = read.csv(paste0(work_path, "/revision/", dataset, "/processing/", run_id, "/df_cell.csv"), row.names=1)
    scrublet_score = read.csv(paste0(work_path, "/revision/", dataset, "/processing/", run_id, "/doublet_scores_observed_cells.csv"), header=F)
    pd_filter$doublet_score = as.vector(scrublet_score$V1)
    saveRDS(pd_filter, paste0(work_path, "/revision/", dataset, "/processing/", run_id, "/df_cell.rds"))

    simulated_scrublet_score = read.csv(paste0(work_path, "/revision/", dataset, "/processing/", run_id, "/doublet_scores_simulated_doublets.csv"), header=F)
    df = data.frame(score = c(as.vector(scrublet_score$V1), as.vector(simulated_scrublet_score$V1)),
                    group = rep(c("observed", "simulated"), times = c(nrow(scrublet_score), nrow(simulated_scrublet_score))))
    p = ggplot(df, aes(score, fill = group)) + geom_histogram(alpha=0.6, binwidth = 0.02)
    pdf(paste0(work_path, "/revision/", dataset, "/processing/", run_id, "/hist_doublet_score.pdf"))
    print(p)
    dev.off()

    pdf(paste0(work_path, "/revision/", dataset, "/processing/", run_id, "/hist_log2_umi.pdf"))
    print(hist(pd_filter$log2_umi, 100))
    dev.off()

    pdf(paste0(work_path, "/revision/", dataset, "/processing/", run_id, "/hist_MT_pct.pdf"))
    print(hist(pd_filter$MT_pct, 100))
    dev.off()

    pdf(paste0(work_path, "/revision/", dataset, "/processing/", run_id, "/hist_RIBO_pct.pdf"))
    print(hist(pd_filter$RIBO_pct, 100))
    dev.off()
}

### filtering cells

for(kk in 1:length(run_list)){
    run_id = run_list[kk]
    print(run_id)

    pd_filter = readRDS(paste0(work_path, "/revision/", dataset, "/processing/", run_id, "/df_cell.rds"))
    df_gene = read.csv(paste0(work_path, "/revision/", dataset, "/processing/", run_id, "/df_gene.csv"), row.names=1)
    count_matrix_filter = Matrix::readMM(paste0(work_path, "/revision/", dataset, "/processing/", run_id, "/gene_count.mtx"))
    count_matrix_filter = t(count_matrix_filter)
    rownames(count_matrix_filter) = as.vector(df_gene$gene_ID)
    colnames(count_matrix_filter) = as.vector(pd_filter$cell_id)

    pd_filter = pd_filter[pd_filter$log2_umi >= 11 &
                          pd_filter$UMI_count <= quantile(pd_filter$UMI_count, 0.995) &
                          pd_filter$doublet_score <= 0.2 &
                          pd_filter$MT_pct <= 10 &
                          pd_filter$MT_pct >= 1 &
                          pd_filter$RIBO_pct <= 40,]

    count_matrix_filter = count_matrix_filter[,as.vector(pd_filter$cell_id)]

    obj = CreateSeuratObject(count_matrix_filter, meta.data = pd_filter)
    saveRDS(obj, paste0(work_path, "/revision/", dataset, "/processing/", run_id, "/obj_", run_id, ".rds"))
}



### rm df_cell.csv df_gene.csv gene_count.mtx
### mkdir doublet_removing; mv *.pdf *.csv ./doublet_removing



############################################################################################
### Merging eight samples, followed by dimension reduction, to check if any quality issue

source("~/work/scripts/utils.R")
mouse_gene = read.table("~/work/tome/code/refdata-gex-GRCm39-2024-A.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape"
data_path = "/net/shendure/vol10/projects/prime_editing/nobackup/data/2025_11_21_revision"

dataset = "monoclonal_RNAseq"

run_list = paste0("Mono_L0", 1:9, "_GRCm39")
run_list = c(run_list, paste0("Mono_L1", 0:2, "_GRCm39"))

count = NULL
pd = NULL
for(run_id in run_list){
    print(run_id)
    obj = readRDS(paste0(work_path, "/revision/", dataset, "/processing/", run_id, "/obj_", run_id, ".rds"))
    count = cbind(count, GetAssayData(obj, slot = "counts"))
    pd = rbind(pd, data.frame(obj[[]]))
}

obj = CreateSeuratObject(count, meta.data = pd)
saveRDS(obj, paste0(work_path, "/revision/", dataset, "/obj.rds"))
saveRDS(count, paste0(work_path, "/revision/", dataset, "/count.rds"))
saveRDS(pd, paste0(work_path, "/revision/", dataset, "/pd.rds"))

### processing the dataset
mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "M")) &
                                mouse_gene$gene_type %in% c("protein_coding", "lncRNA"),]
count = count[rownames(count) %in% as.vector(mouse_gene_sub$gene_ID),]
obj = CreateSeuratObject(count, meta.data = pd)

obj$group = "monoclonal_RNAseq"
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
saveRDS(cds, paste0(work_path, "/revision/", dataset, "/cds.rds"))








