
############################################################################################################
### Processing sci-RNA-seq3 dataset (filtering low quality cells, dimension reduction, cell type annotation)
### Chengxiang Qiu
### Feb-20, 2025

### Of note, we have provided all of the processed data on the GEO, after preprocessing and removing doublets. 
### sci_RNA_seq3.gene_count.mtx
### sci_RNA_seq3.cell_annotation.csv
### sci_RNA_seq3.gene_annotation.csv

### This script is only necessary if you want to reprocess the raw data.
### Read alignment and gene count matrix generation was performed using the pipeline that we developed for sci-RNA-seq3 (https://github.com/JunyueC/sci-RNA-seq3_pipeline)
### After preprocessing, each folder generates a sci_summary.RData file, which can be loaded into R and includes the gene count matrix (gene_count), cell metadata (df_cell), and gene annotation data (df_gene).
### Please contact CX Qiu (cxqiu@uw.edu) if you have any questions about running the pipeline.


####################################
### Step-1: Reading pipeline outputs

source("help_script.R")
mouse_gene = read.table("mouse39-samchoiTAPE.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

data_path = list()
data_path["plate_1"] = "./2023-10-03-sciA-samchoi"
data_path["plate_2"] = "./2023-10-03-sciB-samchoi"
data_path["plate_3"] = "./2023-10-27-sciA-samchoi"
data_path["plate_4"] = "./2023-10-27-sciB-samchoi"
data_path["plate_5"] = "./2023-11-01-sciA-samchoi"
data_path["plate_6"] = "./2023-11-01-sciB-samchoi"

work_path = "Your_work_path"

for(plate_i in 1:6){
    print(plate_i)
    
    load(paste0(data_path[plate_i], "/nobackup/output/report/sci_summary.RData"))
    
    rownames(gene_count) = rownames(df_gene) = df_gene$gene_ID = 
        unlist(lapply(rownames(df_gene), function(x) strsplit(x,"[.]")[[1]][1]))
    df_gene = df_gene %>% select(gene_ID) %>% left_join(mouse_gene %>% select(gene_ID, chr, gene_type, gene_short_name), by = "gene_ID")
    colnames(gene_count) = as.vector(df_cell$sample)
    
    print(dim(gene_count))
    
    df_cell$UMI_count = Matrix::colSums(gene_count)
    gene_count_copy = gene_count
    gene_count_copy@x[gene_count_copy@x > 0] = 1
    df_cell$gene_count = Matrix::colSums(gene_count_copy)
    
    gene_keep = df_gene$chr %in% paste0("chr", c(1:19, "M", "X", "Y"))
    
    keep = df_cell$UMI_count >= 200 & 
        df_cell$gene_count >= 100 & 
        df_cell$unmatched_rate < 0.4
    df_cell = df_cell[keep,]
    df_gene = df_gene[gene_keep,]
    gene_count = gene_count[gene_keep, keep]
    rownames(gene_count) = as.vector(df_gene$gene_ID)
    colnames(gene_count) = rownames(df_cell) = as.vector(df_cell$sample)

    print(dim(gene_count))
    
    saveRDS(gene_count, paste0(work_path, "/data_sci/", "/plate_", plate_i, "/gene_count.rds"))
    saveRDS(df_cell, paste0(work_path, "/data_sci/", "/plate_", plate_i, "/df_cell.rds"))
    
    writeMM(t(gene_count), paste0(work_path, "/data_sci/", "/plate_", plate_i, "/gene_count.mtx"))
    write.csv(df_cell, paste0(work_path, "/data_sci/", "/plate_", plate_i, "/df_cell.csv"))
    rownames(df_gene) = as.vector(df_gene$gene_ID)
    write.csv(df_gene, paste0(work_path, "/data_sci/", "/plate_", plate_i, "/df_gene.csv"))
}


##############################
### Step-2: Detecting doublets

### Applying scrublet to detect doublets
### python Detect_Doublets.py

for(plate_i in 1:6){
    print(plate_i)
    
    df_cell = readRDS(paste0(work_path, "/data_sci/", "/plate_", plate_i, "/df_cell.rds"))
    
    doublet_scores_observed_cells = read.csv(paste0(work_path, "/data_sci/", "/plate_", plate_i, "/doublet_scores_observed_cells.csv"), header=F)
    df_cell$doublet_score = as.vector(doublet_scores_observed_cells$V1)
    df_cell$detected_doublets = df_cell$doublet_score > 0.2
    
    global = read.csv(paste0(work_path, "/data_sci/", "/plate_", plate_i, "/doublet_cluster/global.csv"), header=T)
    main_cluster_list = sort(as.vector(unique(global$louvain)))
    
    res = NULL
    
    for(i in 1:length(main_cluster_list)){
        print(paste0(i, "/", length(main_cluster_list)))
        dat = read.csv(paste0(work_path, "/data_sci/", "/plate_", plate_i, "/doublet_cluster/adata.obs.louvain_", (i-1), ".csv"), header=T)
        print(nrow(dat))
        dat$louvain = as.vector(paste0("cluster_", dat$louvain))
        dat = dat %>%
            left_join(df_cell[,c("sample", "detected_doublets", "doublet_score")], by = "sample")
        
        tmp2 = dat %>%
            group_by(louvain) %>%
            tally() %>%
            dplyr::rename(n_sum = n)
        
        tmp1 = dat %>%
            filter(detected_doublets == "TRUE") %>%
            group_by(louvain) %>%
            tally() %>%
            left_join(tmp2, by = "louvain") %>%
            mutate(frac = n/n_sum) %>%
            filter(frac > 0.15)
        
        dat$doublet_cluster = dat$louvain %in% as.vector(tmp1$louvain) 
        
        p1 = ggplot(dat, aes(umap_1, umap_2, color = louvain)) + geom_point() + theme(legend.position="none") 
        p2 = ggplot(dat, aes(umap_1, umap_2, color = doublet_cluster)) + geom_point()
        p3 = ggplot(dat, aes(umap_1, umap_2, color = detected_doublets)) + geom_point()
        p4 = ggplot(dat, aes(umap_1, umap_2, color = doublet_score)) + geom_point() + scale_color_viridis(option = "plasma")
        p5 = ggplot(dat, aes(umap_1, umap_2, color = UMI_count)) + geom_point() + scale_color_viridis(option = "plasma")
        p6 = ggplot(dat, aes(umap_1, umap_2, color = gene_count)) + geom_point() + scale_color_viridis(option = "plasma")
        
        pdf(paste0(work_path, "/data_sci/", "/plate_", plate_i, "/doublet_cluster/adata.obs.louvain_", (i-1), ".pdf"), 12, 18)
        grid.arrange(p1, p2, p3, p4, p5, p6, nrow=3, ncol=2) 
        dev.off()
        
        dat[,!colnames(dat) %in% c("umap_1", "umap_2")]
        dat$main_louvain = (i-1)
        
        res = rbind(res, dat)
    }
    
    rownames(res) = as.vector(res$sample)
    res = res[rownames(df_cell),]
    df_cell$doublet_cluster = res$doublet_cluster
    
    saveRDS(df_cell, paste0(work_path, "/data_sci/", "/plate_", plate_i, "/df_cell.rds"))

}


##########################################################
### Step-3: Filtering cells based on mito%, ribo%, and UMI

for(plate_i in 1:6){
    print(plate_i)
    
    fd = read.csv(paste0(work_path, "/data_sci/df_gene.csv"), row.names=1)
    rownames(fd) = fd$gene_id = as.vector(fd$gene_ID)
    pd = readRDS(paste0(work_path, "/data_sci/", "/plate_", plate_i, "/df_cell.rds"))
    
    pd$log2_umi = log2(pd$UMI_count)
    pd$EXON_pct = 100 * pd$all_exon / (pd$all_exon + pd$all_intron)
    print(nrow(pd))
    
    ### calculate MT_pct and Ribo_pct per cell
    gene = fd
    MT_gene = as.vector(gene[grep("^mt-",gene$gene_short_name),]$gene_id)
    Rpl_gene = as.vector(gene[grep("^Rpl",gene$gene_short_name),]$gene_id)
    Mrpl_gene = as.vector(gene[grep("^Mrpl",gene$gene_short_name),]$gene_id)
    Rps_gene = as.vector(gene[grep("^Rps",gene$gene_short_name),]$gene_id)
    Mrps_gene = as.vector(gene[grep("^Mrps",gene$gene_short_name),]$gene_id)
    RIBO_gene = c(Rpl_gene, Mrpl_gene, Rps_gene, Mrps_gene)
    
    count = readRDS(paste0(work_path, "/data_sci/", "/plate_", plate_i, "/gene_count.rds"))
    print(sum(colnames(count) != rownames(pd)))
    pd$MT_pct = 100 * Matrix::colSums(count[gene$gene_id %in% MT_gene, ])/Matrix::colSums(count)
    pd$RIBO_pct = 100 * Matrix::colSums(count[gene$gene_id %in% RIBO_gene, ])/Matrix::colSums(count)
    
    print(sum(pd$detected_doublets | pd$doublet_cluster)/nrow(pd))
    pd = pd[!(pd$detected_doublets | pd$doublet_cluster),]
    print(nrow(pd))
    saveRDS(pd, paste0(work_path, "/data_sci/", "/plate_", plate_i, "/pd.rds"))
}

for(plate_i in 1:6){
    print(plate_i)
    
    pd = readRDS(paste0(work_path, "/data_sci/", "/plate_", plate_i, "/pd.rds"))
    
    pd_sub = pd %>% filter(EXON_pct < 85, MT_pct < 10, RIBO_pct < 5, doublet_score < 0.15)
    print(nrow(pd_sub)/nrow(pd))
    
    x1 = mean(pd_sub$log2_umi) - 1.5*sd(pd_sub$log2_umi)
    x2 = mean(pd_sub$log2_umi) + 2*sd(pd_sub$log2_umi)
    pd_sub = pd_sub %>% filter(log2_umi >= x1, log2_umi <= x2)
    pd_sub = pd_sub[pd_sub$UMI_count >= 250,]
    print(nrow(pd_sub)/nrow(pd))
    
    pd = pd[pd$sample %in% pd_sub$sample,]
    
    ### For each cell, we only retain protein-coding genes, lincRNA genes and pseudogenes
    fd = fd[(fd$gene_type %in% c('protein_coding', 'pseudogene', 'lincRNA')) & fd$chr %in% paste0("chr", c(1:19, "M", "X", "Y")),]
    
    count = readRDS(paste0(work_path, "/data_sci/", "/plate_", plate_i, "/gene_count.rds"))
    print(sum(!colnames(count) %in% rownames(pd)))
    count = count[rownames(count) %in% rownames(fd),rownames(pd)]
    
    obj = CreateSeuratObject(count, meta.data = pd)
    print(dim(obj))
    
    saveRDS(obj, paste0(work_path, "/data_sci/", "/plate_", plate_i, "/obj_plate_", plate_i, ".rds"))
}



###########################################################################
### Step-4: Combining six plates and performing regular dimension reduction

count_all = NULL
pd_all = NULL
for(plate_i in 1:6){
    print(plate_i)
    obj_i = readRDS(paste0(work_path, "/data_sci/", "/plate_", plate_i, "/obj_plate_", plate_i, ".rds"))
    
    count = GetAssayData(obj_i, slot = "counts")
    pd = data.frame(obj_i[[]])
    pd$plate = paste0("plate_", plate_i)
    
    count_all = cbind(count_all, count)
    pd_all = rbind(pd_all, pd)
}
mouse_gene_sub = mouse_gene[rownames(count_all),]

### Of note, these data have been provided from the GEO, including 21,970 genes x 247,064 cells

writeMM(count_all, paste0(work_path, "/data_sci/", "sci_RNA_seq3.gene_count.mtx"))
write.csv(pd_all[,c("UMI_count","gene_count","doublet_score","MT_pct","RIBO_pct","plate")], paste0(work_path, "/data_sci/", "sci_RNA_seq3.cell_annotation.csv"))
write.csv(mouse_gene_sub[,c("chr","start","end","gene_type","gene_short_name")], paste0(work_path, "/data_sci/", "/sci_RNA_seq3.gene_annotation.csv"))


df_gene_sub = mouse_gene_sub[!mouse_gene_sub$chr %in% c("chrX", "chrY"),]

for(plate_i in 1:6){
    obj_i = readRDS(paste0(work_path, "/data_sci/", "/plate_", plate_i, "/obj_plate_", plate_i, ".rds"))
    
    count = GetAssayData(obj_i, slot = "counts")
    count = count[rownames(count) %in% rownames(df_gene_sub),]
    obj_i = CreateSeuratObject(count, meta.data = data.frame(obj_i[[]]))
    obj_i$plate = paste0("plate_", plate_i)
    
    if(plate_i == 1){
        obj = obj_i
    } else {
        obj = merge(obj, obj_i)
    }
}

obj$group = "gastruloid"
obj_processed = doClusterSeurat(obj)
obj_processed = FindClusters(object = obj_processed, resolution = 2)
obj_processed = FindClusters(object = obj_processed, resolution = 5)

obj_processed$UMAP_2d_1 = Embeddings(obj_processed, reduction = "umap")[,1]
obj_processed$UMAP_2d_2 = Embeddings(obj_processed, reduction = "umap")[,2]
obj_processed = RunUMAP(object = obj_processed, 
                        reduction = "pca", 
                        dims = 1:30, 
                        min.dist = 0.3, 
                        n.components = 3)

obj_processed$UMAP_1 = Embeddings(obj_processed, reduction = "umap")[,1]
obj_processed$UMAP_2 = Embeddings(obj_processed, reduction = "umap")[,2]
obj_processed$UMAP_3 = Embeddings(obj_processed, reduction = "umap")[,3]

saveRDS(obj_processed, paste0(work_path, "/data_sci/obj_processed.rds"))
saveRDS(data.frame(obj_processed[[]]), paste0(work_path, "/data_sci/obj_processed_pd.rds"))

df = data.frame(obj_processed[[]])
anno = rep(NA, nrow(df))
anno[df$RNA_snn_res.1 %in% c(0,3,9,24)] = "Hindbrain"
anno[df$RNA_snn_res.1 %in% c(1,4,13)] = "Spinal cord"
anno[df$RNA_snn_res.1 %in% c(2,5,10,18)] = "Somites"
anno[df$RNA_snn_res.1 %in% c(6)] = "NMPs"
anno[df$RNA_snn_res.1 %in% c(7,11)] = "PSC-like cells"
anno[df$RNA_snn_res.1 %in% c(8,14,16)] = "Transitional cells"
anno[df$RNA_snn_res.1 %in% c(12)] = "Early neurons"
anno[df$RNA_snn_res.1 %in% c(15)] = "Node-like cells"
anno[df$RNA_snn_res.1 %in% c(17)] = "Mesodermal progenitors"
anno[df$RNA_snn_res.1 %in% c(19)] = "Notochord"
anno[df$RNA_snn_res.1 %in% c(20)] = "Extraembryonic endoderm"
anno[df$RNA_snn_res.1 %in% c(21)] = "Cardiac mesoderm"
anno[df$RNA_snn_res.1 %in% c(22)] = "Definitive endoderm"
anno[df$RNA_snn_res.1 %in% c(23)] = "Endothelial cells"
anno[df$RNA_snn_res.2 %in% c(32)] = "Anterior mesendoderm"
anno[df$RNA_snn_res.5 %in% c(38)] = "Floor plate"
df$celltype = as.vector(anno)
saveRDS(df, paste0(work_path, "/data_sci/obj_processed_pd.rds"))



