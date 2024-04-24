
library(Matrix)
library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra) 
library(viridis)

mouse_gene = read.table("~/work/tome/code/mouse39-samchoiTAPE.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

data_path = list()
data_path["plate_1"] = "/net/shendure/vol10/projects/gastruloid-TAPE-sci/2023-10-03-sciA-samchoi"
data_path["plate_2"] = "/net/shendure/vol10/projects/gastruloid-TAPE-sci/2023-10-03-sciB-samchoi"
data_path["plate_3"] = "/net/shendure/vol10/projects/gastruloid-TAPE-sci/2023-10-27-sciA-samchoi"
data_path["plate_4"] = "/net/shendure/vol10/projects/gastruloid-TAPE-sci/2023-10-27-sciB-samchoi"
data_path["plate_5"] = "/net/shendure/vol10/projects/gastruloid-TAPE-sci/2023-11-01-sciA-samchoi"
data_path["plate_6"] = "/net/shendure/vol10/projects/gastruloid-TAPE-sci/2023-11-01-sciB-samchoi"

work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape"

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


### run scrublet using python to detect doublets
### python step2_run_scrublet.py

###########################################
### after running scrublet using python ###
###########################################

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
    
    ### sum(df_cell$detected_doublets | df_cell$doublet_cluster) = 18425
    ### sum(df_cell$detected_doublets | df_cell$doublet_cluster)/nrow(df_cell) = 0.03150725
    saveRDS(df_cell, paste0(work_path, "/data_sci/", "/plate_", plate_i, "/df_cell.rds"))

}

### plate_1: 11.12%
### plate_2: 8.24%
### plate_3: 7.36%
### plate_4: 8.27%
### plate_5: 10.23%
### plate_6: 10.37%

#########################
### main cell cluster ###
#########################


for(plate_i in 1:6){
    print(plate_i)
    
    fd = read.csv(paste0(work_path, "/data_sci/df_gene.csv"), row.names=1)
    rownames(fd) = fd$gene_id = as.vector(fd$gene_ID)
    pd = readRDS(paste0(work_path, "/data_sci/", "/plate_", plate_i, "/df_cell.rds"))
    
    pd$log2_umi = log2(pd$UMI_count)
    pd$EXON_pct = 100 * pd$all_exon / (pd$all_exon + pd$all_intron)
    print(nrow(pd))
    ### n = 584,786 cells
    
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

library(Seurat)
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


###############################
### check cutoff of UMI #######
###############################

plate_id = "3"

setwd(paste0("~/work/sam_tape/data_sci/plate_", plate_id))
library(dplyr)
library(ggplot2)
library(gridExtra) 

pd = readRDS("pd.rds")
print(dim(pd))

p1 = ggplot(pd, aes(log2_umi)) + geom_histogram(binwidth = 0.1) 

p2 = ggplot(pd, aes(RIBO_pct)) + geom_histogram(binwidth = 0.1)

p3 = ggplot(pd, aes(MT_pct)) + geom_histogram(binwidth = 0.1) 

p4 = ggplot(pd, aes(EXON_pct)) + geom_histogram(binwidth = 0.1) + geom_vline(xintercept = 85) 

pdf("quality_summary.pdf",6,5)
grid.arrange(p1, p2, p3, p4, nrow=2, ncol=2) 
dev.off()

pd_sub = pd %>% filter(EXON_pct < 85, MT_pct < 10, RIBO_pct < 5, doublet_score < 0.15)
print(nrow(pd_sub)/nrow(pd))

x1 = mean(pd_sub$log2_umi) - 1.5*sd(pd_sub$log2_umi)
x2 = mean(pd_sub$log2_umi) + 2*sd(pd_sub$log2_umi)

x1 = max(x1, log2(250))

hist(pd_sub$log2_umi, 500); abline(v = x1); abline(v = x2)
print(sum(pd_sub$log2_umi >= x1 & pd_sub$log2_umi <= x2)/nrow(pd_sub))


#######################################################################
### Combining six plates and performing regular dimension reduction ###
#######################################################################

df_gene_sub = read.csv(paste0(work_path, "/data_sci/df_gene_sub.csv"), row.names=1)
df_gene_sub = df_gene_sub[!df_gene_sub$chr %in% c("chrX", "chrY"),]

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

source("~/work/scripts/tome/utils.R")
work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape"
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

df = readRDS(paste0(work_path, "/data_sci/obj_processed_pd.rds"))
pd_gas = readRDS(paste0(work_path, "/data_sci/pd_gastruloid.rds"))
sum(rownames(pd_gas) %in% rownames(df))
pd_gas = pd_gas[rownames(df),]
df$pijuan_celltype = as.vector(pd_gas$pijuan_celltype)
df$jax_major_trajectory = as.vector(pd_gas$jax_major_trajectory)
df$jax_sub_trajectory = as.vector(pd_gas$jax_sub_trajectory)
saveRDS(df, paste0(work_path, "/data_sci/obj_processed_pd.rds"))

save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup"
fig = plot_ly(df, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~RNA_snn_res.1)
saveWidget(fig, paste0(save_path, "/sam_tape/gastruloid_sci_RNA_snn_res.1.html"), selfcontained = FALSE, libdir = "tmp")
fig = plot_ly(df, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~RNA_snn_res.5)
saveWidget(fig, paste0(save_path, "/sam_tape/gastruloid_sci_RNA_snn_res.5.html"), selfcontained = FALSE, libdir = "tmp")


anno = rep(NA, nrow(df))
anno[df$RNA_snn_res.1 %in% c(0,3,9,24)] = "Hindbrain"
anno[df$RNA_snn_res.1 %in% c(1,4,13)] = "Spinal cord"
anno[df$RNA_snn_res.1 %in% c(2,5,10,18)] = "Somites"
anno[df$RNA_snn_res.1 %in% c(6)] = "NMPs"
anno[df$RNA_snn_res.1 %in% c(7,11)] = "Epiblast"
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

df$UMAP_3d_1 = df$UMAP_1
df$UMAP_3d_2 = df$UMAP_2
df$UMAP_3d_3 = df$UMAP_3
df$UMAP_1 = df$UMAP_2d_1
df$UMAP_2 = df$UMAP_2d_2
df$UMAP_3 = NULL
saveRDS(df, paste0(work_path, "/data_sci/obj_processed_pd.rds"))

fig = plot_ly(df, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~celltype)
saveWidget(fig, paste0(save_path, "/sam_tape/gastruloid_sci_celltype.html"), selfcontained = FALSE, libdir = "tmp")



pd = df

try(ggplot() +
        geom_point(data = pd, aes(x = UMAP_1, y = UMAP_2), size=0.15, color = "black") +
        geom_point(data = pd,
                   aes(x = UMAP_1, y = UMAP_2, color = celltype), size=0.1) +
        ggrepel::geom_text_repel(data = pd %>% group_by(celltype) %>% sample_n(1), aes(x = UMAP_1, y = UMAP_2, label = celltype), color = "black", size = 3, family = "Arial") +
        theme_void() +
        theme(legend.position="none") +
        labs(title = "mGASv4") +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=gastruloid_celltype_color_code) +
        ggsave(paste0(work_path, "/data_sci/UMAP_celltype_title.png"),
               dpi = 300,
               height  = 5, 
               width = 6), silent = TRUE)

try(ggplot() +
        geom_point(data = pd, aes(x = UMAP_1, y = UMAP_2), size=0.15, color = "black") +
        geom_point(data = pd,
                   aes(x = UMAP_1, y = UMAP_2, color = celltype), size=0.1) +
        # ggrepel::geom_text_repel(data = pd %>% group_by(celltype) %>% sample_n(1), aes(x = UMAP_1, y = UMAP_2, label = celltype), color = "black", size = 3, family = "Arial") +
        theme_void() +
        theme(legend.position="none") +
        # labs(title = "mGASv4") +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=gastruloid_celltype_color_code) +
        ggsave(paste0(work_path, "/data_sci/UMAP_celltype.png"),
               dpi = 300,
               height  = 5, 
               width = 6), silent = TRUE)




### plot gene expression if necessary

setwd("~/work/sam_tape/data_sci/")
source("~/work/scripts/tome/utils.R")

cds = readRDS("cds_processed.rds")
pd = readRDS("obj_processed_pd.rds")
reducedDims(cds)$UMAP = as.matrix(pd[colnames(cds),c("UMAP_2d_1", "UMAP_2d_2")])

### Cardiac mesoderm
my_plot_cells(cds, genes = c("Hand2","Tbx20"), how_many_rows = 1) 

### Definitive endoderm
my_plot_cells(cds, genes = c("Sox17", "Trh"), how_many_rows = 1) 

### Early neurons
my_plot_cells(cds, genes = c("Ebf2", "Neurod1"), how_many_rows = 1) 

### Endothelial cells
my_plot_cells(cds, genes = c("Kdr", "Cdh5"), how_many_rows = 1) 

### Epiblast/PGC
my_plot_cells(cds, genes = c("Nanog","Utf1"), how_many_rows = 1) 

### Floor plate
my_plot_cells(cds, genes = c("Foxa2","Shh"), how_many_rows = 1) 

### Node-like cells
my_plot_cells(cds, genes = c("Kcnip4"), how_many_rows = 1) 

### Hindbrain
my_plot_cells(cds, genes = c("Egr2", "Pax5"), how_many_rows = 1) 

### Mesodermal cells (Eomes+)
my_plot_cells(cds, genes = c("Eomes", "Lhx1"), how_many_rows = 1) 

### Mesodermal progenitors (Tbx6+)
my_plot_cells(cds, genes = c("Tbx6", "Hes7"), how_many_rows = 1) 

### NMPs
my_plot_cells(cds, genes = c("T","Cdx2"), how_many_rows = 1) 

### Notochord
my_plot_cells(cds, genes = c("Noto"), how_many_rows = 1) 

### Parietal endoderm
my_plot_cells(cds, genes = c("Lamb1", "Sparc"), how_many_rows = 1) 

### Somites
my_plot_cells(cds, genes = c("Tcf15", "Meox1"), how_many_rows = 1) 

### Spinal cord
my_plot_cells(cds, genes = c("Hoxb4","Hoxc6"), how_many_rows = 1) 

### Transitional cells
my_plot_cells(cds, genes = c("Bhlhe41","Atp6v0b"), how_many_rows = 1) 


gene_list = unique(c("Hand2","Tbx20", "Sox17", "Trh", 
                     "Ebf2", "Neurod1", "Kdr", "Cdh5",
                     "Nanog","Utf1","Foxa2","Shh", "Kcnip4",
                     "Egr2", "Pax5", "Eomes", "Lhx1", 
                     "Tbx6", "Hes7", "T","Cdx2", "Noto",
                     "Lamb1", "Sparc", "Tcf15", "Meox1", "Hoxb4","Hoxc6","Bhlhe41","Atp6v0b"))

my_plot_cells(cds, genes = gene_list, how_many_rows = 5, cell_size = 0.5) + 
    theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
    ggsave("./plot/Gastruloid_sci_1.png", dpi = 300, height = 10, width = 14)

my_plot_cells(cds, genes = gene_list, how_many_rows = 5) + 
    ggsave("./plot/Gastruloid_sci_2.png", dpi = 300, height = 10, width = 15)


###







#################### BACKUP #############################################


#############################
### How about if we align_cds by log2_umi


source("~/work/scripts/tome/utils.R")
work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape"
obj_processed = readRDS(paste0(work_path, "/data_sci/obj_processed.rds"))
obj = CreateSeuratObject(GetAssayData(obj_processed, slot = "counts"),
                         meta.data = data.frame(obj_processed[[]]))
obj$group = obj$RNA_snn_res.1 = obj$seurat_clusters = NULL
obj$UMAP_2d_1 = obj$UMAP_2d_2 = obj$UMAP_1 = obj$UMAP_2 = obj$UMAP_3 = NULL

obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2500)
obj <- ScaleData(object = obj, verbose = FALSE)
obj <- RunPCA(object = obj, npcs = 30, verbose = FALSE)

pca_coor = as.matrix(Embeddings(obj, reduction = "pca"))

if(!"log2_umi" %in% colnames(obj[[]])){
    if(!"UMI_count" %in% colnames(obj[[]])){
        obj$log2_umi = log2(Matrix::colSums(GetAssayData(obj, slot = "counts")))
    } else {
        obj$log2_umi = log2(obj$UMI_count)
    }
}

set.seed(2016)
residual_model_formula_str = "~log2_umi"
X.model_mat <- Matrix::sparse.model.matrix(stats::as.formula(residual_model_formula_str), 
                                           data = data.frame(obj[[]]), drop.unused.levels = TRUE)
fit <- limma::lmFit(Matrix::t(pca_coor), X.model_mat)

beta <- fit$coefficients[, -1, drop = FALSE]
beta[is.na(beta)] <- 0
aligned_coor <- Matrix::t(as.matrix(Matrix::t(pca_coor)) - 
                              beta %*% Matrix::t(X.model_mat[, -1]))

colnames(aligned_coor) = paste0("aligned_", 1:30)
rownames(aligned_coor) = colnames(obj)

obj[['aligned']] = Seurat::CreateDimReducObject(embeddings=as.matrix(aligned_coor), key='aligned_')

obj = FindNeighbors(object = obj, dims = 1:30, reduction = "aligned")
obj = FindClusters(object = obj, resolution = 1)
obj = FindClusters(object = obj, resolution = 2)

obj = RunUMAP(object = obj, reduction = "aligned", dims = 1:30, min.dist = 0.3, n.components = 2)
obj$UMAP_2d_1 = Embeddings(obj, reduction = "umap")[,1]
obj$UMAP_2d_2 = Embeddings(obj, reduction = "umap")[,2]

obj = RunUMAP(object = obj, reduction = "aligned", dims = 1:30, min.dist = 0.3, n.components = 3)
obj$UMAP_1 = Embeddings(obj, reduction = "umap")[,1]
obj$UMAP_2 = Embeddings(obj, reduction = "umap")[,2]
obj$UMAP_3 = Embeddings(obj, reduction = "umap")[,3]
saveRDS(obj, paste0(work_path, "/data_sci/obj_aligned.rds"))
saveRDS(data.frame(obj[[]]), paste0(work_path, "/data_sci/obj_aligned_pd.rds"))

df = data.frame(obj[[]])
df_x = readRDS(paste0(work_path, "/data_sci/obj_processed_pd.rds"))
df$celltype_x = as.vector(df_x$celltype)
saveRDS(df, paste0(work_path, "/data_sci/obj_aligned_pd.rds"))

save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup"
fig = plot_ly(df, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~RNA_snn_res.1)
saveWidget(fig, paste0(save_path, "/sam_tape/gastruloid_sci_aligned_RNA_snn_res.1.html"), selfcontained = FALSE, libdir = "tmp")
fig = plot_ly(df, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~RNA_snn_res.2)
saveWidget(fig, paste0(save_path, "/sam_tape/gastruloid_sci_aligned_RNA_snn_res.2.html"), selfcontained = FALSE, libdir = "tmp")
fig = plot_ly(df, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~celltype_x)
saveWidget(fig, paste0(save_path, "/sam_tape/gastruloid_sci_aligned_celltype_x.html"), selfcontained = FALSE, libdir = "tmp")
fig = plot_ly(df, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~plate)
saveWidget(fig, paste0(save_path, "/sam_tape/gastruloid_sci_aligned_plate.html"), selfcontained = FALSE, libdir = "tmp")

anno = rep(NA, nrow(df))
anno[df$RNA_snn_res.1 %in% c(1,2,9,21)] = "Hindbrain"
anno[df$RNA_snn_res.1 %in% c(0)] = "Spinal cord"
anno[df$RNA_snn_res.1 %in% c(3,4,8)] = "Somites"
anno[df$RNA_snn_res.1 %in% c(6)] = "NMPs"
anno[df$RNA_snn_res.1 %in% c(5)] = "Epi/PGC"
anno[df$RNA_snn_res.1 %in% c(7,10,13)] = "Epithelial cells (not sure)"
anno[df$RNA_snn_res.1 %in% c(12,17)] = "Early neurons"
anno[df$RNA_snn_res.1 %in% c(11)] = "Floor plate"
anno[df$RNA_snn_res.1 %in% c(16)] = "Mesodermal progenitors (Tbx6+)"
anno[df$RNA_snn_res.1 %in% c(14)] = "Notochord"
anno[df$RNA_snn_res.1 %in% c(15)] = "Parietal endoderm"
anno[df$RNA_snn_res.1 %in% c(19)] = "Cardiac mesoderm"
anno[df$RNA_snn_res.1 %in% c(18)] = "Definitive endoderm"
anno[df$RNA_snn_res.1 %in% c(20)] = "Endothelial cells"
anno[df$RNA_snn_res.2 %in% c(30)] = "Mesodermal cells (Eomes+)"
df$celltype = as.vector(anno)

saveRDS(df, paste0(work_path, "/data_sci/obj_aligned_pd.rds"))

fig = plot_ly(df, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~celltype)
saveWidget(fig, paste0(save_path, "/sam_tape/gastruloid_sci_aligned_celltype.html"), selfcontained = FALSE, libdir = "tmp")





