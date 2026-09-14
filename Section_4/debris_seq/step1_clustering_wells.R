
################################################################################
### We started with 50 cells and allowed them to grow into larger populations. 
### Some of these cells were used for co-assays, some for poly-gastruloid formation, 
### and others for mono-gastruloid formation. 

### From the debris-seq data (mono-gastruloids), we obtained 107 wells, 
### which should correspond to the 50 original clusters.

work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape"
source("~/work/scripts/utils.R")

dat = read.csv(paste0(work_path, "/revision/debris_seq/Step1File_DebrisTAPE_v251118.csv"), header=F)
colnames(dat) = c("well", "barcode", "read_count")
dat$log2_read_count = log2(dat$read_count)

p1 = ggplot() + 
    geom_histogram(data = dat, aes(log2_read_count), bins = 100) + 
    plot_tmp

mat = dat %>% select(well, barcode, read_count) %>% dcast(barcode ~ well, fill = 0)
rownames(mat) = mat[,1]
mat = as.matrix(mat[,-1]) ### 1,014 barcodes x 107 wells, w/o filtering

read_count_per_barcode = apply(mat, 1, sum)

p2 = ggplot() + 
    geom_histogram(data = data.frame(log2_read_count_per_barcode = log2(read_count_per_barcode)), aes(log2_read_count_per_barcode), bins = 50) + 
    geom_vline(xintercept = log2(1000), color = "red") +
    plot_tmp

mat_filter = mat[read_count_per_barcode >= 1000,] ### 948 barcodes x 107 wells, after filtering low-count barcodes

read_count_per_well = apply(mat_filter, 2, sum)

p3 = ggplot() + 
    geom_histogram(data = data.frame(log2_read_count_per_well = log2(read_count_per_well)), aes(log2_read_count_per_well), bins = 25) + 
    geom_vline(xintercept = log2(1000), color = "red") +
    plot_tmp

mat_filter = mat_filter[,read_count_per_well >= 1000] ### 948 barcodes x 106 wells, after filtering low-count wells

ggsave(paste0(work_path, "/revision/debris_seq/QC_debris_seq_hist.pdf"), p1 + p2 + p3, height = 4, width = 15)

### normalization by the total counts per well
mat_filter_norm = t(t(mat_filter)/apply(mat_filter, 2, sum))
cor_mat = cor(mat_filter_norm)

pdf(paste0(work_path, "/revision/debris_seq/well_well_corr.pdf"))
p = heatmap.2(as.matrix(cor_mat), 
          col=Colors, 
          scale="none", 
          Rowv = T, 
          Colv = T, 
          key=T, 
          density.info="none", 
          trace="none", 
          cexRow = 0.3, 
          cexCol = 0.3,
          margins = c(5,5))
dev.off()



#####################################################################################
### As shown above, we observed some "doublets" wells, which are needed to be removed

exclude_well = c("Sample_P2_Well_C5", "Sample_P1_Well_C8", 
                 "Sample_P3_Well_A10", "Sample_P1_Well_A6",
                 "Sample_P1_Well_B10", "Sample_P1_Well_C10")

cor_mat_filter = cor_mat[!rownames(cor_mat) %in% exclude_well,
                         !rownames(cor_mat) %in% exclude_well]

pdf(paste0(work_path, "/revision/debris_seq/well_well_corr_filter.pdf"))
p = heatmap.2(as.matrix(cor_mat_filter), 
              col=Colors, 
              scale="none", 
              Rowv = T, 
              Colv = T, 
              key=T, 
              density.info="none", 
              trace="none", 
              cexRow = 0.3, 
              cexCol = 0.3,
              margins = c(5,5))
dev.off()

anno = data.frame(well_id = rev(rownames(cor_mat_filter)[p$rowInd]))
anno$well_cluster = cutree(as.hclust(p$rowDendrogram), k = 30)[as.vector(anno$well_id)]
anno$well_cluster[anno$well_cluster == 27] = 5

cluster_counts = unique(anno$well_cluster)
cluster_map = setNames(seq_along(cluster_counts), cluster_counts)
anno$well_cluster = cluster_map[as.character(anno$well_cluster)]

write.csv(anno, paste0(work_path, "/revision/debris_seq/well_anno.csv"), row.names = F)

anno = anno[match(rownames(cor_mat_filter), anno$well_id), ]
clusters = anno$well_cluster
cluster_levels = sort(unique(anno$well_cluster))
cluster_colors = setNames(
    colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(cluster_levels)),
    cluster_levels)
row_side_colors = cluster_colors[as.character(anno$well_cluster)]

pdf(paste0(work_path, "/revision/debris_seq/well_well_corr_filter.pdf"))
hm = heatmap.2(
    as.matrix(cor_mat_filter),
    col = Colors,
    scale = "none",
    Rowv = TRUE,
    Colv = TRUE,
    RowSideColors = row_side_colors,
    ColSideColors = row_side_colors,
    key = TRUE,
    density.info = "none",
    trace = "none",
    cexRow = 0.3,
    cexCol = 0.3,
    margins = c(6,5)
)
dev.off()

###############################################
### Output Barcodes for individual well cluster

df_mat_filter = melt(mat_filter)
colnames(df_mat_filter) = c("barcode", "well_id", "read_count")
df_mat_filter = df_mat_filter %>% filter(read_count > 0) %>% 
    left_join(anno, by = "well_id") %>% filter(!is.na(well_cluster)) %>% as.data.frame()
write.csv(df_mat_filter, paste0(work_path, "/revision/debris_seq/well_barcode_set.csv"), row.names = F)

df_mat_filter_num = df_mat_filter %>% select(well_cluster, barcode) %>% unique() %>%
    group_by(well_cluster) %>% tally()

### mean = 57, median = 50, ranging from 10 to 164






