
#########################################################################################################
### Different from monoclonal gastruloids, our polyclonal gastruloids have been separated by 10x directly

work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape"
data_path = "/net/shendure/vol10/projects/prime_editing/nobackup/data/2025_11_21_revision"
source("~/work/scripts/utils.R")

pd = readRDS(paste0(work_path, "/revision/monoclonal_polyclonal_RNAseq/obj_batchelor_pd.rds"))
table(pd$dataset)

pd = pd[pd$dataset == "polyclonal_RNAseq",]

sample_list= data.frame(
lane = c("Poly_L01_GRCm39","Poly_L01_GRCm39","Poly_L01_GRCm39","Poly_L01_GRCm39","Poly_L01_GRCm39","Poly_L01_GRCm39","Poly_L01_GRCm39","Poly_L01_GRCm39","Poly_L02_GRCm39","Poly_L02_GRCm39","Poly_L02_GRCm39","Poly_L02_GRCm39","Poly_L02_GRCm39","Poly_L02_GRCm39","Poly_L02_GRCm39","Poly_L02_GRCm39","Poly_L03_GRCm39","Poly_L03_GRCm39","Poly_L03_GRCm39","Poly_L03_GRCm39","Poly_L03_GRCm39","Poly_L03_GRCm39","Poly_L03_GRCm39","Poly_L03_GRCm39","Poly_L04_GRCm39","Poly_L04_GRCm39","Poly_L04_GRCm39","Poly_L04_GRCm39","Poly_L04_GRCm39","Poly_L04_GRCm39","Poly_L04_GRCm39","Poly_L04_GRCm39"),
sample = c("PL01_A01","PL01_A02","PL01_A04","PL01_A05","PL01_A06","PL01_A07","PL01_A10","PL01_A11","PL02_B01","PL02_B02","PL02_B05","PL02_B06","PL02_B08","PL02_B10","PL02_B11","PL02_B12","PL03_C02","PL03_C03","PL03_C04","PL03_C06","PL03_C08","PL03_C09","PL03_C11","PL03_C12","PL04_D01","PL04_D02","PL04_D04","PL04_D06","PL04_D07","PL04_D08","PL04_D09","PL04_D11"),
folder = c("Poly_A01_01","Poly_A02_02","Poly_A04_04","Poly_A05_05","Poly_A06_06","Poly_A07_07","Poly_A10_10","Poly_A11_11","Poly_B01_13","Poly_B02_14","Poly_B05_17","Poly_B06_18","Poly_B08_20","Poly_B10_22","Poly_B11_23","Poly_B12_24","Poly_C02_26","Poly_C03_27","Poly_C04_28","Poly_C06_30","Poly_C08_32","Poly_C09_33","Poly_C11_35","Poly_C12_36","Poly_D01_37","Poly_D02_38","Poly_D04_40","Poly_D06_42","Poly_D07_43","Poly_D08_44","Poly_D09_45","Poly_D11_47")
)

dat_list = list()
for(i in 1:nrow(sample_list)){
    print(paste0(i, "/", nrow(sample_list)))
    file_path <- paste0(data_path, "/", sample_list$lane[i], "/outs/per_sample_outs/", 
                    sample_list$folder[i], "/count/sample_filtered_barcodes.csv")
    if (file.exists(file_path) && file.info(file_path)$size > 0) {
      tmp <- read.csv(file_path, header = FALSE)
      dat_list[[i]] = data.frame(cell_id = paste0(sample_list$lane[i], "_", as.vector(tmp$V2)), 
            lane = sample_list$lane[i], sample = sample_list$sample[i], folder = sample_list$folder[i])
    }
}

dat = do.call(rbind, dat_list)

pd = pd %>% left_join(dat[,c("cell_id", "sample")], by = "cell_id") %>% as.data.frame()

saveRDS(pd, paste0(work_path, "/revision/polyclonal_RNAseq/pd_assigned_samples.rds"))




#######################################################
### comparing cell-type diversity between mono vs. poly

work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape"
data_path = "/net/shendure/vol10/projects/prime_editing/nobackup/data/2025_11_21_revision"
source("~/work/scripts/utils.R")

pd = readRDS(paste0(work_path, "/revision/monoclonal_polyclonal_RNAseq/obj_batchelor_pd.rds"))
table(pd$dataset)

df = pd %>% group_by(celltype, dataset) %>% tally() %>%
    left_join(pd %>% group_by(dataset) %>% tally() %>% rename(total_n = n), by = "dataset") %>%
    mutate(frac = 100 * n/total_n)

df$celltype = factor(df$celltype, levels = names(gastruloid_celltype_color_code)[names(gastruloid_celltype_color_code) %in% df$celltype])

p = ggplot(df, aes(x=celltype, y=frac, fill=celltype)) + 
    geom_bar(stat='identity') +
    scale_fill_manual(values=gastruloid_celltype_color_code) +
    facet_grid(.~dataset) + 
    coord_flip() +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
ggsave("~/share/celltype_composition_datasets.pdf", p, height = 5, width = 5)


#######################################################
### plotting individual polyclonal gastruloids in UMAPs

pd_all = readRDS(paste0(work_path, "/revision/monoclonal_polyclonal_RNAseq/obj_batchelor_pd.rds"))

pd = readRDS(paste0(work_path, "/revision/polyclonal_RNAseq/pd_assigned_samples.rds"))
gastruloid_num = pd[!is.na(pd$sample),] %>% group_by(sample) %>% tally() %>% arrange(-n) %>% as.data.frame()

### at least 200 cells
gastruloid_list = gastruloid_num %>% filter(n >= 200) %>% pull(sample)

plot_list = list()
for (i in gastruloid_list) {
    print(i)
    pd_x = pd %>% filter(sample == i & !is.na(sample))
    plot_list[[i]] = ggplot() +
        geom_point(data = pd_all, aes(x = UMAP_1, y = UMAP_2), color = "grey80",size = 0.1) +
        geom_point(data = pd_x, aes(x = UMAP_1, y = UMAP_2, color = "red"), size = 0.15) +
        theme_void() +
        labs(title = paste0(i, ", n = ", nrow(pd_x), " cells")) +
        theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 9))
}
combined_plot <- wrap_plots(plot_list, nrow = 3)

ggsave(paste0("~/share/polyclonal_gastruloids_umap.png"), combined_plot, width = 28, height = 12, dpi = 300)









##########################################################################
### Assigning cells from scRNA-seq of polyclonal gastruloids to debris-seq

work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape"
data_path = "/net/shendure/vol10/projects/prime_editing/nobackup/data/2025_11_21_revision"
source("~/work/scripts/utils.R")

pd = readRDS(paste0(work_path, "/revision/monoclonal_polyclonal_RNAseq/obj_batchelor_pd.rds"))
table(pd$dataset)
# monoclonal_RNAseq polyclonal_RNAseq
#             56436             26929

well_barcode = read.csv(paste0(work_path, "/revision/debris_seq/well_barcode_set.csv"), header=T)
well_barcode$well_cluster = paste0("cluster_", well_barcode$well_cluster)

well_barcode_content = unique(well_barcode[,c("well_id", "well_cluster")]) 
well_barcode_content_num = well_barcode_content %>% group_by(well_cluster) %>% tally()
well_cluster_list = paste0("cluster_", 1:nrow(well_barcode_content_num))

run_list = c("Poly_L01_GRCm39","Poly_L01_GRCm39","Poly_L01_GRCm39","Poly_L01_GRCm39","Poly_L01_GRCm39","Poly_L01_GRCm39","Poly_L01_GRCm39","Poly_L01_GRCm39","Poly_L02_GRCm39","Poly_L02_GRCm39","Poly_L02_GRCm39","Poly_L02_GRCm39","Poly_L02_GRCm39","Poly_L02_GRCm39","Poly_L02_GRCm39","Poly_L02_GRCm39","Poly_L03_GRCm39","Poly_L03_GRCm39","Poly_L03_GRCm39","Poly_L03_GRCm39","Poly_L03_GRCm39","Poly_L03_GRCm39","Poly_L03_GRCm39","Poly_L03_GRCm39","Poly_L04_GRCm39","Poly_L04_GRCm39","Poly_L04_GRCm39","Poly_L04_GRCm39","Poly_L04_GRCm39","Poly_L04_GRCm39","Poly_L04_GRCm39","Poly_L04_GRCm39")

run_short_list = c("PL01_A01","PL01_A02","PL01_A04","PL01_A05","PL01_A06","PL01_A07","PL01_A10","PL01_A11","PL02_B01","PL02_B02","PL02_B05","PL02_B06","PL02_B08","PL02_B10","PL02_B11","PL02_B12","PL03_C02","PL03_C03","PL03_C04","PL03_C06","PL03_C08","PL03_C09","PL03_C11","PL03_C12","PL04_D01","PL04_D02","PL04_D04","PL04_D06","PL04_D07","PL04_D08","PL04_D09","PL04_D11")

dat_all = NULL
for(kk in 1:length(run_list)){
    print(paste0(kk, "/", length(run_list)))
    dat_i = read.csv(paste0(data_path, "/data/Step1_", run_short_list[kk], "_CellByTapeV3_t2.csv"))[,c(1:4)]
    if(nrow(dat_i) > 0){
        dat_i$cell_id = paste0(run_list[kk], "_", dat_i$cBC)
        dat_i$barcode = dat_i$mBC
        dat_i$sample = run_short_list[kk]
        dat_i = dat_i[dat_i$cell_id %in% pd$cell_id & dat_i$barcode %in% well_barcode$barcode,]
        dat_all = rbind(dat_all, dat_i[,c("cell_id", "barcode", "n_UMI", "sample")])
    }
}

well_barcode = well_barcode[well_barcode$barcode %in% dat_all$barcode,]

### 11451/26929 = 43% cells at least one barcode recoved

##########################################
### Method-2: also consider the read count

dat_1 = well_barcode[,c("well_id", "barcode", "read_count")] %>% dcast(barcode ~ well_id, fill = 0)
rownames(dat_1) = dat_1[,1]; dat_1 = as.matrix(dat_1[,-1])
dat_1 = t(t(dat_1)/apply(dat_1, 2, sum))

dat_2 = dat_all[,c("cell_id", "barcode", "n_UMI")] %>% dcast(barcode ~ cell_id, fill = 0)
rownames(dat_2) = dat_2[,1]; dat_2 = as.matrix(dat_2[,-1])
dat_2 = t(t(dat_2)/apply(dat_2, 2, sum))

mat = t(dat_2) %*% dat_1

mat_regroup = NULL
for(i in well_cluster_list){
  print(i)
  well_include = well_barcode_content %>% filter(well_cluster == i) %>% pull(well_id)
  mat_regroup = cbind(mat_regroup, apply(mat[,colnames(mat) %in% well_include, drop = F], 1, mean))
}

rownames(mat_regroup) = rownames(mat)
colnames(mat_regroup) = well_cluster_list
mat = mat_regroup

max_idx = max.col(mat, ties.method = "first")
max_val = mat[cbind(seq_len(nrow(mat)), max_idx)]
mat2 = mat
mat2[cbind(seq_len(nrow(mat)), max_idx)] = -Inf
second_idx = max.col(mat2, ties.method = "first")
second_val = mat2[cbind(seq_len(nrow(mat)), second_idx)]

df = data.frame(cell_id = rownames(mat),
    max_idx = max_idx, max_val = max_val, second_idx = second_idx, second_val = second_val)
df$ratio = df$max_val/df$second_val
df$well_cluster = paste0("cluster_", df$max_idx)

saveRDS(df, paste0(work_path, "/revision/polyclonal_RNAseq/cell_well_cluster_assignment.rds"))

df_sub = df[df$ratio >=3,]
### n = 3,928 cells assigned; 3928/26929 = 15% of cells were assigned

well_cluster_list_final = df_sub %>% group_by(well_cluster) %>% tally() %>% 
 filter(n >= 200) %>% arrange(-n) %>% pull(well_cluster)
### only 1 cluster; 81% of cells were assigned to cluster_29

### if we narrow down to the 26 clusters as monoclonal gastruloids assignment
well_cluster_list_final = paste0("cluster_", c(8,16,18,20,
                                               3,4,5,21,22,
                                               9,14,15,19,
                                               10,11,17,2,23,24,25,26,27,28,29,6,7))
df_sub = df_sub[df_sub$well_cluster %in% well_cluster_list_final,]
sum(df_sub$well_cluster == "cluster_29")/nrow(df_sub) ### 82%
print(nrow(df_sub)/26929) ### 14%

pd = readRDS(paste0(work_path, "/revision/polyclonal_RNAseq/pd_assigned_samples.rds"))
gastruloid_num = pd[!is.na(pd$sample),] %>% group_by(sample) %>% tally() %>% arrange(-n) %>% as.data.frame()
gastruloid_list = gastruloid_num %>% filter(n >= 200) %>% pull(sample)
pd = pd[!is.na(pd$sample) & pd$sample %in% gastruloid_list,]
### 10708/26929 = 40%
### cell-to-gastruloid assignment by CMO demuxing, which is 40% of profiled cells assigned to 20 polyclonal gastruloids (at least 200 cells)

x = pd %>% select(cell_id, sample) %>% 
    left_join(df_sub %>% select(cell_id, well_cluster), by = 'cell_id')

x$well_cluster[is.na(x$well_cluster)] = "not_assigned"
x$well_cluster[!x$well_cluster %in% c("not_assigned", "cluster_29")] = "other_clusters"

x_p = x %>%
    group_by(sample, well_cluster) %>%
    tally() %>%
    dplyr::rename(cell_num = n) %>%
    ungroup() %>%
    tidyr::complete(sample, well_cluster, fill = list(cell_num = 0))

x_p$well_cluster = factor(x_p$well_cluster, levels = rev(c("not_assigned", "other_clusters", "cluster_29")))

p = ggplot() +
    geom_bar(data = x_p, aes(x = sample, y = cell_num, group = well_cluster, fill = well_cluster), stat="identity", width = 0.6) +
    labs(x = "", y = "") +
    theme_classic(base_size = 12) +
    scale_fill_manual(values = c("cluster_29" = "blue", "other_clusters" = "red", "not_assigned" = "grey80")) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black", angle = 45, hjust = 1), axis.text.y = element_text(color="black"))

ggsave("~/share/poly_cluster_29.pdf", p, width = 10, height = 4)


