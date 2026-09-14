
#####################
### can we plot 26 clonal clusters
### 1) how many monoclonal gastruloids in each
### 2) cell type distributions across clusters from monoclonal gastruloids
### 3) cell type distributions across clusters from polyclonal gastruloids

work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape"
source("~/work/scripts/utils.R")

pd = readRDS(paste0(work_path, "/revision/monoclonal_polyclonal_RNAseq/obj_batchelor_pd.rds"))

anno = read.csv(paste0(work_path, "/revision/debris_seq/well_anno.csv"))
anno$well_cluster = paste0("cluster_", anno$well_cluster)

assign_m = readRDS(paste0(work_path, "/revision/monoclonal_RNAseq/cell_well_cluster_assignment.rds"))

assign_p = readRDS(paste0(work_path, "/revision/polyclonal_RNAseq/cell_well_cluster_assignment.rds"))
assign_p = assign_p[assign_p$ratio >=3,]

########################
### focus on 26 clusters

well_cluster_list_final = paste0("cluster_", c(8,16,18,20,
                                               3,4,5,21,22,
                                               9,14,15,19,
                                               10,11,17,2,23,24,25,26,27,28,29,6,7))

color_map = c(rep("#52c05a", 4), rep("#528ac0", 5), rep("#dc87ba", 4), rep("grey50", 13))
names(color_map) = well_cluster_list_final

### 93 mono-gastruloids
anno = anno[anno$well_cluster %in% well_cluster_list_final,] %>% group_by(well_cluster) %>% tally()

### 45861 cells
assign_m = assign_m %>% filter(well_cluster %in% well_cluster_list_final) %>%
    select(cell_id, well_cluster) %>%
    left_join(pd[,c("cell_id","celltype")], by = "cell_id")

### 3899 cells
assign_p = assign_p %>% filter(well_cluster %in% well_cluster_list_final) %>%
    select(cell_id, well_cluster) %>%
    left_join(pd[,c("cell_id","celltype")], by = "cell_id")


anno$well_cluster = factor(anno$well_cluster, levels = well_cluster_list_final)
p1 = anno %>% ggplot(aes(x=well_cluster, y=n, fill=well_cluster)) +
    geom_bar(stat="identity") +
    labs(x="", y="") +
    theme_classic(base_size = 12) +
    theme(legend.position="none") +
    scale_y_continuous(breaks = seq(0, 14, by = 2), limits = c(0, 14)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values=color_map) +
    theme(axis.text.x = element_text(color="black", angle = 90, hjust = 1), axis.text.y = element_text(color="black"))


assign_m$well_cluster = factor(assign_m$well_cluster, levels = well_cluster_list_final)
assign_m$celltype = factor(assign_m$celltype, levels = names(gastruloid_celltype_color_code)[names(gastruloid_celltype_color_code) %in% assign_m$celltype])

df = assign_m %>%
    group_by(well_cluster, celltype) %>%
    tally() %>%
    dplyr::rename(cell_num = n) %>%
    ungroup() %>%
    tidyr::complete(well_cluster, celltype, fill = list(cell_num = 0))

df_sub = assign_m %>%
    group_by(well_cluster) %>%
    tally() %>%
    dplyr::rename(cell_num_total = n)
df = df %>%
    left_join(df_sub, by = "well_cluster") %>%
    mutate(percentage = cell_num/cell_num_total) %>%
    mutate(cell_num_pred_log2 = log2(cell_num_total + 1) * percentage)

df$well_cluster = factor(df$well_cluster, levels = well_cluster_list_final)
df$celltype = factor(df$celltype, levels = names(gastruloid_celltype_color_code)[names(gastruloid_celltype_color_code) %in% df$celltype])

p2 = ggplot() +
    geom_bar(data = df, aes(x = well_cluster, y = cell_num, group = celltype, fill = celltype), stat="identity", width = 0.6) +
    labs(x = "", y = "") +
    theme_classic(base_size = 12) +
    scale_fill_manual(values = gastruloid_celltype_color_code) +
    theme(legend.position="none") +
#    scale_y_continuous(breaks = seq(0, 14, by = 2), limits = c(0, 14)) +
    theme(axis.text.x = element_text(color="black", angle = 90, hjust = 1), axis.text.y = element_text(color="black"))


assign_p$well_cluster = factor(assign_p$well_cluster, levels = well_cluster_list_final)
assign_p$celltype = factor(assign_p$celltype, levels = names(gastruloid_celltype_color_code)[names(gastruloid_celltype_color_code) %in% assign_p$celltype])

df = assign_p %>%
    group_by(well_cluster, celltype) %>%
    tally() %>%
    dplyr::rename(cell_num = n) %>%
    ungroup() %>%
    tidyr::complete(well_cluster, celltype, fill = list(cell_num = 0))

df_sub = assign_p %>%
    group_by(well_cluster) %>%
    tally() %>%
    dplyr::rename(cell_num_total = n)
df = df %>%
    left_join(df_sub, by = "well_cluster") %>%
    mutate(percentage = cell_num/cell_num_total) %>%
    mutate(cell_num_pred_log2 = log2(cell_num_total + 1) * percentage)

df$well_cluster = factor(df$well_cluster, levels = well_cluster_list_final)
df$celltype = factor(df$celltype, levels = names(gastruloid_celltype_color_code)[names(gastruloid_celltype_color_code) %in% df$celltype])

p3 = ggplot() +
    geom_bar(data = df, aes(x = well_cluster, y = cell_num, group = celltype, fill = celltype), stat="identity", width = 0.6) +
    labs(x = "", y = "") +
    theme_classic(base_size = 12) +
    scale_fill_manual(values = gastruloid_celltype_color_code) +
    theme(legend.position="none") +
#    scale_y_continuous(breaks = seq(0, 14, by = 2), limits = c(0, 14)) +
    theme(axis.text.x = element_text(color="black", angle = 90, hjust = 1), axis.text.y = element_text(color="black"))

ggsave("~/share/cluster_celltype_dis.pdf", p1 / p2 / p3, width = 10, height = 10)



