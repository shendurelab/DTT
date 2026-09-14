
data_path=/net/shendure/vol10/projects/akeith/nobackup/analyses/SR_revisions/ES_multiome/count
work_path=/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape/revision/data_10X_multiome/mito
python "$work_path"/count_chrM.py "$data_path"/arc_L1/outs/atac_possorted_bam.bam "$work_path"/arc_L1_chrM.csv
python "$work_path"/count_chrM.py "$data_path"/arc_L2/outs/atac_possorted_bam.bam "$work_path"/arc_L2_chrM.csv
python "$work_path"/count_chrM.py "$data_path"/arc_L3/outs/atac_possorted_bam.bam "$work_path"/arc_L3_chrM.csv
python "$work_path"/count_chrM.py "$data_path"/arc_L4/outs/atac_possorted_bam.bam "$work_path"/arc_L4_chrM.csv



#############

work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape"
source("~/work/scripts/utils.R")

cell_assignment = readRDS(paste0(work_path, "/revision/data_10X_multiome/cell_well_cluster_assignment.rds"))
dat_1 = data.frame(well_cluster = paste0("cluster_", c(8,16,18,20)), category = "PSC")
dat_2 = data.frame(well_cluster = paste0("cluster_", c(3,4,5,21,22)), category = "Neuroectodermal")
dat_3 = data.frame(well_cluster = paste0("cluster_", c(9,14,15,19)), category = "Mesodermal")
dat_4 = data.frame(well_cluster = paste0("cluster_", c(10,11,17,2,23,24,25,26,27,28,29,6,7)), category = "Mixed")
pd = rbind(dat_1, dat_2, dat_3, dat_4)
cell_assignment = cell_assignment %>% left_join(pd, by = "well_cluster") %>% filter(!is.na(category))
rownames(cell_assignment) = cell_assignment$cell_id

mito = NULL
for(i in 1:4){
    mito_i = read.table(paste0(work_path, "/revision/data_10X_multiome/mito/arc_L", i, "_chrM.csv"), sep="\t", header=T)
    mito_i$cell_id = gsub("-1", paste0("-", i), mito_i$cell_id)
    mito = rbind(mito, mito_i)
}

df = cell_assignment %>% left_join(mito, by = "cell_id")
df$mito_pct = 100*df$chrM_reads/df$total_reads

df %>% group_by(category) %>% summarize(mean_mito_pct = mean(mito_pct), median_mito_pct = median(mito_pct))

res = NULL
category_list = c("PSC", "Neuroectodermal", "Mesodermal", "Mixed")
for(i in 1:3){
    for(j in c((i+1):4)){
        category_i = category_list[i]
        category_j = category_list[j]
        fit = wilcox.test(df$mito_pct[df$category == category_i], df$mito_pct[df$category == category_j])
        res = rbind(res, data.frame(category_1 = category_i, category_2 = category_j, pval = fit$p.val))
    }
}


clone_color_map = c("#52c05a", "#528ac0", "#dc87ba", "grey50")
names(clone_color_map) = c("PSC", "Neuroectodermal", "Mesodermal", "Mixed")

p = ggplot(df, aes(category, mito_pct, fill = category)) + 
    geom_boxplot(outlier.shape = NA) +
    labs(x="", y="Mito%") +
    theme_classic(base_size = 12) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values=clone_color_map) +
#    coord_cartesian(ylim = c(0, 4)) +
    theme(axis.text.x = element_text(color="black", angle = 45), axis.text.y = element_text(color="black")) 
ggsave("~/share/box_plot.pdf", p, width = 4, height = 4)

