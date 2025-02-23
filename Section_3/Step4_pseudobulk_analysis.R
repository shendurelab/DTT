
##################################
### Performing pseudobulk analysis
### Chengxiang Qiu
### Feb-20, 2025

source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

work_path = "Your_work_path"

###########################################################
### Step-1: aggregating cells to create pseudobulk profiles

### read mGASv5 dataset
count = readRDS(paste0(work_path, "/data_mGASv5/count.rds"))
pd_v5 = readRDS(paste0(work_path, "/data_mGASv5/obj_sctransform_pd.rds"))
pd_v5$Cell = paste0(unlist(lapply(as.vector(pd_v5$cell_id), function(x) strsplit(x,"[_]")[[1]][3])),
                    gsub("mGASv5_", "", pd_v5$experiment_id))

### my cell-to-gastruloid assignments, w/o filtering gastruloids with <50 cells assigned (n = 88 wells)
cell_assign_group = readRDS(paste0(work_path, "/data_mGASv5/tree_sanjay/assignment/assignment_well.rds"))
well_include = as.vector(cell_assign_group$well[cell_assign_group$cell_num >= 50])
cell_assign = readRDS(paste0(work_path, "/data_mGASv5/tree_sanjay/assignment/assignment_cell.rds"))
cell_assign = cell_assign[cell_assign$well %in% well_include,]
colnames(cell_assign) = c("Cell","node", "Well", "celltype")

pd_v5_x = pd_v5 %>% left_join(cell_assign[,c("Cell","Well")])
pd_v5$Well = as.vector(pd_v5_x$Well)

sample_list = unique(cell_assign$Well)
count_aggr = NULL
for(i in 1:length(sample_list)){
    sample_i = sample_list[i]
    print(sample_i)
    count_sub = count[,colnames(count) %in% rownames(pd_v5)[pd_v5$Well == sample_i & !is.na(pd_v5$Well)], drop=FALSE]
    count_aggr = cbind(count_aggr,
                       Matrix::rowSums(count_sub))
}
colnames(count_aggr) = sample_list

saveRDS(count_aggr, paste0(work_path, "/data_mGASv5/pca/pseudobulk_count.rds"))


#################################################
### Step-2: Performing PCA on mGASv5 dataset only

obj_sci = readRDS(paste0(work_path, "/data_sci/tape_barcode_3/pca/obj_pseudobulk.rds"))
count_aggr = readRDS(paste0(work_path, "/data_mGASv5/pca/pseudobulk_count.rds"))
count_aggr = count_aggr[rownames(count_aggr) %in% rownames(obj_sci),]

obj = CreateSeuratObject(count_aggr)
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2500)
genes_include = VariableFeatures(obj)

cds = new_cell_data_set(as.matrix(count_aggr))

set.seed(2016)
FM = monocle3:::normalize_expr_data(cds, 
                                    norm_method = "log",
                                    pseudo_count = 1)

FM = FM[genes_include,]

num_dim = 10
scaling = TRUE
set.seed(2016)
irlba_res = my_sparse_prcomp_irlba(Matrix::t(FM), 
                                   n = min(num_dim, min(dim(FM)) - 1), 
                                   center = scaling, 
                                   scale. = scaling)
preproc_res = irlba_res$x
row.names(preproc_res) = colnames(cds)

saveRDS(irlba_res, paste0(work_path, "/data_mGASv5/pca/PCA_pseudobulk.rds"))

prop_var_expl = irlba_res$sdev^2/sum(irlba_res$sdev^2)
print(prop_var_expl)
### 35.8%, 18.3%, 11.9%

#########################
### Plot the PC_1 vs PC_2
df_y = data.frame(sample_id = rownames(preproc_res),
                  PC_1 = preproc_res[,1],
                  PC_2 = preproc_res[,2],
                  PC_3 = preproc_res[,3])

well_list = readRDS(paste0(work_path, "/data_mGASv5/tree/sample_group.rds"))
df_y = df_y %>% left_join(well_list, by = "sample_id")

gas_size = read.table(paste0(work_path, "/data_mGASv5/tree/gastruloid_size.txt"))
colnames(gas_size) = c("sample_id", "size")
df_y = df_y %>% left_join(gas_size, by = "sample_id")

saveRDS(df_y, paste0(work_path, "/data_mGASv5/pca/pca_coor.rds"))

sample_group_color = c("#f26522", "#ec008c", "#2e3192", "#00aeef", "#00a651", "#a8b439", "#ed1c24")
names(sample_group_color) = paste0("group_", c(1:7))

library(gridExtra) 
plot_tmp = theme_classic(base_size = 12) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 

p1 = ggplot() + geom_point(data = df_y, aes(x = PC_1, y = PC_2, color = sample_group)) +
    labs(x = "PC_1 (35.8%)", y = "PC_2 (18.3%)") + 
    scale_color_manual(values=sample_group_color) + plot_tmp
p2 = ggplot() + geom_point(data = df_y, aes(x = PC_2, y = PC_3, color = sample_group)) +
    labs(x = "PC_2 (18.3%)", y = "PC_3 (11.9%)") + 
    scale_color_manual(values=sample_group_color) + plot_tmp
p3 = ggplot() + geom_point(data = df_y, aes(x = PC_3, y = PC_1, color = sample_group)) +
    labs(x = "PC_3 (11.9%)", y = "PC_1 (35.8%)") + 
    scale_color_manual(values=sample_group_color) + plot_tmp

ggsave("~/share/pseudobulk_pca_scatter.pdf", p1 + p2 + p3, width = 12, height = 4)


df = cell_assign %>%
    group_by(Well, celltype) %>% tally() %>%
    dcast(celltype~Well, fill = 0)
rownames(df) = as.vector(df[[1]])
df = df[,-1]
df_pct = 100*t(df)/apply(df,2,sum)
df_pct = df_pct[as.vector(df_y$sample_id),]

df_num = cell_assign %>% group_by(Well) %>% tally()
colnames(df_num) = c("sample_id", "n")
df_num$log2_num = log2(df_num$n + 1)
df_y = df_y %>% left_join(df_num[,c("sample_id","log2_num")], by = "sample_id")


df_y$tmp = df_pct[,colnames(df_pct) == "Epiblast"]
fit = cor.test(df_y$PC_1, df_y$tmp)
p1 = ggplot() + geom_point(data = df_y, aes(x = PC_1, y = tmp)) +
    geom_point(data = df_y[df_y$sample_group == "group_1",], aes(x = PC_1, y = tmp), color = "red", size = 3) +
    labs(x = "PC_1 (35.8%)", y = "% of PSC-like cells", title = paste0("Corr = ", round(fit$estimate, 2), ", P-val = ", sprintf("%.2e", fit$p.val))) + plot_tmp            
fit = cor.test(df_y$PC_2, df_y$tmp)
p2 = ggplot() + geom_point(data = df_y, aes(x = PC_2, y = tmp)) +
    labs(x = "PC_2 (18.3%)", y = "% of PSC-like cells", title = paste0("Corr = ", round(fit$estimate, 2), ", P-val = ", sprintf("%.2e", fit$p.val))) + plot_tmp            
fit = cor.test(df_y$PC_3, df_y$tmp)
p3 = ggplot() + geom_point(data = df_y, aes(x = PC_3, y = tmp)) +
    labs(x = "PC_3 (11.9%)", y = "% of PSC-like cells", title = paste0("Corr = ", round(fit$estimate, 2), ", P-val = ", sprintf("%.2e", fit$p.val))) + plot_tmp            
ggsave("~/share/pseudobulk_pca_corr_PSC.pdf", p1 + p2 + p3, width = 12, height = 4)

df_y$tmp = df_y$size
fit = cor.test(df_y$PC_1, df_y$tmp)
p1 = ggplot() + geom_point(data = df_y, aes(x = PC_1, y = tmp)) +
    labs(x = "PC_1 (35.8%)", y = "Gastruloid size", title = paste0("Corr = ", round(fit$estimate, 2), ", P-val = ", sprintf("%.2e", fit$p.val))) + plot_tmp            
fit = cor.test(df_y$PC_2, df_y$tmp)
p2 = ggplot() + geom_point(data = df_y, aes(x = PC_2, y = tmp)) +
    labs(x = "PC_2 (18.3%)", y = "Gastruloid size", title = paste0("Corr = ", round(fit$estimate, 2), ", P-val = ", sprintf("%.2e", fit$p.val))) + plot_tmp            
fit = cor.test(df_y$PC_3, df_y$tmp)
p3 = ggplot() + geom_point(data = df_y, aes(x = PC_3, y = tmp)) +
    labs(x = "PC_3 (11.9%)", y = "Gastruloid size", title = paste0("Corr = ", round(fit$estimate, 2), ", P-val = ", sprintf("%.2e", fit$p.val))) + plot_tmp            
ggsave("~/share/pseudobulk_pca_corr_size.pdf", p1 + p2 + p3, width = 12, height = 4)

df_y$tmp = df_y$log2_num
fit = cor.test(df_y$PC_1, df_y$tmp)
p1 = ggplot() + geom_point(data = df_y, aes(x = PC_1, y = tmp)) +
    labs(x = "PC_1 (35.8%)", y = "Log2 (cell # + 1)", title = paste0("Corr = ", round(fit$estimate, 2), ", P-val = ", sprintf("%.2e", fit$p.val))) + plot_tmp            
fit = cor.test(df_y$PC_2, df_y$tmp)
p2 = ggplot() + geom_point(data = df_y, aes(x = PC_2, y = tmp)) +
    labs(x = "PC_2 (18.3%)", y = "Log2 (cell # + 1)", title = paste0("Corr = ", round(fit$estimate, 2), ", P-val = ", sprintf("%.2e", fit$p.val))) + plot_tmp            
fit = cor.test(df_y$PC_3, df_y$tmp)
p3 = ggplot() + geom_point(data = df_y, aes(x = PC_3, y = tmp)) +
    labs(x = "PC_3 (11.9%)", y = "Log2 (cell # + 1)", title = paste0("Corr = ", round(fit$estimate, 2), ", P-val = ", sprintf("%.2e", fit$p.val))) + plot_tmp            
ggsave("~/share/pseudobulk_pca_corr_log2_cell_num.pdf", p1 + p2 + p3, width = 12, height = 4)



df_y$tmp = df_pct[,colnames(df_pct) == "Hindbrain"]
fit = cor.test(df_y$PC_1, df_y$tmp)
p1 = ggplot() + geom_point(data = df_y, aes(x = PC_1, y = tmp)) +
    labs(x = "PC_1 (35.8%)", y = "% of Hindbrain", title = paste0("Corr = ", round(fit$estimate, 2), ", P-val = ", sprintf("%.2e", fit$p.val))) + plot_tmp            
fit = cor.test(df_y$PC_2, df_y$tmp)
p2 = ggplot() + geom_point(data = df_y, aes(x = PC_2, y = tmp)) +
    labs(x = "PC_2 (18.3%)", y = "% of Hindbrain", title = paste0("Corr = ", round(fit$estimate, 2), ", P-val = ", sprintf("%.2e", fit$p.val))) + plot_tmp            
fit = cor.test(df_y$PC_3, df_y$tmp)
p3 = ggplot() + geom_point(data = df_y, aes(x = PC_3, y = tmp)) +
    geom_point(data = df_y[df_y$sample_group == "group_6",], aes(x = PC_3, y = tmp), color = "#FFA500", size = 3) +
    labs(x = "PC_3 (11.9%)", y = "% of Hindbrain", title = paste0("Corr = ", round(fit$estimate, 2), ", P-val = ", sprintf("%.2e", fit$p.val))) + plot_tmp            
ggsave("~/share/pseudobulk_pca_corr_Hindbrain.pdf", p1 + p2 + p3, width = 12, height = 4)

df_y$tmp = df_pct[,colnames(df_pct) == "Somites"]
fit = cor.test(df_y$PC_1, df_y$tmp)
p1 = ggplot() + geom_point(data = df_y, aes(x = PC_1, y = tmp)) +
    labs(x = "PC_1 (35.8%)", y = "% of Somites", title = paste0("Corr = ", round(fit$estimate, 2), ", P-val = ", sprintf("%.2e", fit$p.val))) + plot_tmp            
fit = cor.test(df_y$PC_2, df_y$tmp)
p2 = ggplot() + geom_point(data = df_y, aes(x = PC_2, y = tmp)) +
    labs(x = "PC_2 (18.3%)", y = "% of Somites", title = paste0("Corr = ", round(fit$estimate, 2), ", P-val = ", sprintf("%.2e", fit$p.val))) + plot_tmp            
fit = cor.test(df_y$PC_3, df_y$tmp)
p3 = ggplot() + geom_point(data = df_y, aes(x = PC_3, y = tmp)) +
    labs(x = "PC_3 (11.9%)", y = "% of Somites", title = paste0("Corr = ", round(fit$estimate, 2), ", P-val = ", sprintf("%.2e", fit$p.val))) + plot_tmp            
ggsave("~/share/pseudobulk_pca_corr_Somites.pdf", p1 + p2 + p3, width = 12, height = 4)


######################################################################################
### Step-3: For individual PC1-3, plotting its own value vs. its most neighbor's value

df = readRDS(paste0(work_path, "/data_mGASv5/pca/pca_coor.rds"))

### lineage distance in tree
tree = readRDS(paste0(work_path, "/data_mGASv5/tree/tree_well_tape_dominant_edits.rds"))

dat_i = castor::get_all_pairwise_distances(tree, 
                                           only_clades = 1:length(tree$tip.label), 
                                           as_edge_counts = TRUE)
rownames(dat_i) = colnames(dat_i) = 1:length(tree$tip.label)
dat_i = melt(as.matrix(dat_i))
dat_i$A = tree$tip.label[dat_i$Var1]
dat_i$B = tree$tip.label[dat_i$Var2]
dat_i$lineage_dist = dat_i$value
dat_i$Var1 = dat_i$Var2 = dat_i$value = NULL

dat_x = dat_i %>% filter(A != B, A %in% df$sample_id, B %in% df$sample_id) %>% 
    group_by(A) %>% slice_min(order_by = lineage_dist, n = 1)
### n = 112

sort_pair <- function(x, y) {
    apply(cbind(x, y), 1, function(row) paste(sort(row), collapse = "-"))
}
dat_x$pair = sort_pair(dat_x$A, dat_x$B)
dat_x_uniq = dat_x[!duplicated(dat_x$pair), ]
### n = 76


sample_group = readRDS(paste0(work_path, "/data_mGASv5/tree/sample_group.rds"))
sample_group_color = c("#f26522", "#ec008c", "#2e3192", "#00aeef", "#00a651", "#a8b439", "#ed1c24")
names(sample_group_color) = paste0("group_", c(1:7))


### Based on 112 pairs
dat = dat_x[,c("A","B")]
dat = dat %>% mutate(sample_id = A) %>% left_join(df[,c("sample_id","PC_1")], by = "sample_id") %>%
    left_join(sample_group, by = "sample_id") %>%
    rename(PC_A = PC_1) %>% mutate(sample_id = B) %>% left_join(df[,c("sample_id","PC_1")], by = "sample_id") %>%
    rename(PC_B = PC_1) 
fit = cor.test(dat$PC_A, dat$PC_B)
p1 = ggplot() + geom_point(data = dat, aes(x = PC_A, y = PC_B)) +
    geom_point(data = dat[dat$sample_group == "group_1",], aes(x = PC_A, y = PC_B), color = "red", size = 3) +
    labs(x = "PC_1 (Original)", y = "PC_1 (Nearest)", title = paste0("Corr = ", round(fit$estimate, 2), ", P-val = ", sprintf("%.2e", fit$p.val))) + plot_tmp            

dat = dat_x[,c("A","B")]
dat = dat %>% mutate(sample_id = A) %>% left_join(df[,c("sample_id","PC_2")], by = "sample_id") %>%
    left_join(sample_group, by = "sample_id") %>%
    rename(PC_A = PC_2) %>% mutate(sample_id = B) %>% left_join(df[,c("sample_id","PC_2")], by = "sample_id") %>%
    rename(PC_B = PC_2)
fit = cor.test(dat$PC_A, dat$PC_B)
p2 = ggplot() + geom_point(data = dat, aes(x = PC_A, y = PC_B)) +
    labs(x = "PC_2 (Original)", y = "PC_2 (Nearest)", title = paste0("Corr = ", round(fit$estimate, 2), ", P-val = ", sprintf("%.2e", fit$p.val))) + plot_tmp            

dat = dat_x[,c("A","B")]
dat = dat %>% mutate(sample_id = A) %>% left_join(df[,c("sample_id","PC_3")], by = "sample_id") %>%
    left_join(sample_group, by = "sample_id") %>%
    rename(PC_A = PC_3) %>% mutate(sample_id = B) %>% left_join(df[,c("sample_id","PC_3")], by = "sample_id") %>%
    rename(PC_B = PC_3)
fit = cor.test(dat$PC_A, dat$PC_B)
p3 = ggplot() + geom_point(data = dat, aes(x = PC_A, y = PC_B)) +
    geom_point(data = dat[dat$sample_group == "group_6",], aes(x = PC_A, y = PC_B), color = "#FFA500", size = 3) +
    labs(x = "PC_3 (Original)", y = "PC_3 (Nearest)", title = paste0("Corr = ", round(fit$estimate, 2), ", P-val = ", sprintf("%.2e", fit$p.val))) + plot_tmp            

ggsave("~/share/pseudobulk_pca_original_nearest.pdf", p1 + p2 + p3, width = 12, height = 4)

x = dat[dat$PC_A > 20 & dat$PC_B > 20,]
x = unique(c(x$A, x$B))
df_pct[x, colnames(df_pct) == "Hindbrain"]
P1-B5    P1-C8    P2-E1    P2-E3
79.08879 70.21792 38.27160 59.95717

tree = readRDS(paste0(work_path, "/tree_sanjay/tree.rds"))
p_cir = readRDS(paste0(work_path, "/tree_sanjay/tree_plot_cir.rds"))
select_node = c(2213, 2448, 2985, 2986) + 58283
### P1-B5, P2-E1, P1-C8, P2-E3
p_x <- p_cir +
    geom_label2(aes(subset = node %in% select_node, label = node), size = 1) +
    theme(legend.position="none")
ggsave(paste0("~/share/tree_mGASv5_cir_four_nodes.png"), p_x, dpi = 300)

P1-B5     P1-C8     P2-E1     P2-E3
"#A9A9A9" "#FFDAB9" "#F08080" "#EEE8AA"





######################################################################################################################
### Step-4: For individual PC1-3, can we compare the observed values categoried by major clades vs. 5,000 permutations

work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape"
source("~/work/scripts/utils.R")

df = readRDS(paste0(work_path, "/data_mGASv5/pca/pca_coor.rds"))

perm_time = 5000

df_perm = list()

for(perm_i in 1:perm_time){
    print(perm_i)
    df_i = df
    df_i$sample_group = df$sample_group[sample(1:nrow(df))]
    if(perm_i == 1){
        for(i in 1:7){
            df_perm[[i]] = df_i[df_i$sample_group == paste0("group_", i), c(2,3,4)]
        }
    } else {
        for(i in 1:7){
            df_perm[[i]] = rbind(df_perm[[i]], df_i[df_i$sample_group == paste0("group_", i), c(2,3,4)])
        }
    }
}

saveRDS(df_perm, paste0(work_path, "/data_mGASv5/pca/pca_clades_perm.rds"))

res = NULL
for(i in 1:7){
    df_x = as.matrix(df[df$sample_group == paste0("group_",i), c(2,3,4)])
    df_y = as.matrix(df_perm[[i]])
    
    res_i = NULL
    for(j in 1:3){
        fit = wilcox.test(df_x[,j], df_y[,j])
        res_i = c(res_i, fit$p.val)
    }
    res = rbind(res, res_i)
}

res = data.frame(res)
rownames(res) = paste0("clade_", 1:7)
colnames(res) = paste0("PC_", 1:3)

               PC_1      PC_2         PC_3
clade_1 0.007318307 0.8095746 0.4055825852
clade_2 0.406097533 0.6852646 0.0292585424
clade_3 0.777697466 0.8011754 0.4878554616
clade_4 0.100091206 0.8243533 0.3208682996
clade_5 0.926649697 0.7431665 0.2919047290
clade_6 0.002125825 0.3862184 0.0007909652
clade_7 0.148881120 0.9482230 0.3719812718





################################################################
### Step-5: cell type composition bias between gastruloid groups

sample_group = readRDS(paste0(work_path, "/data_mGASv5/tree/sample_group.rds"))

### my cell-to-gastruloid assignments
cell_assign_group = readRDS(paste0(work_path, "/data_mGASv5/tree_sanjay/assignment/assignment_well.rds"))
well_include = as.vector(cell_assign_group$well[cell_assign_group$cell_num >= 50])
cell_assign = readRDS(paste0(work_path, "/data_mGASv5/tree_sanjay/assignment/assignment_cell.rds"))
cell_assign = cell_assign[cell_assign$well %in% well_include,]
colnames(cell_assign) = c("Cell","node", "Well", "celltype")

celltype_num = cell_assign %>% group_by(celltype) %>% tally() %>% arrange(desc(n))

df = cell_assign %>% group_by(Well, celltype) %>% tally() %>% ungroup() %>%
    tidyr::complete(Well, celltype, fill = list(n = 0)) %>%
    left_join(cell_assign %>% group_by(Well) %>% tally() %>% rename(total_n = n), by = "Well") %>%
    mutate(frac = 100*n/total_n, sample_id = Well) %>%
    left_join(sample_group, by = "sample_id") %>%
    filter(celltype %in% as.vector(celltype_num$celltype)[1:8]) %>% as.data.frame()
df$celltype = factor(df$celltype, levels = as.vector(celltype_num$celltype)[1:8])
df$sample_group = factor(df$sample_group, levels = paste0("group_", c(1:7)))

sample_group_color = c("#f26522", "#ec008c", "#2e3192", "#00aeef", "#00a651", "#a8b439", "#ed1c24")
names(sample_group_color) = paste0("group_", c(1:7))

p = ggplot(data = df, aes(x = sample_group, y = frac, fill = sample_group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 0.8) +
    labs(x = "", y = "% of cells") + 
    scale_fill_manual(values = sample_group_color) +
    theme_classic(base_size = 10) +
    theme(legend.position = "none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color = "black"), 
          axis.text.y = element_text(color = "black")) +
    facet_wrap(~celltype, scales = "free_y", nrow = 4, ncol = 2) # Split panels by "celltype" and allow free y-axes

ggsave("~/share/data_mGASV5_celltype_frac_by_group.pdf", p, width = 4, height = 5)




######################################################################################
### If we permute the gastruloid labels, what are the cell type compositions changing?

df = df[,c("sample_id","celltype","frac","sample_group")]
df_x = df[,c("sample_id","celltype","frac")]

perm_time = 5000
df_perm = NULL
for(perm_i in 1:perm_time){
    print(perm_i)
    sample_group_i = sample_group
    sample_group_i$sample_group = sample_group$sample_group[sample(1:nrow(sample_group))]
    df_i = df_x %>%
        left_join(sample_group_i, by = "sample_id") %>%
        filter(celltype %in% as.vector(celltype_num$celltype)[1:8]) %>% as.data.frame()
    df_perm = rbind(df_perm, df_i)
}

df$group = "obs"
df_perm$group = "perm"

dat = rbind(df, df_perm)

dat$celltype = factor(dat$celltype, levels = as.vector(celltype_num$celltype)[1:8])
dat$sample_group = factor(dat$sample_group, levels = paste0("group_", c(1:7)))

saveRDS(dat, paste0(work_path, "/data_mGASv5/tree/clade_celltype_compositions_update.rds"))

sample_group_color = c("#f26522", "#ec008c", "#2e3192", "#00aeef", "#00a651", "#a8b439", "#ed1c24")
names(sample_group_color) = paste0("group_", c(1:7))

p = ggplot() +
    geom_boxplot(data = dat,
                 aes(x = sample_group, y = frac, fill = group), 
                 outlier.shape = NA) +
    geom_jitter(data = subset(dat, group == "obs"), 
                aes(x = as.numeric(sample_group) - 0.2, y = frac), # Shift points left
                color = "black", width = 0.1, size = 0.4) + # Reduce width slightly
    scale_fill_manual(values = c("obs" = "red", "perm" = "grey80")) +
    labs(x = "", y = "% of cells") + 
    theme_classic(base_size = 10) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(color = "black"), 
          axis.text.y = element_text(color = "black")) +
    facet_wrap(~celltype, scales = "free_y", nrow = 4, ncol = 2)
ggsave("~/share/data_mGASV5_celltype_frac_by_group_perm.pdf", p, width = 4, height = 6)


### wilcox-test
x_obs = dat %>% filter(celltype == "Epiblast", sample_group == "group_1", group == "obs") %>% pull(frac) 
x_pem = dat %>% filter(celltype == "Epiblast", sample_group == "group_1", group == "perm") %>% pull(frac) 
wilcox.test(x_obs, x_pem)
### p = 0.00839

x_obs = dat %>% filter(celltype == "Transitional cells", sample_group == "group_1", group == "obs") %>% pull(frac) 
x_pem = dat %>% filter(celltype == "Transitional cells", sample_group == "group_1", group == "perm") %>% pull(frac) 
wilcox.test(x_obs, x_pem)
### p = 0.03735

x_obs = dat %>% filter(celltype == "Hindbrain", sample_group == "group_6", group == "obs") %>% pull(frac) 
x_pem = dat %>% filter(celltype == "Hindbrain", sample_group == "group_6", group == "perm") %>% pull(frac) 
wilcox.test(x_obs, x_pem)
### p = 0.0002421












