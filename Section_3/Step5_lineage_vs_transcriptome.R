
#######################################################
### Comparing transcriptome distance vs. phylo distance
### Chengxiang Qiu
### Feb-20, 2025

source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

work_path = "Your_work_path"

########################################################################################
### Step-1: Elucidean distance in the PCA space between gastruloids vs. lineage distance

### my cell-to-gastruloid assignments
cell_assign_group = readRDS(paste0(work_path, "/data_mGASv5/tree_sanjay/assignment/assignment_well.rds"))
well_include = as.vector(cell_assign_group$well[cell_assign_group$cell_num >= 50])
cell_assign = readRDS(paste0(work_path, "/data_mGASv5/tree_sanjay/assignment/assignment_cell.rds"))
cell_assign = cell_assign[cell_assign$well %in% well_include,]
colnames(cell_assign) = c("Cell","node", "Well", "celltype")

### tree of trees
tree = readRDS(paste0(work_path, "/data_mGASv5/tree/tree_well_tape_dominant_edits.rds"))

### lineage distance
dat_i = castor::get_all_pairwise_distances(tree, 
                                           only_clades = 1:length(tree$tip.label), 
                                           as_edge_counts = TRUE)
rownames(dat_i) = colnames(dat_i) = 1:length(tree$tip.label)
dat_i = melt(as.matrix(dat_i))
dat_i$A = tree$tip.label[dat_i$Var1]
dat_i$B = tree$tip.label[dat_i$Var2]
dat_i$lineage_dist = dat_i$value
dat_i$Var1 = dat_i$Var2 = dat_i$value = NULL

### Elucidean distance on PCA 
### pca coordinates after projecting onto mGASv4 dataset
pca_coor = readRDS(paste0(work_path, "/data_mGASv5/pca/pca_coor.rds"))
irlba_res = readRDS(paste0(work_path, "/data_mGASv5/pca/PCA_pseudobulk.rds"))
preproc_res = irlba_res$x
row.names(preproc_res) = as.vector(pca_coor$sample_id)
pca_coor = as.matrix(preproc_res)

euclidean_distance <- function(x, y) {
    sqrt(sum((x - y)^2))
}

###
gastruloid_list = rownames(pca_coor)
pca_dist = NULL
for(i in 1:(length(gastruloid_list)-1)){
    for(j in (i+1):length(gastruloid_list)){
        pca_dist = rbind(pca_dist, data.frame(A = gastruloid_list[i],
                                              B = gastruloid_list[j],
                                              pca_dist = euclidean_distance(pca_coor[i,], pca_coor[j,])))
    }
}

df = dat_i %>% left_join(pca_dist, by = c("A","B")) %>% filter(!is.na(pca_dist))
df$lineage_dist = factor(df$lineage_dist, levels = 2:24)
### 3828 pairs

p = ggplot(df, aes(lineage_dist, pca_dist, fill = lineage_dist)) + geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2) +
    labs(x="", y="Euclidean distances in the PCA space") +
    theme_classic(base_size = 12) +
    theme(legend.position="none") +
    scale_fill_viridis(discrete=TRUE) +
    theme(axis.text.x = element_text(color="black", angle = 45), axis.text.y = element_text(color="black"))
ggsave("~/share/PCA_dist_vs_lineage_dist.pdf", p)

cor.test(df$pca_dist, as.numeric(df$lineage_dist), method = "spearman")
### rho = 0.05240494; p-value = 0.001181

sample_group = readRDS(paste0(work_path, "/data_mGASv5/tree/sample_group.rds"))
df_sub = df %>% filter(lineage_dist == 2) %>% arrange(pca_dist) %>%
    mutate(sample_id = A) %>% left_join(sample_group, by = "sample_id")
df_sub$sample_id = NULL

### plotting it using log2 scale on y-axis
df$log2_pca_dist = log2(df$pca_dist)

p = ggplot(df, aes(lineage_dist, log2_pca_dist, fill = lineage_dist)) + geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2) +
    labs(x="", y="Log2(Euclidean distances in the PCA space)") +
    theme_classic(base_size = 12) +
    theme(legend.position="none") +
    scale_fill_viridis(discrete=TRUE) +
    theme(axis.text.x = element_text(color="black", angle = 45), axis.text.y = element_text(color="black"))
ggsave("~/share/PCA_dist_log2_vs_lineage_dist.pdf", p)

df_1 = df


##############################
### Step-2, kNN based approach

### read mGASv5 dataset
obj = readRDS(paste0(work_path, "/data_mGASv5/obj_sctransform.rds"))
pca_coor = Embeddings(obj,reduction = "pca")
pd_v5 = readRDS(paste0(work_path, "/data_mGASv5/obj_sctransform_pd.rds"))
pd_v5$Cell = paste0(unlist(lapply(as.vector(pd_v5$cell_id), function(x) strsplit(x,"[_]")[[1]][3])),
                    gsub("mGASv5_", "", pd_v5$experiment_id))

### my cell-to-gastruloid assignments
cell_assign_group = readRDS(paste0(work_path, "/data_mGASv5/tree_sanjay/assignment/assignment_well.rds"))
well_include = as.vector(cell_assign_group$well[cell_assign_group$cell_num >= 50])
cell_assign = readRDS(paste0(work_path, "/data_mGASv5/tree_sanjay/assignment/assignment_cell.rds"))
cell_assign = cell_assign[cell_assign$well %in% well_include,]
colnames(cell_assign) = c("Cell","node", "Well", "celltype")

pd_v5_x = pd_v5 %>% left_join(cell_assign[,c("Cell","Well")])
pd_v5$Well = as.vector(pd_v5_x$Well)

### subset each Well to roughly similar cells
pd_v5_sub = pd_v5 %>% filter(!is.na(Well)) %>% as.data.frame()
rownames(pd_v5_sub) = as.vector(pd_v5_sub$cell_id)
pca_coor_sub = pca_coor[rownames(pd_v5_sub),]


### calculate kNNs
k.param = 15; nn.method = "rann"; nn.eps = 0; annoy.metric = "euclidean"
nn.ranked = Seurat:::NNHelper(
    data = pca_coor_sub,
    k = k.param + 1,
    method = nn.method,
    searchtype = "standard",
    eps = nn.eps,
    metric = annoy.metric)
nn.ranked = Indices(object = nn.ranked)
nn_matrix = nn.ranked[,-1]

resultA = NULL
for(i in 1:k.param){
    print(i)
    resultA = cbind(resultA, as.vector(pd_v5_sub$Well)[as.vector(nn_matrix[,i])])
}

dat_x = data.frame(A = rep(pd_v5_sub$Well, times = ncol(resultA)),
                   B = c(resultA))

dat = dat_x %>% group_by(A, B) %>% tally() %>%
    dcast(A~B, fill=0)
rownames(dat) = dat[,1]
dat = as.matrix(dat[,-1])

dat_norm = t(t(dat)/apply(dat, 2, sum))
dat_norm = dat_norm/apply(dat_norm, 1, sum)

### calculate average of A-B and B-A
dat_norm_mean = (dat_norm + t(dat_norm))/2

mean_pct = melt(dat_norm_mean)
colnames(mean_pct) = c("A", "B", "pct")

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
dat_i$A = gsub("[.]", "", dat_i$A)
dat_i$B = gsub("[.]", "", dat_i$B)

dat_x = dat_i %>% filter(A != B, A %in% rownames(dat_norm_mean), B %in% rownames(dat_norm_mean)) %>% 
    group_by(A) %>% slice_min(order_by = lineage_dist, n = 1)

sort_pair <- function(x, y) {
    apply(cbind(x, y), 1, function(row) paste(sort(row), collapse = "-"))
}
dat_x$pair = sort_pair(dat_x$A, dat_x$B)
dat_x_uniq = dat_x[!duplicated(dat_x$pair), ]

x = NULL
for(i in 1:nrow(dat_x_uniq)){
    x = c(x, mean_pct$pct[mean_pct$A == dat_x_uniq$A[i] & mean_pct$B == dat_x_uniq$B[i]])
}
dat_x_uniq$mean_pct = x*100
### n = 76


dat_y = dat_i %>% filter(A != B, A %in% rownames(dat_norm_mean), B %in% rownames(dat_norm_mean))
dat_y$pair = sort_pair(dat_y$A, dat_y$B)
dat_y_uniq = dat_y[!duplicated(dat_y$pair), ]

set.seed(1234)
dat_y_uniq_sample = dat_y_uniq[sample(1:nrow(dat_y_uniq), 1000),]

y = NULL
for(i in 1:nrow(dat_y_uniq_sample)){
    y = c(y, mean_pct$pct[mean_pct$A == dat_y_uniq_sample$A[i] & mean_pct$B == dat_y_uniq_sample$B[i]])
}
dat_y_uniq_sample$mean_pct = y*100


df = data.frame(mean_pct = c(dat_x_uniq$mean_pct, dat_y_uniq_sample$mean_pct),
                group = c(rep("Nearest neighbors", nrow(dat_x_uniq)), rep("Random", nrow(dat_y_uniq_sample))))
p = ggplot(df, aes(group, mean_pct, fill = group)) + geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 0.3) +
    labs(x="", y="Mean % of kNN cells from each other") +
    theme_classic(base_size = 12) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black", angle = 45), axis.text.y = element_text(color="black"))
ggsave("~/share/knn_pct_vs_lineage_dist_boxplot.pdf", p, width=2.5, height=5)

### p-val = 0.01097
t.test(dat_x_uniq$mean_pct, dat_y_uniq_sample$mean_pct)


### making boxplot, with lineage distance as categories.

y = NULL
for(i in 1:nrow(dat_y_uniq)){
    y = c(y, mean_pct$pct[mean_pct$A == dat_y_uniq$A[i] & mean_pct$B == dat_y_uniq$B[i]])
}
dat_y_uniq$mean_pct = y*100
dat_y_uniq$lineage_dist = factor(dat_y_uniq$lineage_dist, levels = 2:24)

p = ggplot(dat_y_uniq, aes(lineage_dist, mean_pct, fill = lineage_dist)) + geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2) +
    labs(x="", y="Mean % of kNN cells from each other") +
    theme_classic(base_size = 12) +
    theme(legend.position="none") +
    scale_fill_viridis(discrete=TRUE) +
    theme(axis.text.x = element_text(color="black", angle = 45), axis.text.y = element_text(color="black"))
ggsave("~/share/knn_pct_category_by_lineage_dist.pdf", p)

cor.test(dat_y_uniq$mean_pct, as.numeric(dat_y_uniq$lineage_dist), method = "spearman")
### rho = -0.09373124; p-value = 6.242e-09

sample_group = readRDS(paste0(work_path, "/data_mGASv5/tree/sample_group.rds"))
df_sub = dat_y_uniq %>% filter(lineage_dist == 2) %>% arrange(desc(mean_pct)) %>%
    mutate(sample_id = A) %>% left_join(sample_group, by = "sample_id")
df_sub$sample_id = NULL
df_sub$pair = NULL

### log2-scale
dat_y_uniq$log2_mean_pct = log2(dat_y_uniq$mean_pct + 1)
p = ggplot(dat_y_uniq, aes(lineage_dist, log2_mean_pct, fill = lineage_dist)) + geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2) +
    labs(x="", y="Log2(Mean % of kNN cells from each other + 1)") +
    theme_classic(base_size = 12) +
    theme(legend.position="none") +
    scale_fill_viridis(discrete=TRUE) +
    theme(axis.text.x = element_text(color="black", angle = 45), axis.text.y = element_text(color="black"))
ggsave("~/share/knn_pct_log2_category_by_lineage_dist.pdf", p)

df_2 = dat_y_uniq

####################################################################
### correlation between step-1 (pca dist) and step-2 (mean % of kNN)

df$pair = sort_pair(df$A, df$B)
x = df %>% select(pair, pca_dist, lineage_dist) %>% left_join(dat_y_uniq %>% select(pair, mean_pct), by = "pair")

x$log2_mean_pct = log2(x$mean_pct + 1)
x$log2_pca_dist = log2(x$pca_dist)

p = ggplot() +
    geom_point(data = x, aes(x = log2_pca_dist, y = log2_mean_pct, color = lineage_dist)) +
    labs(x = "Log2(Euclidean distances in the PCA space + 1)", y = "Log2(Mean % of kNN cells from each other + 1)") +
    theme_classic(base_size = 12) +
    theme(legend.position="none") +
    scale_color_viridis(discrete=TRUE) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
ggsave("~/share/PCA_vs_MeanPCT.pdf", p, height = 5, width = 5)

cor.test(x$log2_mean_pct, x$log2_pca_dist, method = "spearman")
### rho = -0.5813561; p-value < 2.2e-16






#########################################################################
### Repeat the above analysis but excluding the nearest neighboring pairs

########################################################################################
### Step-3: Elucidean distance in the PCA space between gastruloids vs. lineage distance

### my cell-to-gastruloid assignments
cell_assign_group = readRDS(paste0(work_path, "/data_mGASv5/tree_sanjay/assignment/assignment_well.rds"))
well_include = as.vector(cell_assign_group$well[cell_assign_group$cell_num >= 50])
cell_assign = readRDS(paste0(work_path, "/data_mGASv5/tree_sanjay/assignment/assignment_cell.rds"))
cell_assign = cell_assign[cell_assign$well %in% well_include,]
colnames(cell_assign) = c("Cell","node", "Well", "celltype")

### tree of trees
tree = readRDS(paste0(work_path, "/data_mGASv5/tree/tree_well_tape_dominant_edits.rds"))

### lineage distance
dat_i = castor::get_all_pairwise_distances(tree, 
                                           only_clades = 1:length(tree$tip.label), 
                                           as_edge_counts = TRUE)
rownames(dat_i) = colnames(dat_i) = 1:length(tree$tip.label)
dat_i = melt(as.matrix(dat_i))
dat_i$A = tree$tip.label[dat_i$Var1]
dat_i$B = tree$tip.label[dat_i$Var2]
dat_i$lineage_dist = dat_i$value
dat_i$Var1 = dat_i$Var2 = dat_i$value = NULL

### Elucidean distance on PCA 
### pca coordinates after projecting onto mGASv4 dataset
pca_coor = readRDS(paste0(work_path, "/data_mGASv5/pca/pca_coor.rds"))
irlba_res = readRDS(paste0(work_path, "/data_mGASv5/pca/PCA_pseudobulk.rds"))
preproc_res = irlba_res$x
row.names(preproc_res) = as.vector(pca_coor$sample_id)
pca_coor = as.matrix(preproc_res)

euclidean_distance <- function(x, y) {
    sqrt(sum((x - y)^2))
}

###
gastruloid_list = rownames(pca_coor)
pca_dist = NULL
for(i in 1:(length(gastruloid_list)-1)){
    for(j in (i+1):length(gastruloid_list)){
        pca_dist = rbind(pca_dist, data.frame(A = gastruloid_list[i],
                                              B = gastruloid_list[j],
                                              pca_dist = euclidean_distance(pca_coor[i,], pca_coor[j,])))
    }
}

df = dat_i %>% left_join(pca_dist, by = c("A","B")) %>% filter(!is.na(pca_dist))
### 3828 pairs (88 x 88 "unique" pairs)


### excluding those nearest neighbors
dat_x = dat_i %>% filter(A != B, A %in% c(df$A, df$B), B %in% c(df$A, df$B)) %>% 
    group_by(A) %>% slice_min(order_by = lineage_dist, n = 1)

sort_pair <- function(x, y) {
    apply(cbind(x, y), 1, function(row) paste(sort(row), collapse = "-"))
}
dat_x$pair = sort_pair(dat_x$A, dat_x$B) ### 76 unique "nearest" pairs

df$pair = sort_pair(df$A, df$B)

df = df[!df$pair %in% dat_x$pair,] ### 3752 pairs left
df$lineage_dist = factor(df$lineage_dist, levels = 4:23)


cor.test(df$pca_dist, as.numeric(df$lineage_dist), method = "spearman")
### rho = 0.04553602; p-value = 0.005275

### plotting it using log2 scale on y-axis
df$log2_pca_dist = log2(df$pca_dist)

p_1 = ggplot(df, aes(lineage_dist, log2_pca_dist, fill = lineage_dist)) + geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2) +
    labs(x="After excluding 76 nearest neighboring pairs", y="Log2(Euclidean distances in the PCA space)", title = "Rho = 0.0455, P-val = 5e-3") +
    theme_classic(base_size = 12) +
    theme(legend.position="none") +
    scale_fill_viridis(discrete=TRUE) +
    theme(axis.text.x = element_text(color="black", angle = 45), axis.text.y = element_text(color="black"))

df_1 = df


##############################
### Step-4: kNN based approach

### Present the heterogenity of transcriptomes aligning with their lineages

work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape"
source("~/work/scripts/utils.R")

### read mGASv5 dataset
obj = readRDS(paste0(work_path, "/data_mGASv5/obj_sctransform.rds"))
pca_coor = Embeddings(obj,reduction = "pca")
pd_v5 = readRDS(paste0(work_path, "/data_mGASv5/obj_sctransform_pd.rds"))
pd_v5$Cell = paste0(unlist(lapply(as.vector(pd_v5$cell_id), function(x) strsplit(x,"[_]")[[1]][3])),
                    gsub("mGASv5_", "", pd_v5$experiment_id))

### my cell-to-gastruloid assignments
cell_assign_group = readRDS(paste0(work_path, "/data_mGASv5/tree_sanjay/assignment/assignment_well.rds"))
well_include = as.vector(cell_assign_group$well[cell_assign_group$cell_num >= 50])
cell_assign = readRDS(paste0(work_path, "/data_mGASv5/tree_sanjay/assignment/assignment_cell.rds"))
cell_assign = cell_assign[cell_assign$well %in% well_include,]
colnames(cell_assign) = c("Cell","node", "Well", "celltype")

pd_v5_x = pd_v5 %>% left_join(cell_assign[,c("Cell","Well")])
pd_v5$Well = as.vector(pd_v5_x$Well)

### subset each Well to roughly similar cells
pd_v5_sub = pd_v5 %>% filter(!is.na(Well)) %>% as.data.frame()
rownames(pd_v5_sub) = as.vector(pd_v5_sub$cell_id)
pca_coor_sub = pca_coor[rownames(pd_v5_sub),]

### calculate kNNs
k.param = 15; nn.method = "rann"; nn.eps = 0; annoy.metric = "euclidean"
nn.ranked = Seurat:::NNHelper(
    data = pca_coor_sub,
    k = k.param + 1,
    method = nn.method,
    searchtype = "standard",
    eps = nn.eps,
    metric = annoy.metric)
nn.ranked = Indices(object = nn.ranked)
nn_matrix = nn.ranked[,-1]

resultA = NULL
for(i in 1:k.param){
    print(i)
    resultA = cbind(resultA, as.vector(pd_v5_sub$Well)[as.vector(nn_matrix[,i])])
}

dat_x = data.frame(A = rep(pd_v5_sub$Well, times = ncol(resultA)),
                   B = c(resultA))

dat = dat_x %>% group_by(A, B) %>% tally() %>%
    dcast(A~B, fill=0)
rownames(dat) = dat[,1]
dat = as.matrix(dat[,-1])

dat_norm = t(t(dat)/apply(dat, 2, sum))
dat_norm = dat_norm/apply(dat_norm, 1, sum)

### calculate average of A-B and B-A
dat_norm_mean = (dat_norm + t(dat_norm))/2

mean_pct = melt(dat_norm_mean)
colnames(mean_pct) = c("A", "B", "pct")

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
dat_i$A = gsub("[.]", "", dat_i$A)
dat_i$B = gsub("[.]", "", dat_i$B)

dat_y = dat_i %>% filter(A != B, A %in% rownames(dat_norm_mean), B %in% rownames(dat_norm_mean))
dat_y$pair = sort_pair(dat_y$A, dat_y$B)
dat_y_uniq = dat_y[!duplicated(dat_y$pair), ]

### excluding those nearest neighbors
dat_x = dat_i %>% filter(A != B, A %in% c(df$A, df$B), B %in% c(df$A, df$B)) %>% 
    group_by(A) %>% slice_min(order_by = lineage_dist, n = 1)

sort_pair <- function(x, y) {
    apply(cbind(x, y), 1, function(row) paste(sort(row), collapse = "-"))
}
dat_x$pair = sort_pair(dat_x$A, dat_x$B) ### 76 unique "nearest" pairs

dat_y_uniq = dat_y_uniq[!dat_y_uniq$pair %in% dat_x$pair,] ### 3752 pairs left


### making boxplot, with lineage distance as categories.
y = NULL
for(i in 1:nrow(dat_y_uniq)){
    y = c(y, mean_pct$pct[mean_pct$A == dat_y_uniq$A[i] & mean_pct$B == dat_y_uniq$B[i]])
}
dat_y_uniq$mean_pct = y*100
dat_y_uniq$lineage_dist = factor(dat_y_uniq$lineage_dist, levels = 4:23)

cor.test(dat_y_uniq$mean_pct, as.numeric(dat_y_uniq$lineage_dist), method = "spearman")
### rho = -0.08443167; p-value = 2.228e-07

### log2-scale
dat_y_uniq$log2_mean_pct = log2(dat_y_uniq$mean_pct + 1)
p_2 = ggplot(dat_y_uniq, aes(lineage_dist, log2_mean_pct, fill = lineage_dist)) + geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2) +
    labs(x="After excluding 76 nearest neighboring pairs", y="Log2(Mean % of kNN cells from each other + 1)", title = "Rho = -0.0844, P-val = 2e-7") +
    theme_classic(base_size = 12) +
    theme(legend.position="none") +
    scale_fill_viridis(discrete=TRUE) +
    theme(axis.text.x = element_text(color="black", angle = 45), axis.text.y = element_text(color="black"))

df_2 = dat_y_uniq

### correlation between step-1 (pca dist) and step-2 (mean % of kNN)

x = df_1 %>% select(pair, pca_dist, lineage_dist) %>% 
    left_join(df_2 %>% select(pair, mean_pct), by = "pair")

x$log2_mean_pct = log2(x$mean_pct + 1)
x$log2_pca_dist = log2(x$pca_dist)

cor.test(x$log2_mean_pct, x$log2_pca_dist, method = "spearman")
### rho = -0.5788041; p-value < 2.2e-16

p_3 = ggplot() +
    geom_point(data = x, aes(x = log2_pca_dist, y = log2_mean_pct, color = lineage_dist)) +
    labs(x = "Log2(Euclidean distances in the PCA space + 1)", y = "Log2(Mean % of kNN cells from each other + 1)", title = "Rho = -0.579, P-val < 2e-16") +
    theme_classic(base_size = 12) +
    theme(legend.position="none") +
    scale_color_viridis(discrete=TRUE) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))

ggsave("~/share/Lineage_vs_transcriptome_after_remove_NN.pdf", p_1 + p_2 + p_3, width = 15, height = 5)







