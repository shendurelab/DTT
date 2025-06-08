
##################################################################
### Cell-to-Well assignments in the "eight gastruloids" experiment
### Chengxiang Qiu
### Feb-20, 2025

### For this "eight gastruloids" experiment, performing cell-to-gastruloid assignments for clones separately

### Clone-05: Well14, Well17, Well25
### Clone-25: Well01, Well03, Well21
### Clone-32: Well16, Well28

### Debris-seq DTT data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/sam_tape
### 1) clone05: Step1File_DebrisSeqTAPE_mGASv3_240329.csv
### 2) clone25: Step1File_DebrisSeqTAPE_mGASv5_C25_v240331.csv
### 3) clone32: step1file_debrisseqtape_mgasv5_c32_v240331.csv

### scRNA-seq DTT data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/sam_tape
### 1) mGASv3_Lane1_CellByTape_10X_bamExtractV3_t100_Site1collapsed2_v240330.csv
### 2) mGASv3_Lane2_CellByTape_10X_bamExtractV3_t100_Site1collapsed2_v240331.csv

### Some intermediate profiles are provided in case you prefer not to regenerate them
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/sam_tape/data_mGASv3
### 1) cell_id.txt ### line 277, 289
### 2) cell_assign.rds ### line 292

clone_list = data.frame(clone_id = c(rep("Clone-05", 3), rep("Clone-25", 3), rep("Clone-32", 2)),
                        well_id = c("Well14", "Well17", "Well25", "Well01", "Well03", "Well21", "Well16", "Well28"))

##################################################################
### Step-1: investigating the overlap of TapeBCs within each clone

source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

work_path = "Your_work_path"

clone_i = "Clone-05"
dat = read.csv(paste0(work_path, "/data_mGASv3/TapeBC_data/Step1File_DebrisSeqTAPE_mGASv3_240329.csv"), header=F)
dat = dat[,c(1,2,3,4,5,6,7,8)]
sample = gsub('[.|_]', '', dat$V1)
dat$sample = as.vector(sample)
dat$barcode = paste0(substr(dat$V2, 1, 5), 
                     substr(dat$V2, 8, 12), '-',
                     substr(dat$V3, 1, 3), '-',
                     substr(dat$V4, 1, 3), '-',
                     substr(dat$V5, 1, 3), '-',
                     substr(dat$V6, 1, 3))
dat_1 = dat[,c("sample", "barcode")]

clone_i = "Clone-25"
dat = read.csv(paste0(work_path, "/data_mGASv3/TapeBC_data/Step1File_DebrisSeqTAPE_mGASv5_C25_v240331.csv"), header=F)
dat = dat[,c(1,2,3,4,5,6,7,8)]
sample = gsub('[.|_]', '', dat$V1)
dat$sample = as.vector(sample)
dat$barcode = paste0(substr(dat$V2, 1, 5), 
                     substr(dat$V2, 8, 12), '-',
                     substr(dat$V3, 1, 3), '-',
                     substr(dat$V4, 1, 3), '-',
                     substr(dat$V5, 1, 3), '-',
                     substr(dat$V6, 1, 3))
dat_2 = dat[,c("sample", "barcode")]

clone_i = "Clone-32"
dat = read.csv(paste0(work_path, "/data_mGASv3/TapeBC_data/step1file_debrisseqtape_mgasv5_c32_v240331.csv"), header=F)
dat = dat[,c(1,2,3,4,5,6,7,8)]
sample = gsub('[.|_]', '', dat$V1)
dat$sample = as.vector(sample)
dat$barcode = paste0(substr(dat$V2, 1, 5), 
                     substr(dat$V2, 8, 12), '-',
                     substr(dat$V3, 1, 3), '-',
                     substr(dat$V4, 1, 3), '-',
                     substr(dat$V5, 1, 3), '-',
                     substr(dat$V6, 1, 3))
dat_3 = dat[,c("sample", "barcode")]

dat = rbind(dat_1, dat_2, dat_3)
write.table(dat, paste0(work_path, "/data_mGASv3/cell_to_gastruloid/debris_sample_barcode.txt"), row.names=F, col.names=F, sep="\t", quote=F)



#########################################
### Step-2: Creating TapeBC x Cell matrix

### python Cell_to_gastruloid.py

dat = read.table(paste0(work_path, "/data_mGASv3/cell_to_gastruloid/Tape_cell.txt"), as.is=T)
dat = dcast(dat, V1 ~ V2, fill=0)
rownames(dat) = dat[,1]
dat = dat[,-1]
count = as(as.matrix(dat), "sparseMatrix") 
### 470 TapeBC x 9,358 cells (out of 9,929 cells in total)

count_nonzero = count
count_nonzero@x[count_nonzero@x > 1] = 1

pd = data.frame(cell_id = colnames(count),
                bc_count = colSums(count_nonzero),
                read_count = colSums(count))

try(ggplot() +
        geom_point(data = pd %>% mutate(log2_bc_count = log2(bc_count), log2_read_count = log2(read_count)),
                   aes(x = log2_bc_count, y = log2_read_count), size=0.1) +
        theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5)) +
        ggsave(paste0(work_path, "/data_mGASv3/cell_to_gastruloid/scatter_bc_count.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)


saveRDS(count, paste0(work_path, "/data_mGASv3/cell_to_gastruloid/count.rds"))
saveRDS(pd, paste0(work_path, "/data_mGASv3/cell_to_gastruloid/pd.rds"))


#######################################################################
### Step-3: Assinging individual cells to their best matched gastruloid

count = readRDS(paste0(work_path, "/data_mGASv3/cell_to_gastruloid/count.rds"))
pd = readRDS(paste0(work_path, "/data_mGASv3/cell_to_gastruloid/pd.rds"))

keep = pd$read_count >= 10 & pd$bc_count >= 2 ### Of note, in the sci-experiment, I used bc_count >= 3 as cutoff
count_sub = count[,keep]
pd_sub = pd[keep,]
### 470 TapeBCs x 8547 cells

count_sub_nonzero = count_sub
count_sub_nonzero@x[count_sub_nonzero@x > 1] = 1

debris_set = read.table(paste0(work_path, "/data_mGASv3/cell_to_gastruloid/debris_sample_barcode.txt"), as.is=T)
names(debris_set) = c("Well", "BC")
fd = data.frame(BC = rownames(count_sub)) %>% left_join(debris_set, by = "BC") %>% 
    group_by(BC, Well) %>% tally() %>% dcast(Well~BC, fill=0)
rownames(fd) = fd[,1]
fd = as.matrix(fd[,-1])
fd = fd[,rownames(count_sub_nonzero)]

dat = fd %*% count_sub_nonzero
dat_norm = t(t(dat)/apply(dat, 2, sum))

saveRDS(dat, paste0(work_path, "/data_mGASv3/cell_to_gastruloid/well_cell_match_num.rds"))

top_each_column = apply(dat_norm, 2, max)
second_each_column = apply(dat_norm, 2, function(x) sort(x, decreasing=T)[2])
top_each_column_index = rownames(dat_norm)[apply(dat_norm, 2, which.max)]

best_match = data.frame(well = top_each_column_index,
                        ratio = top_each_column/second_each_column)
best_match$if_sig = best_match$ratio >= 1.5
table(best_match$if_sig)
### 8409 vs. 138

saveRDS(best_match, paste0(work_path, "/data_mGASv3/cell_to_gastruloid/best_match.rds"))


#####################################################################
### Step-4: Performing dimension reduction on the TapeBC x Cell matrix

count = readRDS(paste0(work_path, "/data_mGASv3/cell_to_gastruloid/count.rds"))
pd = readRDS(paste0(work_path, "/data_mGASv3/cell_to_gastruloid/pd.rds"))

keep = pd$read_count >= 10 & pd$bc_count >= 2 ### Of note, in the sci-experiment, I used bc_count >= 3 as cutoff
count_sub = count[,keep]
pd_sub = pd[keep,]
### 470 TapeBCs x 8547 cells

count_sub = count_sub[rowSums(count_sub) >= 10,]
### 440 TapeBCs x 8547 cells

obj = CreateSeuratObject(count_sub, meta.data = pd_sub)
obj = NormalizeData(obj, normalization.method = "RC", scale.factor = 10000)
obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = nrow(obj))
obj = ScaleData(object = obj, verbose = FALSE)
obj = RunPCA(object = obj, npcs = 30, verbose = FALSE)
obj = FindNeighbors(object = obj, dims = 1:30, reduction = "pca", k.param = 10)
obj = FindClusters(object = obj, resolution = 1)
obj = RunUMAP(object = obj, reduction = "pca", dims = 1:30, min.dist = 0.3, n.neighbors = 10, n.components = 2)

df = data.frame(obj[[]])
df$UMAP_1 = Embeddings(obj, reduction = "umap")[,1]
df$UMAP_2 = Embeddings(obj, reduction = "umap")[,2]

best_match = readRDS(paste0(work_path, "/data_mGASv3/cell_to_gastruloid/best_match.rds"))
print(sum(rownames(best_match) == rownames(df)))

df$well = as.vector(best_match$well)
df$if_sig = as.vector(best_match$if_sig)
df$well[!df$if_sig] = NA

df_mean = df %>% filter(!is.na(well)) %>%
    group_by(well) %>% summarize(UMAP_1_mean = mean(UMAP_1), UMAP_2_mean = mean(UMAP_2))

try(ggplot() +
        geom_point(data = df, aes(x = UMAP_1, y = UMAP_2), color = "black", size=0.3) +
        geom_point(data = df, aes(x = UMAP_1, y = UMAP_2), color = "grey80", size=0.2) +
        geom_point(data = df %>% filter(!is.na(well)), aes(x = UMAP_1, y = UMAP_2, color = well), size=0.2) +
        ggrepel::geom_text_repel(data = df_mean, aes(x = UMAP_1_mean, y = UMAP_2_mean, label = well), color = "black", size = 3, family = "Arial") +
        theme_void() +
        theme(legend.position="none") +
        scale_color_manual(values=gastruloid_celltype_color_code) +
        theme(plot.title = element_text(hjust = 0.5)) +
        ggsave(paste0(work_path, "/data_mGASv3/cell_to_gastruloid/UMAP_cell_cluster.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)

saveRDS(obj, paste0(work_path, "/data_mGASv3/cell_to_gastruloid/obj.rds"))
saveRDS(df, paste0(work_path, "/data_mGASv3/cell_to_gastruloid/obj_pd.rds"))
saveRDS(Embeddings(obj, reduction = "pca"), paste0(work_path, "/data_mGASv3/cell_to_gastruloid/obj_pca_coor.rds"))



#######################################################################################################
### Step-5: For cells which are not assigned by the first step, we assign them based on their Knn cells

pca_coor = Embeddings(obj, reduction = "pca")
pca_coor_1 = pca_coor[df$if_sig,]
pca_coor_2 = pca_coor[!df$if_sig,]
df_1 = df[df$if_sig,]
df_2 = df[!df$if_sig,]

k.param = 10; nn.method = "rann"; nn.eps = 0; annoy.metric = "euclidean"
nn.ranked = Seurat:::NNHelper(
    data = pca_coor_1,
    query = pca_coor_2,
    k = k.param,
    method = nn.method,
    searchtype = "standard",
    eps = nn.eps,
    metric = annoy.metric)
nn.ranked = Indices(object = nn.ranked)
nn_matrix = nn.ranked

df_2$well = df_2$if_sig = NULL
resultA = NULL
for(i in 1:k.param){
    print(i)
    resultA = cbind(resultA, as.vector(df_1$well)[as.vector(nn_matrix[,i])])
}
resultB = lapply(1:nrow(resultA), function(x){
    return(names(sort(table(resultA[x,]), decreasing = T))[1])
})
df_2$well = unlist(resultB)
resultB = lapply(1:nrow(resultA), function(x){
    return(sort(table(resultA[x,]), decreasing = T)[1])
})
df_2$top = unlist(resultB)
resultB = lapply(1:nrow(resultA), function(x){
    return(sort(table(resultA[x,]), decreasing = T)[2])
})
df_2$second = unlist(resultB)

df_2$if_sig = df_2$top == 10 | (!is.na(df_2$second) & df_2$top/df_2$second >= 1.5)
df_2$top = df_2$second = NULL
df_2$well[!df_2$if_sig] = NA

tmp = readRDS(paste0(work_path, "/data_mGASv3/cell_to_gastruloid/well_cell_match_num.rds"))
tmp = melt(as.matrix(tmp))
colnames(tmp) = c("well","cell_id","num")
df_2_tmp = df_2 %>% left_join(tmp, by = c("cell_id","well"))
df_2_tmp = df_2_tmp[df_2_tmp$num == 0 & !is.na(df_2_tmp$num),]
df_2$if_sig[df_2$cell_id %in% as.vector(df_2_tmp$cell_id)] = FALSE
df_2$well[df_2$cell_id %in% as.vector(df_2_tmp$cell_id)] = NA
### only 2 cells have not been assigned :)

df_x = rbind(df_1, df_2)
df_x = df_x[rownames(df),]

df_x = df_x[df_x$if_sig,]
df_x = df_x[,c("cell_id", "well")]

saveRDS(df_x, paste0(work_path, "/data_mGASv3/cell_to_gastruloid/cell_assign_before_filter.rds"))

pd = read.table(paste0(work_path, "/data_mGASv3/cell_id.txt"), as.is=T, sep="\t")
colnames(pd) = c("cell_id","cell","celltype")
df_x_out = df_x %>% left_join(pd %>% select(cell_id, celltype), by = "cell_id")
saveRDS(df_x_out, paste0(work_path, "/data_mGASv3/cell_to_gastruloid/cell_assign_before_filter.rds"))


df_x_num = table(df_x$well)
df_x = df_x[df_x$well %in% names(df_x_num)[df_x_num >= 50],]

saveRDS(df_x, paste0(work_path, "/data_mGASv3/cell_to_gastruloid/cell_assign.rds"))


pd = read.table(paste0(work_path, "/data_mGASv3/cell_id.txt"), as.is=T, sep="\t")
colnames(pd) = c("cell_id","cell","celltype")
df_x = df_x %>% left_join(pd %>% select(cell_id, celltype), by = "cell_id")
saveRDS(df_x, paste0(work_path, "/data_mGASv3/cell_to_gastruloid/cell_assign.rds"))
### 8545 cells have been assigned

### Well01 Well03 Well14 Well16 Well17 Well21 Well25 Well28
###   2696   1151   1741    449    625    993    312    578

###################################################################################################
### Step-6: Validation: calculate the jaccard similarity between each gastruloid and each clonotype

cell_assign = readRDS(paste0(work_path, "/data_mGASv3/cell_to_gastruloid/cell_assign.rds"))
count = readRDS(paste0(work_path, "/data_mGASv3/cell_to_gastruloid/count.rds"))

clonotype_list = unique(cell_assign$well)
clonotype_BC = list()

for(clonotype_i in clonotype_list){
    print(clonotype_i)
    count_i = count[,colnames(count) %in% as.vector(cell_assign$cell_id[cell_assign$well == clonotype_i]), drop=FALSE]
    count_i_pct = rowMeans(count_i > 0)
    clonotype_BC[[clonotype_i]] = names(count_i_pct[count_i_pct >= 0.05])
}


jaccard_dist = function(x, y){
    a = length(intersect(x, y))
    b = length(union(x, y))
    return(a/b)
}

dat = NULL
for(i in 1:(length(clonotype_list)-1)){
    print(i)
    for(j in (i+1):length(clonotype_list)){
        dat = rbind(dat, data.frame(A = clonotype_list[i],
                                    B = clonotype_list[j],
                                    dist = jaccard_dist(clonotype_BC[clonotype_list[i]], clonotype_BC[clonotype_list[j]])))
    }
}

debris_set = read.table(paste0(work_path, "/data_mGASv3/cell_to_gastruloid/debris_sample_barcode.txt"), as.is=T)
names(debris_set) = c("Well", "BC")

### Of note, the resulting Jac-dis is very low, 
### because only ~25% of TapeBCs of Debris-seq have been detected in the scRNA-seq;
### If we want to make the Jac-dis looks bigger, we could use subset of TapeBCs of Debris-seq

debris_set = debris_set[debris_set$BC %in% rownames(count),]

well_BC = list()
well_list = unique(cell_assign$well)
for(well_i in well_list){
    well_BC[[well_i]] = as.vector(debris_set$BC[debris_set$Well == well_i])
    print(length(well_BC[[well_i]]))
}

dat_x = NULL
for(i in clonotype_list){
    print(i)
    for(j in well_list){
        dat_x = rbind(dat_x, data.frame(clonotype = i,
                                        well = j,
                                        dist = jaccard_dist(clonotype_BC[[i]], well_BC[[j]])))
    }
}
dat_x = dat_x[order(dat_x$dist, decreasing = T),]

dat_mat = dat_x %>% dcast(clonotype~well)
rownames(dat_mat) = dat_mat[,1]
dat_mat = dat_mat[,-1]

row_names = c("Well14","Well17","Well25","Well01","Well03","Well21","Well16","Well28")
dat_mat = dat_mat[row_names, row_names]

### making heatmap with Jaccard distance

library("gplots")
library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
pdf(paste0(work_path, "/data_mGASv3/cell_to_gastruloid/Jaccard_similarity_heatmap.pdf"), 8, 5)
heatmap.2(as.matrix(dat_mat), 
          col=Colors, 
          scale="none", 
          Rowv = F, 
          Colv = F, 
          key=T, 
          density.info="none", 
          trace="none", 
          cexRow = 1, 
          cexCol = 1,
          margins = c(5,5))
dev.off()

write.table(row_names, paste0(work_path, "/data_mGASv3/cell_to_gastruloid/Jaccard_similarity_heatmap.rownames.txt"), row.names=F, col.names=F, quote=F, sep="\t")


#######################################################################################################
### Step-7: After assignment, making cell-type-composition across wells, and UMAPs for individual wells

cell_assign = readRDS(paste0(work_path, "/cell_to_gastruloid/cell_assign.rds"))
pd = readRDS(paste0(work_path, "/obj_processed_pd.rds"))
pd$cell_id = paste0(unlist(lapply(rownames(pd), function(x) strsplit(x,"[_]")[[1]][3])),
                    gsub("mGASv3_mGas", "", as.vector(pd$experiment_id)))
df = pd %>% left_join(cell_assign[,c("cell_id","well")], by = "cell_id")
df_sub = df[!is.na(df$well),]

try(ggplot() +
        geom_point(data = df, aes(x = UMAP_1, y = UMAP_2), size=0.2, color = "white") +
        geom_point(data = df[sample(1:nrow(df)),], aes(x = UMAP_1, y = UMAP_2), size=0.3, color = "black") +
        geom_point(data = df[sample(1:nrow(df)),],
                   aes(x = UMAP_1, y = UMAP_2, color = well), size=0.2) +
        # ggrepel::geom_text_repel(data = pd %>% group_by(celltype) %>% sample_n(1), aes(x = UMAP_1, y = UMAP_2, label = celltype), color = "black", size = 3, family = "Arial") +
        theme_void() +
        theme(legend.position="none") +
        # labs(title = "mGASv3") +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values=eight_wells_color_code) +
        ggsave(paste0("~/share/UMAP_well.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)


Well_list = unique(df$well[!is.na(df$well)])

for(i in Well_list){
    print(i)
    try(ggplot() +
            geom_point(data = df, aes(x = UMAP_1, y = UMAP_2), size=0.2, color = "grey90") +
            geom_point(data = df[df$well == i & !is.na(df$well),],
                       aes(x = UMAP_1, y = UMAP_2, color = celltype), size=1) +
            # ggrepel::geom_text_repel(data = pd %>% group_by(celltype) %>% sample_n(1), aes(x = UMAP_1, y = UMAP_2, label = celltype), color = "black", size = 3, family = "Arial") +
            theme_void() +
            theme(legend.position="none") +
            # labs(title = "mGASv3") +
            theme(plot.title = element_text(hjust = 0.5)) +
            scale_color_manual(values=gastruloid_celltype_color_code) +
            ggsave(paste0("~/share/", "UMAP_", i,"_color_celltype.png"),
                   dpi = 300,
                   height  = 5, 
                   width = 5), silent = TRUE)
    
}


### cell number across cell types within each well

pd_sub = df %>% filter(!is.na(well)) %>%
    group_by(well, celltype) %>% tally() %>% mutate(cell_num = n)
pd_sub$celltype = factor(pd_sub$celltype, levels = names(gastruloid_celltype_color_code)[names(gastruloid_celltype_color_code) %in% pd_sub$celltype])
pd_sub$well = factor(pd_sub$well, levels = c("Well14","Well17","Well25","Well01","Well03","Well21","Well16","Well28"))


p = pd_sub %>%
    ggplot(aes(x=well, y=n, fill=celltype)) +
    geom_bar(stat="identity") +
    labs(x="", y="", title="") +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    scale_fill_manual(values=gastruloid_celltype_color_code) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    ggsave("~/share/cell_num_Well.pdf",
           height  = 5, 
           width = 8)



