
#######################################
### Performing cell-to-well assignments
### Chengxiang Qiu
### Feb-20, 2025

###################################################################################################
### Step-1: From debris-seq, removing barcodes if they are overlapped between different gastruloids

### Of note, we took the "BC + Site1" as uniqual DTT id, 
### and for each DTT, we only retained the Site2-6 editing combinations with the most reads

### Step1File_DebrisSeqTAPE_mGASv4.csv can be downloaded from our website directly
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/sam_tape

source("help_script.R")
work_path = "Your_work_path"

dat = read.csv("Step1File_DebrisSeqTAPE_mGASv4.csv", header=F)
dat = dat[,c(1,2,3,4,5,6,7,8)]
sample = gsub('[.|_]', '', dat$V1)
dat$sample = as.vector(sample)
dat$barcode = paste0(substr(dat$V2, 1, 5), 
                     substr(dat$V2, 8, 12), '-',
                     substr(dat$V3, 1, 3), '-',
                     substr(dat$V4, 1, 3), '-',
                     substr(dat$V5, 1, 3), '-',
                     substr(dat$V6, 1, 3))
dat = dat[,c("sample", "barcode")]
dat = dat[!dat$sample %in% c("P3-E10", "P3-D9"),] ### P3-D9 and P3-E10 have the identical TapeBCs
dat = unique(dat)

write.table(dat, paste0(work_path, "/data_sci/tape_barcode_2/debris_sample_barcode.txt"), row.names=F, col.names=F, sep="\t", quote=F)

### here, to be more specific for cell-to-gastruloid assignments, we excluded TapeBCs which are overlapped by different gastruloids
dat_x = table(dat$barcode)
exclude_barcode = names(dat_x)[dat_x > 1]
dat = dat[!dat$barcode %in% exclude_barcode,]

write.table(dat, paste0(work_path, "/data_sci/tape_barcode_2/debris_sample_barcode_uniq.txt"), row.names=F, col.names=F, sep="\t", quote=F)
### 4,995 TapeBCs across 138 wells


#########################################
### Step-2: Creating TapeBC x Cell matrix

### python Cell_to_gastruloid.py

dat = read.table(paste0(work_path, "/data_sci/tape_barcode_3/Tape_cell_uniq.txt"), as.is=T)
dat = dcast(dat, V1 ~ V2, fill=0)
rownames(dat) = dat[,1]
dat = dat[,-1]
count = as(as.matrix(dat), "sparseMatrix") 
### TapeBC x Cell matrix
### 4,609 TapeBCs x 213,171 cells

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
        ggsave(paste0(work_path, "/data_sci/tape_barcode_3/scatter_bc_count.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)


saveRDS(count, paste0(work_path, "/data_sci/tape_barcode_3/count.rds"))
saveRDS(pd, paste0(work_path, "/data_sci/tape_barcode_3/pd.rds"))


#######################################################################
### Step-3: Assinging individual cells to their best matched gastruloid

count = readRDS(paste0(work_path, "/data_sci/tape_barcode_3/count.rds"))
pd = readRDS(paste0(work_path, "/data_sci/tape_barcode_3/pd.rds"))

keep = pd$bc_count >= 3
count_sub = count[,keep]
pd_sub = pd[keep,]

### After filtering cells with fewer than 3 TapeBCs detected:
### 4,609 TapeBCs x 154,988 cells

### of note, here I didn't consider read num, but just "yes" or "no" matching
count_sub_nonzero = count_sub
count_sub_nonzero@x[count_sub_nonzero@x > 1] = 1

debris_set = read.table(paste0(work_path, "/data_sci/tape_barcode_2/debris_sample_barcode_uniq.txt"), as.is=T)
names(debris_set) = c("Well", "BC")
fd = data.frame(BC = rownames(count_sub)) %>% left_join(debris_set, by = "BC") %>% 
    group_by(BC, Well) %>% tally() %>% dcast(Well~BC, fill=0)
rownames(fd) = fd[,1]
fd = as.matrix(fd[,-1])
fd = fd[,rownames(count_sub_nonzero)]

### Well x TapeBCs %*% TapeBCs x Cell = Well x Cell matrix
dat = fd %*% count_sub_nonzero
dat_norm = t(t(dat)/apply(dat, 2, sum))

top_each_column = apply(dat_norm, 2, max)
second_each_column = apply(dat_norm, 2, function(x) sort(x, decreasing=T)[2])
top_each_column_index = rownames(dat_norm)[apply(dat_norm, 2, which.max)]

best_match = data.frame(well = top_each_column_index,
                        ratio = top_each_column/second_each_column)
best_match$if_sig = best_match$ratio >= 1.5
### n = 129,853 vs. 25,135

saveRDS(best_match, paste0(work_path, "/data_sci/tape_barcode_3/best_match.rds"))


######################################################################
### Step-4: Performing dimension reduction on the TapeBC x Cell matrix

count = readRDS(paste0(work_path, "/data_sci/tape_barcode_3/count.rds"))
pd = readRDS(paste0(work_path, "/data_sci/tape_barcode_3/pd.rds"))

keep = pd$bc_count >= 3
count_sub = count[,keep]
pd_sub = pd[keep,]
count_sub = count_sub[rowSums(count_sub) >= 10,]
### 3,727 TapeBCs x 154,988 cells

obj = CreateSeuratObject(count_sub, meta.data = pd_sub)
obj = NormalizeData(obj, normalization.method = "RC", scale.factor = 10000)
obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = nrow(obj))
obj = ScaleData(object = obj, verbose = FALSE)
obj = RunPCA(object = obj, npcs = 100, verbose = FALSE)
obj = FindNeighbors(object = obj, dims = 1:100, reduction = "pca", k.param = 10)
obj = FindClusters(object = obj, resolution = 1)
obj = RunUMAP(object = obj, reduction = "pca", dims = 1:100, min.dist = 0.1, n.neighbors = 10, n.components = 2)

df = data.frame(obj[[]])
df$UMAP_1 = Embeddings(obj, reduction = "umap")[,1]
df$UMAP_2 = Embeddings(obj, reduction = "umap")[,2]

try(ggplot() +
        geom_point(data = df, aes(x = UMAP_1, y = UMAP_2), color = "black", size=0.2) +
        geom_point(data = df, aes(x = UMAP_1, y = UMAP_2, color = RNA_snn_res.1), size=0.1) +
        #ggrepel::geom_text_repel(data = df %>% group_by(RNA_snn_res.1) %>% sample_n(1), aes(x = UMAP_1, y = UMAP_2, label = RNA_snn_res.1), color = "black", size = 3, family = "Arial") +
        theme_void() +
        theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5)) +
        ggsave(paste0(work_path, "/data_sci/tape_barcode_3/UMAP_cell_cluster.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)

saveRDS(obj, paste0(work_path, "/data_sci/tape_barcode_3/obj.rds"))
saveRDS(df, paste0(work_path, "/data_sci/tape_barcode_3/obj_pd.rds"))
saveRDS(Embeddings(obj, reduction = "pca"), paste0(work_path, "/data_sci/tape_barcode_3/obj_pca_coor.rds"))

best_match = readRDS(paste0(work_path, "/data_sci/tape_barcode_3/best_match.rds"))
print(sum(rownames(best_match) == rownames(df)))

df$well = as.vector(best_match$well)
df$if_sig = as.vector(best_match$if_sig)
df$well[!df$if_sig] = NA


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

saveRDS(df_2, paste0(work_path, "/data_sci/tape_barcode_3/best_match_step2.rds"))

df_x = rbind(df_1, df_2)
df_x = df_x[rownames(df),]

df_x = df_x[df_x$if_sig,]
df_x = df_x[,c("cell_id", "well")]
df_x_num = table(df_x$well)
### 149,216 cells across 131 wells

df_x_null = df_x[df_x$well %in% names(df_x_num)[df_x_num < 100],]

df_x = df_x[df_x$well %in% names(df_x_num)[df_x_num >= 100],]
### 148,724 cells across 121 wells

saveRDS(df_x, paste0(work_path, "/data_sci/tape_barcode_3/cell_assign.rds"))

df = readRDS(paste0(work_path, "/data_sci/tape_barcode_3/obj_pd.rds"))
df = df %>% left_join(df_x, by = "cell_id")
try(ggplot() +
        geom_point(data = df[is.na(df$well),], aes(x = UMAP_1, y = UMAP_2), color = "grey80", size=0.2) +
        geom_point(data = df[!is.na(df$well),], aes(x = UMAP_1, y = UMAP_2, color = well), size=0.1) +
        #ggrepel::geom_text_repel(data = df %>% group_by(RNA_snn_res.1) %>% sample_n(1), aes(x = UMAP_1, y = UMAP_2, label = RNA_snn_res.1), color = "black", size = 3, family = "Arial") +
        theme_void() +
        theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5)) +
        ggsave(paste0(work_path, "/data_sci/tape_barcode_3/UMAP_best_match.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)

pd = readRDS(paste0(work_path, "/data_sci/obj_processed_pd.rds"))
df_x = df_x %>% left_join(pd %>% select(cell_id = sample, celltype), by = "cell_id")
saveRDS(df_x, paste0(work_path, "/data_sci/tape_barcode_3/cell_assign.rds"))

print(nrow(df_x)/nrow(pd)) ### 148724/247064 = 60.2%
well_table = table(df_x$well)
### mean = 1,229; median = 1,066; 100 - 5,883


### jay's question: In aggregate, are cells assigned to the 17 wells with <100 cells biased towards any particular cell type?
df_x_null = df_x_null %>% left_join(pd %>% select(cell_id = sample, celltype), by = "cell_id")


#########################################################################################
### Step-6: Calculating the jaccard similarity between each gastruloid and each clonotype

cell_assign = readRDS(paste0(work_path, "/data_sci/tape_barcode_3/cell_assign.rds"))
count = readRDS(paste0(work_path, "/data_sci/tape_barcode_3/count.rds"))

clonotype_list = unique(cell_assign$well) ### n = 121
clonotype_BC = list()

BC_num = data.frame()
for(clonotype_i in clonotype_list){
    print(clonotype_i)
    count_i = count[,colnames(count) %in% as.vector(cell_assign$cell_id[cell_assign$well == clonotype_i])]
    count_i_pct = rowMeans(count_i > 0)
    clonotype_BC[[clonotype_i]] = names(count_i_pct[count_i_pct >= 0.05])
    BC_num = rbind(BC_num, data.frame(well = clonotype_i,
                                      BC_num = length(unique(clonotype_BC[[clonotype_i]]))))
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

debris_set = read.table(paste0(work_path, "/data_sci/tape_barcode_2/debris_sample_barcode_uniq.txt"), as.is=T)
names(debris_set) = c("Well", "BC")

BC_num_x = debris_set %>% group_by(Well) %>% tally()


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

row_names = sort(rownames(dat_mat))
dat_mat = dat_mat[row_names, row_names]

### making heatmap with Jaccard distance

library("gplots")
library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
pdf(paste0(work_path, "/data_sci/tape_barcode_3/Jaccard_similarity_heatmap.pdf"), 8, 5)
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

write.table(row_names, paste0(work_path, "/data_sci/tape_barcode_3/Jaccard_similarity_heatmap.rownames.txt"), row.names=F, col.names=F, quote=F, sep="\t")




##################################################
### Step-7: Can we do a cumulative curve analysis? 

### For example, if we rank the gastruloids from highest to lowest in terms of 
### how many ExE cells they have, how many gastruloids do we need to explain X% of the observed ExE cells. My guess is that there will be at least some cell types that are largely just in a handful of gastruloids.

cell_assign = readRDS(paste0(work_path, "/data_sci/tape_barcode_3/cell_assign.rds"))

df = cell_assign %>% group_by(celltype, well) %>% tally() %>%
    left_join(cell_assign %>% group_by(celltype) %>% tally() %>% rename(total_n = n), by = "celltype") %>%
    mutate(frac = 100*n/total_n) %>% select(celltype, well, frac) %>% as.data.frame() %>%
    group_by(celltype) %>% arrange(desc(frac), .by_group = T) %>% as.data.frame() %>%
    group_by(celltype) %>% summarize(cum_frac = cumsum(frac))

celltype_num = table(df$celltype)
x = NULL
for(i in 1:length(celltype_num)){
    x = c(x, 1:celltype_num[i])
}
df$x_axis = as.vector(x)

df$celltype = factor(df$celltype, levels = names(gastruloid_celltype_color_code)[names(gastruloid_celltype_color_code) %in% df$celltype])

ggplot(df, aes(x = x_axis, y = cum_frac, color = celltype)) + 
    #geom_point() + 
    geom_line(linewidth = 1) +
    labs(title = "", x = "gastruloids", y = "cumulative % of cells") +
    theme_classic(base_size = 15) +
    #scale_x_continuous(limits = c(1, 456), breaks = c(1:456)) +
    theme(legend.position="none") +
    scale_color_manual(values=gastruloid_celltype_color_code) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    ggsave(paste0(work_path, "/data_sci/tape_barcode_3/cumulative_frac_gastruloids_for_each_celltype.pdf"),
           dpi = 300,
           height  = 4, 
           width = 6)








