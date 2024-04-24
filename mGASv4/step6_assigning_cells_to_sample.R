### aligning cells to debris-seq wells

work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape/data_sci"
source("~/work/scripts/jax/atac_seq/utils_atac.R")

##################################
### from debris-seq, removing barcodes if they are existed in over 10 gastruloid samples

dat = read.csv(paste0(work_path, "/tape_barcode_2/Step1File_DebrisSeqTAPE_mGASv4.csv"), header=F)
dat = dat[,c(1,2)]
sample = gsub('[.|_]', '', dat$V1)
dat$sample = as.vector(sample)
#dat$barcode = unlist(lapply(as.vector(dat$V2), function(x) strsplit(x,"[-]")[[1]][1])) 
dat$barcode = as.vector(dat$V2)
dat = dat[,c("sample", "barcode")]
dat = unique(dat)

pdf("~/share/debris_sample.pdf")
hist(table(dat$sample), 50)
dev.off()

pdf("~/share/debris_barcode.pdf")
hist(table(dat$barcode), 50)
dev.off()

dat_x = table(dat$barcode)
exclude_barcode = names(dat_x)[dat_x >= 10]
dat = dat[!dat$barcode %in% exclude_barcode,]

write.table(dat, paste0(work_path, "/tape_barcode_2/debris_sample_barcode.txt"), row.names=F, col.names=F, sep="\t", quote=F)



##################################
### Mapping cell and sample by the counts of BCs shared between them

dat = read.table(paste0(work_path, "/tape_barcode_2/cell_debris_well_1_mismatch.txt"), sep="\t", header=T)
rownames(dat) = as.vector(dat$cell)
dat = dat[,-1]
dat = as.matrix(dat)

cell_barcode = read.table(paste0(work_path, "/tape_barcode_2/tape_barcode.txt"), sep="\t", header=F)
names(cell_barcode) = c("cell", "bc", "edit", "num")
cell_barcode_num = cell_barcode %>% select(cell, bc) %>% unique() %>% group_by(cell) %>% tally()

### 2d histogram of number of TapeBC recovered per cell vs. maximum TapeBC matching one gastruloid (debris-seq)

max_count = apply(dat, 1, max)
n = ncol(dat)
second_count = apply(dat, 1, function(x) sort(x,partial=n-1)[n-1])

df = data.frame(cell = rownames(dat),
                max_count = max_count,
                ratio_top_second = (max_count+1)/(second_count+1)) %>% 
    left_join(cell_barcode_num, by = "cell") %>% rename(BC_num = n)

df$max_count_log2 = log2(df$max_count + 1)
df$BC_num_log2 = log2(df$BC_num + 1)

library(hexbin)
pdf("~/share/2D_dist_BC_num_max_count.pdf")
hexbinplot(max_count_log2~BC_num_log2, data=df, colramp=colorRampPalette(hcl.colors(12)))
dev.off()

library(hexbin)
pdf("~/share/2D_dist_BC_num_ratio_top_second.pdf")
hexbinplot(ratio_top_second~BC_num_log2, data=df, colramp=colorRampPalette(hcl.colors(12)))
dev.off()

saveRDS(df, paste0(work_path, "/tape_barcode_2/cell_filter_table.RDS"))

print(sum(df$ratio_top_second > 1 & df$BC_num >= 3)/nrow(df))
# 67.8%

### output gastruloid labels for those "confident" cells, used for semi-supervised UMAP

df_con = subset(df, ratio_top_second > 1.5 & BC_num >= 5) ### 48%
dat_con = dat[rownames(dat) %in% as.vector(df_con$cell),]

res = data.frame(cell = rownames(dat_con),
                 sample_id = apply(dat_con, 1, which.max))

cell_list = read.table(paste0(work_path, "/tape_barcode_2/cell_cell_list.txt"), as.is=T)
colnames(cell_list) = "cell"
cell_list = cell_list %>% left_join(res, by = "cell")
cell_list$sample_id[is.na(cell_list$sample_id)] = -1

write.csv(cell_list, paste0(work_path, "/tape_barcode_2/cell_cell_list_label.csv"), row.names = F)


### output gastruloid labels for those "confident" cells, using slightly looser cutoffs

df_con = subset(df, ratio_top_second > 1 & BC_num >= 3) ### 67.8%
dat_con = dat[rownames(dat) %in% as.vector(df_con$cell),]

res = data.frame(cell = rownames(dat_con),
                 sample_id = apply(dat_con, 1, which.max))

cell_list = read.table(paste0(work_path, "/tape_barcode_2/cell_cell_list.txt"), as.is=T)
colnames(cell_list) = "cell"
cell_list = cell_list %>% left_join(res, by = "cell")
cell_list$sample_id[is.na(cell_list$sample_id)] = -1

write.csv(cell_list, paste0(work_path, "/tape_barcode_2/cell_cell_list_label_2.csv"), row.names = F)




##################################
### I am not sure if the semi-supervised UMAP idea works or not, so I will try a different way -
### simply calculating the mean distance between each cell w/o assigned and cells assigned for individual sample

dat = Matrix::readMM(paste0(work_path, "/tape_barcode_2/cell_pairs.mtx"))
cell_list = read.csv(paste0(work_path, "/tape_barcode_2/cell_cell_list_label_2.csv")) ### 67.8% has been assigned

df = readRDS(paste0(work_path, "/tape_barcode_2/cell_filter_table.RDS"))
df_sub = df[df$BC_num >= 3,]

dat = dat + t(dat)
rownames(dat) = colnames(dat) = as.vector(cell_list$cell)

dat_sub = dat[rownames(dat) %in% as.vector(cell_list$cell[cell_list$sample_id == -1]),
              colnames(dat) %in% as.vector(cell_list$cell[cell_list$sample_id != -1])]
dat_sub = dat_sub[rownames(dat_sub) %in% as.vector(df_sub$cell),]

dat_sub_dis = dat_sub
dat_sub_dis@x  = 1/(dat_sub@x + 1) - 1

rownames(cell_list) = as.vector(cell_list$cell)
cell_list = cell_list[colnames(dat_sub_dis),]

sample_list = sort(unique(cell_list$sample_id))

res = NULL
for(sample_i in sample_list){
    print(sample_i)
    res = cbind(res, Matrix::rowMeans(dat_sub_dis[, cell_list$sample_id == sample_i, drop=FALSE]))
}

res_row_min = apply(res, 1, which.min)

max_count = apply(res, 1, min)
second_count = apply(res, 1, function(x) x[order(x)[2]])

res_row_min = res_row_min[max_count != second_count]

cell_list_2 = data.frame(cell = names(res_row_min),
                         sample_id = sample_list[as.vector(res_row_min)],
                         assigning = 2)

cell_list$assigning = 1

cell_list = rbind(cell_list, cell_list_2)

dat = read.table(paste0(work_path, "/tape_barcode_2/cell_debris_well_1_mismatch.txt"), sep="\t", header=T)
dat_x = names(dat)[2:ncol(dat)]

cell_list$sample = dat_x[as.vector(cell_list$sample_id)]

saveRDS(cell_list, paste0(work_path, "/tape_barcode_2/cell_assigning_res.rds"))

### 247,064 cells recaptured from sci-RNA-seq3




