
#############################
### Creating single-cell tree
### Chengxiang Qiu
### Feb-20, 2025

source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

work_path = "Your_work_path"

### The cell x DTT matrix (after the below Step-1 filtering) can be downloaded directly from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/sam_tape
### loci.csv

#################################################################
### Step-1: Filtering some tapes which have too many integrations

dat = read.csv("./loci_wells_allLanes.csv", header=T, row.names=1)
tape = data.frame(tape_site = colnames(dat),
                  integration = paste0(unlist(lapply(colnames(dat), function(x) strsplit(x,"[.]")[[1]][1])),".",
                                       unlist(lapply(colnames(dat), function(x) strsplit(x,"[.]")[[1]][2]))),
                  tape_bc = unlist(lapply(colnames(dat), function(x) strsplit(x,"[.]")[[1]][1])))

### Filtering-1: checking some tape_bc outliers (same barcode show up many times)

tape_freq = tape %>% group_by(tape_bc) %>% tally() %>% mutate(freq = n/6)

### table(tape_freq$freq) ### clone05
### 1   2   3   4   5   6   7 231
### 7  13  26  23   3   1   1   1

keep = tape$tape_bc != "GTTTCAAAATAA"
dat_filter = dat[,keep]
tape = tape[keep,]
### 65465 cells x 1386 tape_sites (231 integrations; 74 BCs)

### Filtering-2: low informative integrations

N = ncol(dat_filter)/6
edit_list = NULL
cell_num_list = NULL
for(i in 1:N){
    x = dat_filter[,c(((i-1)*6+1):(i*6))]
    x = paste(x[,1], x[,2], x[,3], x[,4], x[,5], x[,6], sep="_")
    x = x[x != "None_None_None_None_None_None"]
    y = x[x != "ETY_None_None_None_None_None"]
    edit_list = c(edit_list, length(unique(y)))
    cell_num_list = c(cell_num_list, length(x))
}
df_integration = data.frame(integration = gsub(".Site1", "", colnames(dat_filter)[grep("Site1", colnames(dat_filter))]),
                            edit_num = edit_list,
                            cell_num = cell_num_list)

p = df_integration %>% ggplot() + geom_point(aes(x = edit_num, y = cell_num)) + scale_y_log10() +
    geom_hline(yintercept = 10, linetype = "longdash", color = "red", size = 1.5) + 
    labs(x = "# of distinct editing patterns", y = "# of cells detected the Tape") +
    theme_classic(base_size = 12) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
ggsave(paste0("~/share/Tape_quality.pdf"), p, height = 5, width = 5)

integration_include = as.vector(df_integration$integration[df_integration$cell_num > 10])
tape_include = as.vector(tape$tape_site[tape$integration %in% integration_include])
dat_filter = dat[,colnames(dat) %in% tape_include]
tape_filter = tape[tape$tape_site %in% tape_include,] ### 231 integrations; 74 distinct tape BCs

### Filtering-3: low informtive cells

dat_filter_site1 = dat_filter[,colnames(dat_filter)[grep("Site1", colnames(dat_filter))]]
df_site1 = melt(as.matrix(dat_filter_site1)) %>% filter(value != "None") %>% group_by(Var1) %>% tally()
cutoff_low = mean(df_site1$n) - 1.5*sd(df_site1$n)
cutoff_high = mean(df_site1$n) + 2*sd(df_site1$n)
p = df_site1 %>% ggplot() + geom_histogram(aes(n), bins = 30) +
    geom_vline(xintercept = cutoff_low, linetype = "longdash", color = "red", size = 1.5) + 
    geom_vline(xintercept = cutoff_high, linetype = "longdash", color = "red", size = 1.5) + 
    labs(x = "# of Tapes observed per cell", y = "Count") +
    theme_classic(base_size = 12) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
ggsave(paste0("~/share/Cell_quality.pdf"), p, height = 5, width = 5)

cell_include = as.vector(df_site1$Var1[df_site1$n >= cutoff_low & df_site1$n <= cutoff_high])
dat_filter = dat_filter[cell_include,]

### output the filtered matrix
write.csv(dat_filter, paste0(work_path, "/tree_sanjay/loci.csv"), quote=F)
### 58,283 cells x 1,386 tape_sites (231 tape integrations)

df = melt(as.matrix(dat_filter)) %>% filter(value != "None", value != "ETY") %>%
    group_by(Var1) %>% tally()
print(paste0(round(mean(df$n), 2), " +/- ", round(sd(df$n), 2))) ### 213.73 +/- 89.26
p = df %>% ggplot() + geom_histogram(aes(n), bins = 50) +
    labs(x = "# of edited sites per cell", y = "Count") +
    theme_classic(base_size = 12) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
ggsave(paste0("~/share/Edit_sites.pdf"), p, height = 5, width = 5)



######################################################################
### Step-2: Calculating cell cell distance and reconstructing the tree

### Calculating cell cell distance
### python Calculate_cell_cell_distance.py

### save the dis and num matrix to "*rds" format to save space
x = read.csv(paste0(work_path, "/tree_sanjay/cell_cell_dis.csv"), header=F)
y = read.csv(paste0(work_path, "/tree_sanjay/cell_cell_num.csv"), header=F)
saveRDS(x, paste0(work_path, "/tree_sanjay/cell_cell_dis.rds"))
saveRDS(y, paste0(work_path, "/tree_sanjay/cell_cell_num.rds"))

dat = as.matrix(x/y)
dat[is.nan(dat)] = 0

x = read.csv(paste0(work_path, "/tree_sanjay/cell_list.csv"), header=F)
rownames(dat) = colnames(dat) = as.vector(x$V1)

tree = as.phylo(hclust(as.dist(dat), "average"))
saveRDS(tree, paste0(work_path, "/tree_sanjay/tree.rds"))

p_rec = ggtree(tree)
saveRDS(p_rec, paste0(work_path, "/tree_sanjay/tree_plot_rec.rds"))
ggsave(paste0("~/share/tree_mGASv5.png"), p_rec, dpi = 300)

p_cir = ggtree(tree, layout="fan", size=0.15, open.angle=5)
saveRDS(p_cir, paste0(work_path, "/tree_sanjay/tree_plot_cir.rds"))
ggsave(paste0("~/share/tree_mGASv5_cir.png"), p_cir, dpi = 300)

pd = read.table(paste0(work_path, "/cell_id.txt"), as.is=T, sep="\t")
colnames(pd) = c("cell_id","cell","celltype"); pd$cell = NULL
rownames(pd) = as.vector(pd$cell_id)

sample_group_color = c("#f26522", "#ec008c", "#2e3192", "#00aeef", "#00a651", "#a8b439", "#ed1c24", "grey80")
names(sample_group_color) = c(paste0("group_", c(1:7)), "not_assigned")

### Of note, here I colored cells by their assigned wells/gastruloids, you could find how this step was performed from the below section
assign = readRDS(paste0(work_path, "/tree_sanjay/assignment/assignment_cell.rds"))
colors_96 <- c(
    "#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF", "#800000", "#808000",
    "#008000", "#800080", "#008080", "#000080", "#FFA500", "#FFC0CB", "#A52A2A", "#D2691E",
    "#5F9EA0", "#7FFF00", "#DC143C", "#006400", "#8B0000", "#556B2F", "#FF8C00", "#9932CC",
    "#8FBC8F", "#483D8B", "#2F4F4F", "#9400D3", "#FF1493", "#00BFFF", "#696969", "#1E90FF",
    "#B22222", "#228B22", "#FF69B4", "#CD5C5C", "#4B0082", "#20B2AA", "#778899", "#32CD32",
    "#00FF7F", "#7CFC00", "#B0C4DE", "#ADFF2F", "#FF4500", "#DA70D6", "#F4A460", "#BA55D3",
    "#FFD700", "#7B68EE", "#48D1CC", "#C71585", "#191970", "#0000CD", "#8B4513", "#FA8072",
    "#E9967A", "#9400D3", "#9370DB", "#3CB371", "#7B68EE", "#00FA9A", "#48D1CC", "#C71585",
    "#FF6347", "#4682B4", "#D8BFD8", "#008B8B", "#B8860B", "#DAA520", "#BDB76B", "#A9A9A9",
    "#F08080", "#FFDAB9", "#EEE8AA", "#98FB98", "#6B8E23", "#DB7093", "#CD853F", "#FFB6C1",
    "#8A2BE2", "#00CED1", "#1E90FF", "#FF4500", "#FF8C00", "#FFD700", "#ADFF2F", "#32CD32",
    "#8B008B", "#800080", "#FF00FF", "#DC143C", "#8B0000", "#A52A2A", "#B22222", "#D2691E")
names(colors_96) = unique(assign$well)

clade_node = data.frame(group = paste0("group_", c(1:7)),
                        node = paste0("node_", c(5,4,6,32,160,408,407)))
dat_top1000 = readRDS(paste0(work_path, "/tree_sanjay/ancestor_matrix_top1000.rds"))
df = data.frame(cell_id = as.vector(tree$tip.label), group = "not_assigned")
for(i in 1:7){
    cell_include = rownames(dat_top1000)[dat_top1000[,clade_node$node[i]] == 1]
    df$group[df$cell_id %in% cell_include] = clade_node$group[i]
}

dat_label = pd[as.vector(tree$tip.label),]
dat_label$group = as.vector(df$group)
dat_label = dat_label %>% left_join(assign %>% select(cell_id = cell, well), by = "cell_id")
dat_label$well[is.na(dat_label$well)] = 'not_assigned'

p_x <- p_cir +
    geom_fruit(data=dat_label, geom=geom_tile,
               mapping=aes(y=cell_id, fill=group),
               offset = 0.07, pwidth = 0.5) +
    geom_fruit(data=dat_label, geom=geom_tile,
               mapping=aes(y=cell_id, fill=well),
               offset = 0.14, pwidth = 0.5) +
    geom_fruit(data=dat_label, geom=geom_tile,
               mapping=aes(y=cell_id, fill=celltype),
               offset = 0.14, pwidth = 0.5) +
    scale_fill_manual(values = c(gastruloid_celltype_color_code, sample_group_color, colors_96)) +
    theme(legend.position="none")
ggsave(paste0("~/share/tree_mGASv5_cir_celltype.png"), p_x, dpi = 300)


### label the seven clades:

tree = readRDS(paste0(work_path, "/tree_sanjay/tree.rds"))
p_cir = readRDS(paste0(work_path, "/tree_sanjay/tree_plot_cir.rds"))
select_node = c(5,4,6,32,160,408,407) + 58283
p_x <- p_cir +
    geom_label2(aes(subset = node %in% select_node, label = node - 58283), size = 1) +
    theme(legend.position="none")
ggsave(paste0("~/share/tree_mGASv5_cir_seven_nodes.png"), p_x, dpi = 300)


### output the ancestor-children matrix
dat = extract_ancestor_nodes_fast(tree)
saveRDS(dat, paste0(work_path, "/tree_sanjay/ancestor_matrix.rds"))

sparse_mat = Matrix(dat, sparse = TRUE)
writeMM(sparse_mat, paste0(work_path, "/tree_sanjay/ancestor_matrix.mtx"))
write.csv(rownames(sparse_mat), paste0(work_path, "/tree_sanjay/ancestor_matrix.row.csv"))
write.csv(colnames(sparse_mat), paste0(work_path, "/tree_sanjay/ancestor_matrix.col.csv"))

### also output the matrix with top 1000 ancestor nodes (ordering by depth)
node_depths = data.frame(node = 1:ncol(dat), 
                         depth = node.depth.edgelength(tree)[(length(tree$tip.label)+1):(length(tree$tip.label)+tree$Nnode)])
node_depths_top1000 = node_depths[node_depths$depth < quantile(node_depths$depth, 1000/nrow(node_depths)),]

dat_top1000 = dat[,as.numeric(node_depths_top1000$node)]
saveRDS(dat_top1000, paste0(work_path, "/tree_sanjay/ancestor_matrix_top1000.rds"))

sparse_mat = Matrix(dat_top1000, sparse = TRUE)
writeMM(sparse_mat, paste0(work_path, "/tree_sanjay/ancestor_matrix_top1000.mtx"))



#############################################################
### Step-3: for bootstrapping, save them in RDS to save space

work_path = "PATH_1"
save_path = "PATH_2"

args = commandArgs(trailingOnly=TRUE)
kk = as.numeric(args[1])
print(paste0("Processing:", kk))

x = read.csv(paste0(save_path, "/tree_sanjay/bootstrap/cell_cell_dis_", kk, ".csv"), header=F)
y = read.csv(paste0(save_path, "/tree_sanjay/bootstrap/cell_cell_num_", kk, ".csv"), header=F)

dat = as.matrix(x/y)
dat[is.nan(dat)] = 0

rm(x, y)
gc()

x = read.csv(paste0(work_path, "/tree_sanjay/cell_list.csv"), header=F)
rownames(dat) = colnames(dat) = as.vector(x$V1)

tree = as.phylo(hclust(as.dist(dat), "average"))
saveRDS(tree, paste0(save_path, "/tree_sanjay/bootstrap/tree_bootstrap_", kk, ".rds"))

tree = readRDS(paste0(save_path, "/tree_sanjay/bootstrap/tree_bootstrap_", kk, ".rds"))
dat_anc = extract_ancestor_nodes_fast(tree)
saveRDS(dat_anc, paste0(save_path, "/tree_sanjay/bootstrap/ancestor_matrix_", kk, ".rds"))

sparse_mat = Matrix(dat_anc, sparse = TRUE)
writeMM(sparse_mat, paste0(save_path, "/tree_sanjay/bootstrap/ancestor_matrix_", kk, ".mtx"))


####################################################################
### Step-4: Performing TBE to investigate the robustness of the tree

### python Calculate_TBE.py

tree_orig = readRDS(paste0(work_path, "/tree_sanjay/tree.rds"))
dat_top1000 = readRDS(paste0(work_path, "/tree_sanjay/ancestor_matrix_top1000.rds"))
N_cell = nrow(dat_top1000)

N_p = NULL
for(i in 1:ncol(dat_top1000)){
    N_p = c(N_p, min(sum(dat_top1000[,i] == 1), sum(dat_top1000[,i] == 0)))
}

res = NULL
for(kk in 1:100){
    print(kk)
    res_kk = read.csv(paste0(save_path, "/tree_sanjay/bootstrap/TBE_", kk, ".csv"))
    res = cbind(res, as.vector(res_kk$score))
}

score = apply(as.matrix(res), 1, mean)

x2 = NULL
for(i in 1:nrow(res)){
    if(score[i] == 0){
        x2 = rbind(x2, data.frame(node = colnames(dat_top1000)[i], TBE = 1 - score[i]))
    } else {
        x2 = rbind(x2, data.frame(node = colnames(dat_top1000)[i], TBE = 1 - score[i]/(N_p[i] - 1)))
    }
}
saveRDS(x2, paste0(work_path, "/tree_sanjay/ancestor_top1000_TBE.rds"))

print(sum(x2$TBE >= 0.7)/nrow(x2))
### 45%

node_depths = data.frame(node = paste0("node_", 1:tree_orig$Nnode), 
                         depth = node.depth.edgelength(tree_orig)[(length(tree_orig$tip.label)+1):(length(tree_orig$tip.label)+tree_orig$Nnode)])
x2 = x2 %>% left_join(node_depths, by = "node")

cell_num = data.frame(cell_num = apply(dat_top1000, 2, sum), 
                      node = colnames(dat_top1000))
cell_num$log2_cell_num = log2(cell_num$cell_num)
x2 = x2 %>% left_join(cell_num, by = "node")

p = x2 %>% ggplot() + geom_point(aes(x = TBE, y = depth)) +
    geom_vline(xintercept = 0.7, color = "red", size = 1.5) + 
    labs(x = "Transfer bootstrap expectation (TBE)", y = "Distance from root", title = "Top 1000 ancestor nodes ranked by depth") +
    theme_classic(base_size = 12) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
ggsave(paste0("~/share/scatter_TBE_depth.pdf"), p, width = 5, height = 5)

p = x2 %>% ggplot() + geom_point(aes(x = TBE, y = log2_cell_num)) +
    geom_vline(xintercept = 0.7, color = "red", size = 1.5) + 
    labs(x = "Transfer bootstrap expectation (TBE)", y = "Log2 cell number", title = "Top 1000 ancestor nodes ranked by depth") +
    theme_classic(base_size = 12) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
ggsave(paste0("~/share/scatter_TBE_cell_num.pdf"), p, width = 5, height = 5)



############################################################################
### Step-5: assigning clades in the single-cell with potential founder cells

### debris-seq data
well = readRDS(paste0(work_path, "/tree/well_tape_dominant_edits.rds"))
well$total = NULL
well$combination = paste0(well$TargetBC, "-", well$Site, "-", well$edit)
well = well %>% group_by(Well, combination) %>% tally()

### single-cell data
dat = read.csv(paste0(work_path, "/tree_sanjay/loci.csv"), row.names=1, header=T)
tape = data.frame(tape_site = colnames(dat),
                  integration = paste0(unlist(lapply(colnames(dat), function(x) strsplit(x,"[.]")[[1]][1])),".",
                                       unlist(lapply(colnames(dat), function(x) strsplit(x,"[.]")[[1]][2]))),
                  TargetBC = unlist(lapply(colnames(dat), function(x) strsplit(x,"[.]")[[1]][1])),
                  Site = unlist(lapply(colnames(dat), function(x) strsplit(x,"[.]")[[1]][3])))
dat_x = melt(as.matrix(dat))
dat_x = dat_x %>% filter(value != 'None', value != 'ETY')
dat_x = dat_x %>% rename(cell=Var1, tape_site=Var2, edit=value) %>% 
    left_join(tape[,c("TargetBC","Site","tape_site")], by = "tape_site")
dat_x$combination = paste0(dat_x$TargetBC, "-", dat_x$Site, "-", dat_x$edit)
dat_x = dat_x %>% group_by(cell, combination) %>% tally()

### combining two data
overlap_combination = unique(well$combination[well$combination %in% dat_x$combination])

well = well[well$combination %in% overlap_combination,]
well$n = 1
well_x = dcast(data = well, Well~combination, fill = 0)
rownames(well_x) = as.vector(well_x[,1])
well_x = as.matrix(well_x[,-1])
saveRDS(well_x, paste0(work_path, "/tree_sanjay/assignment/well_x_edit.rds"))

dat_x = dat_x[dat_x$combination %in% overlap_combination,]
dat_x$n = 1
dat_x = dcast(data = dat_x, cell~combination, fill = 0)
rownames(dat_x) = as.vector(dat_x[,1])
dat_x = as.matrix(dat_x[,-1])
saveRDS(dat_x, paste0(work_path, "/tree_sanjay/assignment/cell_x_edit.rds"))


######################################################################
### If we focus on the seven major clades in the debris-seq based tree

source("~/work/scripts/utils.R")
work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape/data_mGASv5"

well_x = readRDS(paste0(work_path, "/tree_sanjay/assignment/well_x_edit.rds"))
well_x = Matrix(well_x, sparse = TRUE)
dat_x  = readRDS(paste0(work_path, "/tree_sanjay/assignment/cell_x_edit.rds"))
dat_x = Matrix(dat_x, sparse = TRUE)

### clades cell data
dat = readRDS(paste0(work_path, "/tree_sanjay/ancestor_matrix.sparseMatrix.rds"))
print(sum(colnames(well_x) == colnames(dat_x)))
print(sum(rownames(dat) == rownames(dat_x)))

sample_group = readRDS(paste0(work_path, "/tree/sample_group.rds"))
sample_group = dcast(sample_group %>% group_by(sample_id, sample_group) %>% tally(), sample_group~sample_id, fill = 0)
rownames(sample_group) = sample_group[,1]
sample_group = as.matrix(sample_group[,-1])

print(sum(colnames(sample_group) == rownames(well_x)))
debris_clade_x = sample_group %*% well_x
debris_clade_x[debris_clade_x != 0] = 1

keep = apply(debris_clade_x, 2, sum) == 1
debris_clade_x_exc = debris_clade_x[,keep]
debris_clade_combin_num = apply(debris_clade_x_exc, 1, sum)
print(table(debris_clade_combin_num))

dat_x_norm = t(t(dat_x)/colSums(dat_x))
clade_x_norm = t(dat) %*% dat_x_norm

res = NULL
debris_clade_list = NULL
for(debris_clade_i in 1:nrow(debris_clade_x_exc)){
    print(debris_clade_i)
    combination_yes = colnames(debris_clade_x_exc)[debris_clade_x_exc[debris_clade_i,] == 1]
    combination_no  = colnames(debris_clade_x_exc)[debris_clade_x_exc[debris_clade_i,] == 0]
    if(length(combination_yes) > 0){
        dat_yes = rowMeans(clade_x_norm[,colnames(clade_x_norm) %in% combination_yes])
        dat_no  = rowMeans(clade_x_norm[,colnames(clade_x_norm) %in% combination_no])
        res = cbind(res, log2(as.vector(dat_yes)/as.vector(dat_no) + 1))
        debris_clade_list = c(debris_clade_list, rownames(debris_clade_x_exc)[debris_clade_i])
    }
}
rownames(res) = rownames(clade_x_norm)
colnames(res) = debris_clade_list
res = data.frame(res)
res$cell_num = as.vector(colSums(dat))


assign = NULL
for(ii in 1:(ncol(res)-1)){
    res_tmp = res[res[,ii] > 2,]
    res_tmp = res_tmp[order(res_tmp$cell_num, decreasing=T),]
    for(jj in 1:nrow(res_tmp)){
        x = res_tmp[jj, ii]
        y = res_tmp[jj, -ii]
        y = max(y[1:(length(y)-1)])
        if(x/y > 2){
            assign = rbind(assign, res_tmp[jj,])
            break
        }
    }
}
saveRDS(assign, paste0(work_path, "/tree_sanjay/assignment/assignment_major.rds"))

library("gplots")
library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
pdf(paste0("~/share/assignment_major.pdf"), 8, 5)
heatmap.2(as.matrix(assign[,-ncol(assign)]), 
          col=Colors, 
          scale="row", 
          Rowv = F, 
          Colv = F, 
          key=T, 
          density.info="none", 
          trace="none", 
          cexRow = 1, 
          cexCol = 1,
          margins = c(5,5))
dev.off()


### assigining ancestor nodes (single-cell tree) to seven major clades (debris-seq tree)
clade_node = data.frame(group = paste0("group_", c(1:7)),
                        node = paste0("node_", c(5,4,6,32,160,408,407)))

pd = read.table(paste0(work_path, "/cell_id.txt"), as.is=T, sep="\t")
colnames(pd) = c("cell_id","cell","celltype"); pd$cell = NULL
rownames(pd) = as.vector(pd$cell_id)
pd = pd[rownames(dat_top1000),]

df = NULL
for(i in 1:nrow(clade_node)){
    pd_x = pd[dat_top1000[,colnames(dat_top1000) == clade_node$node[i]] == 1,]
    pd_x$node = paste0(clade_node$group[i], '-', clade_node$node[i])
    df = rbind(df, pd_x)
}
df = df %>% group_by(celltype, node) %>% tally()
df$celltype = factor(df$celltype, levels = names(gastruloid_celltype_color_code)[names(gastruloid_celltype_color_code) %in% df$celltype])

# Stacked + percent
p = df %>%
    ggplot(aes(fill=celltype, y=n, x=node)) + 
    geom_bar(position="fill", stat="identity", width = 0.8) +
    scale_fill_manual(values=gastruloid_celltype_color_code) +
    labs(x="", y="% of cells", title="") +
    theme_classic(base_size = 15) +
    #theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black", angle = 45), axis.text.y = element_text(color="black"))

ggsave("~/share/seven_clades_node_celltype_frac.pdf", p, height  = 5, width = 8)


##################################################
### can we repeat this for individual major clade?

well_x = readRDS(paste0(work_path, "/tree_sanjay/assignment/well_x_edit.rds"))
well_x = Matrix(well_x, sparse = TRUE)
dat_x  = readRDS(paste0(work_path, "/tree_sanjay/assignment/cell_x_edit.rds"))
dat_x = Matrix(dat_x, sparse = TRUE)
dat = readRDS(paste0(work_path, "/tree_sanjay/ancestor_matrix.sparseMatrix.rds"))
print(sum(colnames(well_x) == colnames(dat_x)))
print(sum(rownames(dat) == rownames(dat_x)))

clade_node = data.frame(group = paste0("group_", c(1:7)),
                        node = paste0("node_", c(5,4,6,32,160,408,407)))

sample_group = readRDS(paste0(work_path, "/tree/sample_group.rds"))

df_all = NULL
assign_x_all = NULL

library("gplots")
library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
exclude_well = c("P1.B10","P2.C9",
                 "P1.A5","P1.A4","P2.C3","P2.D12",
                 "P1.D9","P1.D8",
                 "P1.D6","P2.A12","P2.D8")

for(cnt in 1:7){

    well_include = as.vector(sample_group$sample_id[sample_group$sample_group == paste0("group_", cnt)])
    well_x_sub = well_x[well_include,]
    
    cell_include = dat[,colnames(dat) == clade_node$node[cnt]] == 1
    dat_x_sub = dat_x[cell_include,]
    dat_x_sub = dat_x_sub[,colSums(dat_x_sub) != 0]
    dat_tmp = dat[!cell_include,]
    dat_sub = dat[cell_include, colSums(dat_tmp) == 0]
    
    well_x_sub = well_x_sub[,colnames(well_x_sub) %in% colnames(dat_x_sub)]
    keep = apply(well_x_sub, 2, sum) == 1
    well_x_sub_exc = well_x_sub[,keep]
    debris_clade_combin_num = apply(well_x_sub_exc, 1, sum)
    print(table(debris_clade_combin_num))
    
    dat_x_norm = t(t(dat_x_sub)/colSums(dat_x_sub))
    clade_x_norm = t(dat_sub) %*% dat_x_norm
    
    res = NULL
    debris_clade_list = NULL
    for(debris_clade_i in 1:nrow(well_x_sub_exc)){
        print(debris_clade_i)
        combination_yes = colnames(well_x_sub_exc)[well_x_sub_exc[debris_clade_i,] == 1]
        combination_no  = colnames(well_x_sub_exc)[well_x_sub_exc[debris_clade_i,] == 0]
        if(length(combination_yes) > 0){
            dat_yes = rowMeans(clade_x_norm[,colnames(clade_x_norm) %in% combination_yes,drop=F])
            dat_no  = rowMeans(clade_x_norm[,colnames(clade_x_norm) %in% combination_no,drop=F])
            res = cbind(res, log2(as.vector(dat_yes)/as.vector(dat_no) + 1))
            debris_clade_list = c(debris_clade_list, rownames(well_x_sub_exc)[debris_clade_i])
        }
    }
    rownames(res) = rownames(clade_x_norm)
    colnames(res) = debris_clade_list
    res = data.frame(res)
    res$cell_num = as.vector(colSums(dat_sub))
    
    assign = NULL
    for(ii in 1:(ncol(res)-1)){
        if(!colnames(res)[ii] %in% exclude_well){
            res_tmp = res[res[,ii] > 2,]
            res_tmp = res_tmp[order(res_tmp$cell_num, decreasing=T),]
            for(jj in 1:nrow(res_tmp)){
                x = res_tmp[jj, ii]
                y = res_tmp[jj, -ii]
                y = max(y[1:(length(y)-1)])
                if(colnames(res)[ii] %in% c('P1.C4')) {
                    cutoff = 3
                } else {
                    cutoff = 2
                }
                if(x/y > cutoff){
                    assign = rbind(assign, res_tmp[jj,])
                    break
                }
            }
        }
    }
    assign = assign[,!colnames(assign) %in% exclude_well]
    
    pdf(paste0("~/share/assignment_group_", cnt, ".pdf"), 8, 5)
    heatmap.2(as.matrix(assign[,-ncol(assign)]), 
              col=Colors, 
              scale="row", 
              Rowv = F, 
              Colv = F, 
              key=T, 
              density.info="none", 
              trace="none", 
              cexRow = 1, 
              cexCol = 1,
              margins = c(5,5))
    dev.off()
    
    assign_x = data.frame(well = gsub('[.]','-',colnames(assign)[-ncol(assign)]),
                          node = rownames(assign),
                          group = paste0("group_", cnt),
                          cell_num = as.vector(assign$cell_num))
    assign_x_all = rbind(assign_x_all, assign_x)
    
    df = NULL
    for(ii in 1:nrow(assign)){
        df = rbind(df, data.frame(cell = rownames(dat)[dat[,rownames(assign)[ii]] == 1],
                                  node = rownames(assign)[ii],
                                  well = gsub('[.]','-',colnames(assign)[ii])))
    }
    df_all = rbind(df_all, df)
}

### 107 well (P2-B10 missing distinct combinations)
exclude_well = c("P1-B10","P2-C9",
                 "P1-A5","P1-A4","P2-C3","P2-D12",
                 "P1-D9","P1-D8",
                 "P1-D6","P2-A12","P2-D8")
### 11 wells are excluded due to ambiguous assigned
assign_x_all = assign_x_all[!assign_x_all$well %in% exclude_well,]
df_all = df_all[!df_all$well %in% exclude_well,]

### 88 with >= 50 cells

saveRDS(assign_x_all, paste0(work_path, "/tree_sanjay/assignment/assignment_well.rds"))
saveRDS(df_all, paste0(work_path, "/tree_sanjay/assignment/assignment_cell.rds"))

dat = readRDS(paste0(work_path, "/tree_sanjay/ancestor_matrix.sparseMatrix.rds"))
dat_sub = dat[,colnames(dat) %in% as.vector(assign_x_all$node)]
saveRDS(dat_sub, paste0(work_path, "/tree_sanjay/ancestor_matrix_sub.sparseMatrix.rds"))
writeMM(dat_sub, paste0(work_path, "/tree_sanjay/ancestor_matrix_sub.sparseMatrix.mtx"))


### plot cell number of each well, and its correlation with gastruloid size
assign_x_all = readRDS(paste0(work_path, "/tree_sanjay/assignment/assignment_well.rds"))
gas_size = read.table(paste0(work_path, "/tree/gastruloid_size.txt"))
colnames(gas_size) = c("well", "size")

df = gas_size %>% left_join(assign_x_all[,c("well","cell_num")], by = "well") %>% as.data.frame()
df$cell_num[is.na(df$cell_num)] = 0
df$log2_cell_num = log2(df$cell_num + 1)
df$log2_size = log2(df$size + 1)

p1 = df %>% ggplot() +
    geom_histogram(aes(log2_cell_num), bins = 20) + geom_vline(xintercept = log2(51), color = "red") + plot_tmp
p2 = df %>% ggplot() +
    geom_point(aes(x = log2_size, y = log2_cell_num)) + plot_tmp
ggsave("~/share/gastruloid_cell_num_size.pdf", p1+p2, width = 8, height = 4)

cor.test(df$log2_size, df$log2_cell_num)$estimate ### 0.7198781
cor.test(df$log2_size, df$log2_cell_num)$p.value ### 1.652389e-18

summary(df$cell_num)
Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
0.0    74.5   265.0   481.1   650.2  5319.0

sd(df$cell_num)
684.7488

### if we only consider 96 wells assigned
df_sub = df[df$cell_num != 0,]

summary(df_sub$cell_num)
Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
2.0   134.8   365.0   541.3   748.5  5319.0

sd(df_sub$cell_num)
703.7126


#############################################################
### Step-6: obtaining TBE for those nodes assigned with wells

tree_orig = readRDS(paste0(work_path, "/tree_sanjay/tree.rds"))
dat_sub = readRDS(paste0(work_path, "/tree_sanjay/ancestor_matrix_sub.sparseMatrix.rds"))
N_cell = nrow(dat_sub)

N_p = NULL
for(i in 1:ncol(dat_sub)){
    N_p = c(N_p, min(sum(dat_sub[,i] == 1), sum(dat_sub[,i] == 0)))
}

res = NULL
for(kk in 1:100){
    print(kk)
    res_kk = read.csv(paste0(save_path, "/tree_sanjay/bootstrap/TBE_", kk, ".csv"))
    res = cbind(res, as.vector(res_kk$score))
}

score = apply(as.matrix(res), 1, mean)

x2 = NULL
for(i in 1:nrow(res)){
    if(score[i] == 0){
        x2 = rbind(x2, data.frame(node = colnames(dat_sub)[i], TBE = 1 - score[i]))
    } else {
        x2 = rbind(x2, data.frame(node = colnames(dat_sub)[i], TBE = 1 - score[i]/(N_p[i] - 1)))
    }
}

node_depths = data.frame(node = paste0("node_", 1:tree_orig$Nnode), 
                         depth = node.depth.edgelength(tree_orig)[(length(tree_orig$tip.label)+1):(length(tree_orig$tip.label)+tree_orig$Nnode)])
x2 = x2 %>% left_join(node_depths, by = "node")

assign_x_all = readRDS(paste0(work_path, "/tree_sanjay/assignment/assignment_well.rds"))
assign_x_all = assign_x_all %>% left_join(x2, by = "node")

node_id = (tree$Nnode + 1) + as.numeric(gsub("node_","",as.vector(assign_x_all$node)))
division_time = NULL
for(i in 1:length(node_id)){
    division_time = c(division_time, how_many_division(tree_orig, node_id[i]))
}
assign_x_all$division_round = division_time
saveRDS(assign_x_all, paste0(work_path, "/tree_sanjay/ancestor_sub_TBE.rds"))
write.csv(assign_x_all, paste0(work_path, "/tree_sanjay/ancestor_sub_TBE.csv"), row.names=F)






