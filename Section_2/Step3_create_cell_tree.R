
###################################################
### Reconstruct the phylo tree for individual clone
### Chengxiang Qiu
### Feb-20, 2025

#################################################################
### Step-1: Filtering some tapes which have too many integrations

source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

work_path = "Your_work_path"

### scRNA-seq DTT data (cell x DTT matrix format, constructed by first splitting 
### some TapeBCs to account for piggyBAC excision and re-integration events) can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/sam_tape
### 1) clone05: loci_clone05.csv
### 2) clone25: loci_clone25.csv
### 3) clone32: loci_clone32.csv


###############
### clone05 ###
###############

clone_i = "clone05"

dat = read.csv(paste0(work_path, "/loci_", clone_i, ".csv"), header=T, row.names=1)
tape = data.frame(tape_site = colnames(dat),
                  integration = paste0(unlist(lapply(colnames(dat), function(x) strsplit(x,"[.]")[[1]][1])),".",
                                       unlist(lapply(colnames(dat), function(x) strsplit(x,"[.]")[[1]][2]))),
                  tape_bc = unlist(lapply(colnames(dat), function(x) strsplit(x,"[.]")[[1]][1])))


### Filtering-1: checking some tape_bc outliers (same barcode show up many times)

tape_freq = tape %>% group_by(tape_bc) %>% tally() %>% mutate(freq = n/6)

### table(tape_freq$freq) ### clone05
### 1  2  3  4
### 8 48 13  1

dat_filter = dat
### 2678 cells x 882 tape_sites (147 integrations; 70 BCs)

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
ggsave(paste0("~/share/Tape_quality_", clone_i, ".pdf"), p, height = 5, width = 5)

integration_include = as.vector(df_integration$integration[df_integration$cell_num > 10])
tape_include = as.vector(tape$tape_site[tape$integration %in% integration_include])
dat_filter = dat[,colnames(dat) %in% tape_include]
tape_filter = tape[tape$tape_site %in% tape_include,] ### 146 integrations; 70 distinct tape BCs

### Filtering-3: low informtive cells

dat_filter_site1 = dat_filter[,colnames(dat_filter)[grep("Site1", colnames(dat_filter))]]
df_site1 = melt(as.matrix(dat_filter_site1)) %>% filter(value != "None") %>% group_by(Var1) %>% tally()
cutoff_low = mean(df_site1$n) - 1.5*sd(df_site1$n)
cutoff_high = mean(df_site1$n) + 2*sd(df_site1$n)
p = df_site1 %>% ggplot() + geom_histogram(aes(n), bins = 50) +
    geom_vline(xintercept = cutoff_low, linetype = "longdash", color = "red", size = 1.5) + 
    geom_vline(xintercept = cutoff_high, linetype = "longdash", color = "red", size = 1.5) + 
    labs(x = "# of Tapes observed per cell", y = "Count") +
    theme_classic(base_size = 12) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
ggsave(paste0("~/share/Cell_quality_", clone_i, ".pdf"), p, height = 5, width = 5)

cell_include = as.vector(df_site1$Var1[df_site1$n >= cutoff_low & df_site1$n <= cutoff_high])
dat_filter = dat_filter[cell_include,]

### output the filtered matrix
write.csv(dat_filter, paste0(work_path, "/tree_sanjay/loci_", clone_i, ".csv"), quote=F)
### 2438 cells x 876 tape_sites (146 tape integrations)

df = melt(as.matrix(dat_filter)) %>% filter(value != "None", value != "ETY") %>%
    group_by(Var1) %>% tally()
print(paste0(round(mean(df$n), 2), " +/- ", round(sd(df$n), 2))) ### 142.70 +/- 31.12
p = df %>% ggplot() + geom_histogram(aes(n), bins = 100) +
    labs(x = "# of edited sites per cell", y = "Count") +
    theme_classic(base_size = 12) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
ggsave(paste0("~/share/Edit_sites_", clone_i, ".pdf"), p, height = 5, width = 5)




###############
### clone25 ###
###############

clone_i = "clone25"

dat = read.csv(paste0(work_path, "/loci_", clone_i, ".csv"), header=T, row.names=1)
tape = data.frame(tape_site = colnames(dat),
                  integration = paste0(unlist(lapply(colnames(dat), function(x) strsplit(x,"[.]")[[1]][1])),".",
                                       unlist(lapply(colnames(dat), function(x) strsplit(x,"[.]")[[1]][2]))),
                  tape_bc = unlist(lapply(colnames(dat), function(x) strsplit(x,"[.]")[[1]][1])))


### Filtering-1: checking some tape_bc outliers (same barcode show up many times)

tape_freq = tape %>% group_by(tape_bc) %>% tally() %>% mutate(freq = n/6)

### table(tape_freq$freq) ### clone25
### 1  2  3  4
### 73 27  4  1

dat_filter = dat
### 4840 cells x 858 tape_sites (143 integrations; 105 BCs)

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
ggsave(paste0("~/share/Tape_quality_", clone_i, ".pdf"), p, height = 5, width = 5)

integration_include = as.vector(df_integration$integration[df_integration$cell_num > 10])
tape_include = as.vector(tape$tape_site[tape$integration %in% integration_include])
dat_filter = dat[,colnames(dat) %in% tape_include]
tape_filter = tape[tape$tape_site %in% tape_include,] ### 140 integrations; 105 distinct tape BCs

### Filtering-3: low informtive cells

dat_filter_site1 = dat_filter[,colnames(dat_filter)[grep("Site1", colnames(dat_filter))]]
df_site1 = melt(as.matrix(dat_filter_site1)) %>% filter(value != "None") %>% group_by(Var1) %>% tally()
cutoff_low = mean(df_site1$n) - 1.5*sd(df_site1$n)
cutoff_high = mean(df_site1$n) + 2*sd(df_site1$n)
p = df_site1 %>% ggplot() + geom_histogram(aes(n), bins = 30) +
    geom_vline(xintercept = cutoff_low, linetype = "longdash", color = "red", size = 1.5) + 
    geom_vline(xintercept = cutoff_high, linetype = "longdash", color = "red", size = 1.5) + 
    labs(x = "# of Tapes per cell", y = "Count") +
    theme_classic(base_size = 12) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
ggsave(paste0("~/share/Cell_quality_", clone_i, ".pdf"), p, height = 5, width = 5)

cell_include = as.vector(df_site1$Var1[df_site1$n >= cutoff_low & df_site1$n <= cutoff_high])
dat_filter = dat_filter[cell_include,]

### output the filtered matrix
write.csv(dat_filter, paste0(work_path, "/tree_sanjay/loci_", clone_i, ".csv"), quote=F)
### 4418 cells x 840 tape_sites (140 tape integrations)

df = melt(as.matrix(dat_filter)) %>% filter(value != "None", value != "ETY") %>%
    group_by(Var1) %>% tally()
print(paste0(round(mean(df$n), 2), " +/- ", round(sd(df$n), 2))) ### 144.72 +/- 51.72
p = df %>% ggplot() + geom_histogram(aes(n), bins = 100) +
    labs(x = "# of edited sites per cell", y = "Count") +
    theme_classic(base_size = 12) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
ggsave(paste0("~/share/Edit_sites_", clone_i, ".pdf"), p, height = 5, width = 5)





###############
### clone32 ###
###############

clone_i = "clone32"

dat = read.csv(paste0(work_path, "/loci_", clone_i, ".csv"), header=T, row.names=1)
tape = data.frame(tape_site = colnames(dat),
                  integration = paste0(unlist(lapply(colnames(dat), function(x) strsplit(x,"[.]")[[1]][1])),".",
                                       unlist(lapply(colnames(dat), function(x) strsplit(x,"[.]")[[1]][2]))),
                  tape_bc = unlist(lapply(colnames(dat), function(x) strsplit(x,"[.]")[[1]][1])))


### Filtering-1: checking some tape_bc outliers (same barcode show up many times)

tape_freq = tape %>% group_by(tape_bc) %>% tally() %>% mutate(freq = n/6)

### table(tape_freq$freq) ### clone32
### 1  2
### 56 11

dat_filter = dat
### 1027 cells x 468 tape_sites (78 integrations; 67 BCs)

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
ggsave(paste0("~/share/Tape_quality_", clone_i, ".pdf"), p, height = 5, width = 5)

integration_include = as.vector(df_integration$integration[df_integration$cell_num > 10])
tape_include = as.vector(tape$tape_site[tape$integration %in% integration_include])
dat_filter = dat[,colnames(dat) %in% tape_include]
tape_filter = tape[tape$tape_site %in% tape_include,] ### 78 integrations; 67 distinct tape BCs

### Filtering-3: low informtive cells

dat_filter_site1 = dat_filter[,colnames(dat_filter)[grep("Site1", colnames(dat_filter))]]
df_site1 = melt(as.matrix(dat_filter_site1)) %>% filter(value != "None") %>% group_by(Var1) %>% tally()
cutoff_low = mean(df_site1$n) - 1.5*sd(df_site1$n)
cutoff_high = mean(df_site1$n) + 2*sd(df_site1$n)
p = df_site1 %>% ggplot() + geom_histogram(aes(n), bins = 30) +
    geom_vline(xintercept = cutoff_low, linetype = "longdash", color = "red", size = 1.5) + 
    geom_vline(xintercept = cutoff_high, linetype = "longdash", color = "red", size = 1.5) + 
    labs(x = "# of Tapes per cell", y = "Count") +
    theme_classic(base_size = 12) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
ggsave(paste0("~/share/Cell_quality_", clone_i, ".pdf"), p, height = 5, width = 5)

cell_include = as.vector(df_site1$Var1[df_site1$n >= cutoff_low & df_site1$n <= cutoff_high])
dat_filter = dat_filter[cell_include,]

### output the filtered matrix
write.csv(dat_filter, paste0(work_path, "/tree_sanjay/loci_", clone_i, ".csv"), quote=F)
### 960 cells x 468 tape_sites (78 tape integrations)

df = melt(as.matrix(dat_filter)) %>% filter(value != "None", value != "ETY") %>%
    group_by(Var1) %>% tally()
print(paste0(round(mean(df$n), 2), " +/- ", round(sd(df$n), 2))) ### 93.23 +/- 30.23
p = df %>% ggplot() + geom_histogram(aes(n), bins = 50) +
    labs(x = "# of edited sites per cell", y = "Count") +
    theme_classic(base_size = 12) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
ggsave(paste0("~/share/Edit_sites_", clone_i, ".pdf"), p, height = 5, width = 5)





##################################################
### Step-2: making and plotting tree on each clone

### python Calculate_cell_cell_distance.py

clone_i = "clone32"

dat_1 = read.csv(paste0(work_path, "/tree_sanjay/cell_cell_orig_distance_", clone_i, ".csv"), header=F)
dat_2 = read.csv(paste0(work_path, "/tree_sanjay/cell_cell_comp_", clone_i, ".csv"), header=F)
dat = as.matrix(dat_1/dat_2)
dat[is.nan(dat)] = 0

x = read.csv(paste0(work_path, "/tree_sanjay/cell_cell_", clone_i, ".cell_list.csv"), header=F)
rownames(dat) = colnames(dat) = as.vector(x$V1)

cell_assign = readRDS(paste0(work_path, "/cell_to_gastruloid/cell_assign.rds"))
colnames(cell_assign) = c("Cell", "Well", "celltype")
rownames(cell_assign) = as.vector(cell_assign$Cell)
cell_assign = cell_assign[row.names(dat),]

tree = as.phylo(hclust(as.dist(dat), "average"))
saveRDS(tree, paste0(work_path, "/tree_sanjay/", clone_i, "_tree.rds"))

### plotting the tree
dat_label = cell_assign[as.vector(tree$tip.label),]

if(clone_i == "clone05"){
    x1 = 0.055; x2 = 0.1; x3 = 0.3
} else if (clone_i == "clone25") {
    x1 = 0.055; x2 = 0.105; x3 = 0.3
} else {
    x1 = 0.06; x2 = 0.11; x3 = 0.3
}

p_rec <- ggtree(tree)

p_x <- p_rec +
    #geom_tiplab(aes(label = label), size = 1) +
    geom_fruit(data=dat_label[,c("Cell","celltype","Well")], geom=geom_tile,
               mapping=aes(y=Cell, fill=Well),
               offset = x1, pwidth = x3) +
    geom_fruit(data=dat_label[,c("Cell","celltype","Well")], geom=geom_tile,
               mapping=aes(y=Cell, fill=celltype),
               offset = x2, pwidth = x3) +
    scale_fill_manual(values = c(gastruloid_celltype_color_code, eight_wells_color_code))

try(p_x + theme(legend.position="none") +
        ggsave(paste0("~/share/tree_", clone_i, "_rec.pdf"),
               height  = 10, 
               width = 10), silent = T)


# The circular layout tree.
p_cir <- ggtree(tree, layout="fan", size=0.15, open.angle=5)

p_x <- p_cir +
    geom_fruit(data=dat_label[,c("Cell","celltype","Well")], geom=geom_tile,
               mapping=aes(y=Cell, fill=Well),
               offset = x1, pwidth = x3) +
    geom_fruit(data=dat_label[,c("Cell","celltype","Well")], geom=geom_tile,
               mapping=aes(y=Cell, fill=celltype),
               offset = x2, pwidth = x3) +
    scale_fill_manual(values = c(gastruloid_celltype_color_code, eight_wells_color_code))

try(p_x + theme(legend.position="none") +
        ggsave(paste0("~/share/tree_", clone_i, "_cir.pdf"),
               height  = 10, 
               width = 10), silent = T)



###########################################################
### Step-3: making a giant cell x tape table for each clone

clone_i = "clone32"

tree = readRDS(paste0(work_path, "/tree_sanjay/", clone_i, "_tree.rds"))
cell_include = as.vector(tree$tip.label)

p_rec <- ggtree(tree)

tip_order <- p_rec$data %>%
    dplyr::filter(isTip) %>%     # Filter for tip nodes
    dplyr::arrange(y) %>%        # Arrange by y-axis position (plot order)
    dplyr::pull(label)  

dat_all = read.csv(paste0(work_path, "/tree_sanjay/loci_", clone_i, ".csv"), header=T, row.names=1)
dat = as.matrix(dat_all[rownames(dat_all) %in% cell_include,])
N = ncol(dat)/6
tape_list = unique(paste0(unlist(lapply(colnames(dat), function(x) strsplit(x,"[.]")[[1]][1])), ".",
                          unlist(lapply(colnames(dat), function(x) strsplit(x,"[.]")[[1]][2]))))

dat_tape = NULL
for(i in 1:N){
    x = c(dat[,c(((i-1)*6 + 1):((i-1)*6+6))])
    y = rep(0, length(x))
    y[x != "None" & x != "ETY"] = 1
    dat_tape = cbind(dat_tape, y)
}
colnames(dat_tape) = tape_list

### calculate tape-tape distance and order tapes by h-clustering
dis_mat = matrix(0, N, N)
for(i in 1:N){
    print(i)
    for(j in 1:N){
        dis_mat[i,j] = 1 - sum(dat_tape[,i]*dat_tape[,j] == 1)/sqrt(sum(dat_tape[,i])*sum(dat_tape[,j]) + 1)
    }
}
rownames(dis_mat) = colnames(dis_mat) = tape_list

tree_tape = as.phylo(hclust(as.dist(dis_mat), "average"))
p_tape <- ggtree(tree_tape)
tape_order <- p_tape$data %>%
    dplyr::filter(isTip) %>%     # Filter for tip nodes
    dplyr::arrange(y) %>%        # Arrange by y-axis position (plot order)
    dplyr::pull(label)     

### ordering tapes by % that are empty
dat_tape_num = apply(dat_tape, 2, sum)
dat_tape_num = dat_tape_num[order(dat_tape_num, decreasing=T)]
tape_order = names(dat_tape_num)

#####################################
### rows are cells, columns are sites, colors are distinct editing patterns

dat_all = read.csv(paste0(work_path, "/tree_sanjay/loci_", clone_i, ".csv"), header=T, row.names=1)
cell_assign = readRDS(paste0(work_path, "/cell_to_gastruloid/cell_assign.rds"))
colnames(cell_assign) = c("Cell", "Well", "celltype")
rownames(cell_assign) = as.vector(cell_assign$Cell)

tree = readRDS(paste0(work_path, "/tree_sanjay/", clone_i, "_tree.rds"))
cell_include = as.vector(tree$tip.label)

dat = as.matrix(dat_all[rownames(dat_all) %in% cell_include,])
df = melt(as.matrix(dat))
colnames(df) = c("cell", "tape_site", "value")
df$tape = paste0(unlist(lapply(as.vector(df$tape_site), function(x) strsplit(x,"[.]")[[1]][1])), '.',
                 unlist(lapply(as.vector(df$tape_site), function(x) strsplit(x,"[.]")[[1]][2])) )

df$cell_tape = paste0(df$cell, "_", df$tape)
df$site = unlist(lapply(as.vector(df$tape_site), function(x) strsplit(x,"[.]")[[1]][3]))
unobs = df %>% filter(value == "None", site == "Site1") %>% pull(cell_tape)
df_unobs = df %>% filter(cell_tape %in% unobs)
df_obs = df %>% filter(!cell_tape %in% unobs)
df_unobs$value = "Unobserved"
df = rbind(df_obs, df_unobs)

df_x = df %>% group_by(tape, tape_site) %>% tally() %>% select(-n)
df_y = data.frame(tape = unique(df_x$tape), tape_site = paste0(unique(df_x$tape), ".Site7"))
df_x = rbind(df_x, df_y)
df_x$tape = factor(df_x$tape, levels = tape_order)
df_x = df_x[order(df_x$tape, df_x$tape_site),]
df_x$x_axis = 1:nrow(df_x)

cell_order = data.frame(cell = tip_order,
                        y_axis = 1:length(tip_order))

df = df %>% left_join(df_x[,c("tape_site", "x_axis")], by = "tape_site") %>%
    left_join(cell_order, by = "cell") 

df_x = unique(df[,c("tape_site", "x_axis")])
df_x$tape = paste0(unlist(lapply(as.vector(df_x$tape_site), function(x) strsplit(x,"[.]")[[1]][1])), ".",
                   unlist(lapply(as.vector(df_x$tape_site), function(x) strsplit(x,"[.]")[[1]][2])))
df_x = df_x %>% group_by(tape) %>% summarize(mean_x_axis = mean(x_axis))
df_y = unique(df[,c("cell", "y_axis")]) %>% arrange(y_axis)

p = df %>% 
    ggplot() + 
    geom_tile(aes(x = x_axis, y = y_axis, fill = value), color = "grey80") + 
    theme_void() +
    scale_fill_manual(values=edit_color_plate) +
    scale_x_continuous(
        breaks = as.vector(df_x$mean_x_axis),
        labels = as.vector(df_x$tape)) + 
    scale_y_continuous(
        breaks = as.vector(df_y$y_axis),
        labels = as.vector(df_y$cell))

p + theme(legend.position="none") +
    ggsave(paste0("~/share/tree_tape_", clone_i, ".png"),
           dpi = 300,
           height  = 40, 
           width = 30, limitsize = FALSE)

edit_include = as.vector(unique(df$value))
df_edit = data.frame(edit = edit_include,
                     edit_color = edit_color_plate[edit_include],
                     x = sample(1:10, size = length(edit_include), replace = T),
                     y = sample(1:10, size = length(edit_include), replace = T))
df_edit$edit = factor(df_edit$edit, levels = names(edit_color_plate)[names(edit_color_plate) %in% edit_include])
p_edit = df_edit %>% ggplot() + geom_tile(aes(x = x, y = y, fill = edit), color = "grey80") + scale_fill_manual(values=edit_color_plate)
ggsave(paste0("~/share/tree_tape_", clone_i, "_color_plate.pdf"), p_edit, width=5, height=5)

df_x = df_x[order(df_x$mean_x_axis),]
write.table(df_x[,1], paste0("~/share/tree_tape_", clone_i, "_colnames.txt"), row.names=F, col.names=F, sep="\t", quote=F)



##########################################################
### Step-4: Combining three clones to make a combined tree

clone_i = "clone05"
dat_1 = read.csv(paste0(work_path, "/tree_sanjay/cell_cell_orig_distance_", clone_i, ".csv"), header=F)
dat_2 = read.csv(paste0(work_path, "/tree_sanjay/cell_cell_comp_", clone_i, ".csv"), header=F)
dat = as.matrix(dat_1/dat_2)
dat[is.nan(dat)] = 0
x = read.csv(paste0(work_path, "/tree_sanjay/cell_cell_", clone_i, ".cell_list.csv"), header=F)
rownames(dat) = colnames(dat) = as.vector(x$V1)
matrix1 = dat

clone_i = "clone32"
dat_1 = read.csv(paste0(work_path, "/tree_sanjay/cell_cell_orig_distance_", clone_i, ".csv"), header=F)
dat_2 = read.csv(paste0(work_path, "/tree_sanjay/cell_cell_comp_", clone_i, ".csv"), header=F)
dat = as.matrix(dat_1/dat_2)
dat[is.nan(dat)] = 0
x = read.csv(paste0(work_path, "/tree_sanjay/cell_cell_", clone_i, ".cell_list.csv"), header=F)
rownames(dat) = colnames(dat) = as.vector(x$V1)
matrix2 = dat

clone_i = "clone25"
dat_1 = read.csv(paste0(work_path, "/tree_sanjay/cell_cell_orig_distance_", clone_i, ".csv"), header=F)
dat_2 = read.csv(paste0(work_path, "/tree_sanjay/cell_cell_comp_", clone_i, ".csv"), header=F)
dat = as.matrix(dat_1/dat_2)
dat[is.nan(dat)] = 0
x = read.csv(paste0(work_path, "/tree_sanjay/cell_cell_", clone_i, ".cell_list.csv"), header=F)
rownames(dat) = colnames(dat) = as.vector(x$V1)
matrix3 = dat

### combine the three matrix
combined_matrix = matrix(12, 
                         sum(nrow(matrix1), nrow(matrix2), nrow(matrix3)), 
                         sum(ncol(matrix1), ncol(matrix2), ncol(matrix3)))

combined_matrix[1:nrow(matrix1), 1:ncol(matrix1)] = matrix1
combined_matrix[(nrow(matrix1) + 1):(nrow(matrix1) + nrow(matrix2)), 
                (ncol(matrix1) + 1):(ncol(matrix1) + ncol(matrix2))] = matrix2
combined_matrix[(nrow(matrix1) + nrow(matrix2) + 1):(nrow(matrix1) + nrow(matrix2) + nrow(matrix3)), 
                (ncol(matrix1) + ncol(matrix2) + 1):(ncol(matrix1) + ncol(matrix2) + ncol(matrix3))] = matrix3

rownames(combined_matrix) = c(rownames(matrix1), rownames(matrix2), rownames(matrix3))
colnames(combined_matrix) = c(colnames(matrix1), colnames(matrix2), colnames(matrix3))

saveRDS(combined_matrix, paste0(work_path, "/tree_sanjay/combined_matrix.rds"))

### plot the combined tree
dat = readRDS(paste0(work_path, "/tree_sanjay/combined_matrix.rds"))
cell_assign = readRDS(paste0(work_path, "/cell_to_gastruloid/cell_assign.rds"))
colnames(cell_assign) = c("Cell", "Well", "celltype")
rownames(cell_assign) = as.vector(cell_assign$Cell)
cell_assign = cell_assign[row.names(dat),]

tree = as.phylo(hclust(as.dist(dat), "average"))
saveRDS(tree, paste0(work_path, "/tree_sanjay/", "combined", "_tree.rds"))

### plotting the tree
dat_label = cell_assign[as.vector(tree$tip.label),]

# The circular layout tree.
p_cir <- ggtree(tree, layout="fan", size=0.15, open.angle=5)

p_x <- p_cir +
    geom_fruit(data=dat_label[,c("Cell","celltype","Well")], geom=geom_tile,
               mapping=aes(y=Cell, fill=Well),
               offset = 0.05, pwidth = 0.5) +
    geom_fruit(data=dat_label[,c("Cell","celltype","Well")], geom=geom_tile,
               mapping=aes(y=Cell, fill=celltype),
               offset = 0.09, pwidth = 0.5) +
    scale_fill_manual(values = c(gastruloid_celltype_color_code, eight_wells_color_code))

try(p_x + theme(legend.position="none") +
        ggsave(paste0("~/share/tree_", "combined", "_cir.pdf"),
               height  = 10, 
               width = 10), silent = T)



###########################################################
### Step-5: Creating single tree for individual gastruloids


clone_list = c("clone05", "clone25", "clone32")

for(clone_i in clone_list){
    dat_1 = read.csv(paste0(work_path, "/tree_sanjay/cell_cell_orig_distance_", clone_i, ".csv"), header=F)
    dat_2 = read.csv(paste0(work_path, "/tree_sanjay/cell_cell_comp_", clone_i, ".csv"), header=F)
    dat = as.matrix(dat_1/dat_2)
    dat[is.nan(dat)] = 0
    
    x = read.csv(paste0(work_path, "/tree_sanjay/cell_cell_", clone_i, ".cell_list.csv"), header=F)
    rownames(dat) = colnames(dat) = as.vector(x$V1)
    
    cell_assign = readRDS(paste0(work_path, "/cell_to_gastruloid/cell_assign.rds"))
    colnames(cell_assign) = c("Cell", "Well", "celltype")
    rownames(cell_assign) = as.vector(cell_assign$Cell)
    cell_assign = cell_assign[row.names(dat),]
    
    well_list = unique(cell_assign$Well)
    
    for(well_i in well_list){
        print(well_i)
        
        keep = cell_assign$Well == well_i
        dat_sub = dat[keep, keep]
        print(dim(dat_sub))
        
        tree = as.phylo(hclust(as.dist(dat_sub), "average"))
        saveRDS(tree, paste0(work_path, "/tree_sanjay/", clone_i, "_", well_i, "_tree.rds"))
        
        ### plotting the tree
        dat_label = cell_assign[as.vector(tree$tip.label),]
        
        if(well_i == "Well14"){
            x1 = 0.075; x2 = 0.135; x3 = 0.3
        } else if (well_i == "Well17") {
            x1 = 0.065; x2 = 0.105; x3 = 0.3
        } else if (well_i == "Well25") {
            x1 = 0.065; x2 = 0.115; x3 = 0.3
        } else if (well_i == "Well25") {
            x1 = 0.065; x2 = 0.115; x3 = 0.3
        } else if (well_i == "Well01") {
            x1 = 0.065; x2 = 0.115; x3 = 0.3
        } else if (well_i == "Well03") {
            x1 = 0.075; x2 = 0.14; x3 = 0.3
        } else if (well_i == "Well21") {
            x1 = 0.065; x2 = 0.11; x3 = 0.3
        } else if (well_i == "Well16") {
            x1 = 0.11; x2 = 0.2; x3 = 0.3
        } else {
            x1 = 0.095; x2 = 0.172; x3 = 0.3
        }
        
        # The circular layout tree.
        p_cir <- ggtree(tree, layout="fan", size=0.15, open.angle=5)
        
        p_x <- p_cir +
            geom_fruit(data=dat_label[,c("Cell","celltype","Well")], geom=geom_tile,
                       mapping=aes(y=Cell, fill=Well),
                       offset = x1, pwidth = x3) +
            geom_fruit(data=dat_label[,c("Cell","celltype","Well")], geom=geom_tile,
                       mapping=aes(y=Cell, fill=celltype),
                       offset = x2, pwidth = x3) +
            scale_fill_manual(values = c(gastruloid_celltype_color_code, eight_wells_color_code))
        
        try(p_x + theme(legend.position="none") +
                ggsave(paste0("~/share/tree_", clone_i, "_", well_i, "_cir.pdf"),
                       height  = 10, 
                       width = 10), silent = T)
    }
}

# Well  cell_num
# Well14  1656
# Well17   489
# Well25   293
# Well01  2548
# Well03   951
# Well21   919
# Well16   414
# Well28   546




##################################################################################################
### Step-6: making a heatmap to prove cell lineages are consistent with their assigned gastruloids


clone_list = c("clone05", "clone25", "clone32")

well_node = NULL
well_node = rbind(well_node, data.frame(Well = "Well14", node = 7))
well_node = rbind(well_node, data.frame(Well = "Well17", node = 2))
well_node = rbind(well_node, data.frame(Well = "Well17", node = 6))
well_node = rbind(well_node, data.frame(Well = "Well25", node = 4))

well_node = rbind(well_node, data.frame(Well = "Well01", node = 4))
well_node = rbind(well_node, data.frame(Well = "Well01", node = 7))
well_node = rbind(well_node, data.frame(Well = "Well03", node = 6))
well_node = rbind(well_node, data.frame(Well = "Well21", node = 2))

well_node = rbind(well_node, data.frame(Well = "Well16", node = 3))
well_node = rbind(well_node, data.frame(Well = "Well28", node = 2))

###
res = NULL
res_total = NULL
for(clone_i in clone_list){
    print(clone_i)
    tree = readRDS(paste0(work_path, "/tree_sanjay/", clone_i, "_tree.rds"))
    dat = extract_ancestor_nodes(tree)
    cell_assign = readRDS(paste0(work_path, "/cell_to_gastruloid/cell_assign.rds"))
    colnames(cell_assign) = c("Cell", "Well", "celltype")
    rownames(cell_assign) = as.vector(cell_assign$Cell)
    cell_assign = cell_assign[rownames(dat),]
    res_total = rbind(res_total, cell_assign %>% group_by(Well) %>% tally())
    
    well_list = unique(cell_assign$Well)
    
    for(well_i in well_list){
        dat_i = apply(dat[,as.vector(well_node$node[well_node$Well == well_i]),drop=FALSE], 1, sum)
        tmp = table(cell_assign$Well[dat_i == 1])
        res_i = data.frame(clade = well_i,
                          well = names(tmp),
                          num = as.vector(tmp))
        res = rbind(res, res_i)
    }
}

print(res_total) ### how many cells in each gastruloid

res_frac = res %>% filter(clade == well) %>%
    left_join(res %>% group_by(clade) %>% summarize(total_num = sum(num)), by = "clade") %>%
    mutate(frac = round(100*num/total_num, 2)) %>% select(clade, frac)
# clade   frac
# Well14  99.88
# Well17  98.19
# Well25  98.98
# Well21  99.45
# Well03  98.85
# Well01  99.73
# Well28 100.00
# Well16  99.52

df = res %>% dcast(clade~well, fill=0)
rownames(df) = df[,1]
df = df[,-1]
well_list = c("Well14", "Well17", "Well25", "Well01", "Well03", "Well21", "Well16", "Well28")
df = df[well_list, well_list]

library("gplots")
library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
pdf("~/share/tree_well_heatmap.pdf", 10, 5)
heatmap.2(as.matrix(df), 
          col=Colors,
          cellnote = df,
          notecol = "black",
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



#############################################################################
### Step-7: Identifying a subset of Tapes which determine the 10 major clades

result = list()
clone_list = c("clone05", "clone25", "clone32")

for(clone_i in clone_list){
    print(clone_i)
    
    ### Part-1: reading well and clade
    tree = readRDS(paste0(work_path, "/tree_sanjay/", clone_i, "_tree.rds"))
    dat = extract_ancestor_nodes(tree)
    
    cell_assign = readRDS(paste0(work_path, "/cell_to_gastruloid/cell_assign.rds"))
    colnames(cell_assign) = c("Cell", "Well", "celltype")
    rownames(cell_assign) = as.vector(cell_assign$Cell)
    cell_assign = cell_assign[rownames(dat),]
    
    ### read 10 clades
    well_node = NULL
    if(clone_i == "clone05"){
        well_node = rbind(well_node, data.frame(well = "Well14", clade = 7))
        well_node = rbind(well_node, data.frame(well = "Well17", clade = 2))
        well_node = rbind(well_node, data.frame(well = "Well17", clade = 6))
        well_node = rbind(well_node, data.frame(well = "Well25", clade = 4))
    } else if (clone_i == "clone25"){
        well_node = rbind(well_node, data.frame(well = "Well01", clade = 4))
        well_node = rbind(well_node, data.frame(well = "Well01", clade = 7))
        well_node = rbind(well_node, data.frame(well = "Well03", clade = 6))
        well_node = rbind(well_node, data.frame(well = "Well21", clade = 2))
    } else {
        well_node = rbind(well_node, data.frame(well = "Well16", clade = 3))
        well_node = rbind(well_node, data.frame(well = "Well28", clade = 2))
    }
    node_list = as.vector(well_node$clade)
    group = rep(0, nrow(dat))
    group_well = rep(0, nrow(dat))
    for(i in node_list){
        group[dat[,i] == 1] = i
        group_well[dat[,i] == 1] = well_node$well[well_node$clade == i]
    }
    
    clade = data.frame(cell = rownames(dat), clade = group, well = group_well)
    
    ### Part-2: go through the cell x tape matrix
    dat_filter = read.csv(paste0(work_path, "/tree_sanjay/loci_", clone_i, ".csv"), header=T, row.names=1)
    print(sum(rownames(dat) != rownames(dat_filter)))
    
    tape_list = gsub(".Site1", "", colnames(dat_filter)[grep("Site1",colnames(dat_filter))])
    N = length(tape_list)
    
    res_matrix = NULL
    res_table = NULL
    
    for(i in 1:N){
        print(paste0(i,"/",N))
        dat_i = dat_filter[,c(((i-1)*6+1):(i*6))]
        
        tape_edited = data.frame(cell = rownames(dat_i), 
                                 first_edit = dat_i[,1]) %>% filter(first_edit != "ETY", first_edit != "None")
        tape_edited_num = tape_edited %>% left_join(clade, by = "cell") %>% group_by(clade) %>% tally() %>% rename(total_n = n)

        df = dat_i[,1,drop=FALSE]
        for(j in 2:6){
            df = cbind(df, paste0(df[,(j-1)], "_", dat_i[,j]))
        }
        df = data.frame(cell = rep(rownames(dat_i), 6),
                        edit = c(as.matrix(df)))
        df = df[!(grepl("ETY", df$edit) | grepl("None", df$edit)),] 
        
        df_num = df %>% left_join(clade, by = "cell") %>% group_by(clade, well, edit) %>% tally() %>% 
            left_join(tape_edited_num, by = "clade") %>% mutate(frac = n/total_n)
        
        df_num_top = df_num[df_num$frac > 0.95,] %>% group_by(clade) %>% arrange(desc(edit), .by_group = TRUE)
        if(nrow(df_num_top) == 0){
            next
        }
        
        res_table_tmp = NULL
        for(j in 1:nrow(df_num_top)){
            well_j = df_num_top$well[j]; clade_j = df_num_top$clade[j]; edit_j = df_num_top$edit[j]
            if(sum(df_num$frac[df_num$well != well_j & df_num$edit == edit_j] >= 0.05) == 0){
                edit_tmp_list = as.vector(res_table_tmp$edit[res_table_tmp$clade == clade_j])
                if(sum(grepl(edit_j, edit_tmp_list)) == 0){
                    res_table_tmp = rbind(res_table_tmp, data.frame(clade = clade_j, edit = edit_j, tape = tape_list[i]))
                    x = rep(0, nrow(dat_i))
                    x[!rownames(dat_i) %in% as.vector(tape_edited$cell)] = NA
                    x[rownames(dat_i) %in% as.vector(df$cell[df$edit == edit_j])] = 1
                    res_matrix = cbind(res_matrix, x)
                }
            }
        }
        res_table = rbind(res_table, res_table_tmp)
    }
    
    result[[clone_i]] = res_table
    
    ### Part-3: making the heatmap
    res_table$id = paste0(res_table$clade, "_", res_table$edit, "_", res_table$tape)
    res_table_uniq = res_table %>% group_by(edit, tape) %>% slice_max(order_by = clade, n = 1)
    keep = res_table$id %in% paste0(res_table_uniq$clade, "_", res_table_uniq$edit, "_", res_table_uniq$tape)
    
    res_table = res_table[keep,]
    dat_res = res_matrix[,keep]
    colnames(dat_res) = paste0(res_table$tape, ":", res_table$edit)
    rownames(dat_res) = rownames(dat_filter)
    
    d = as.matrix(dist(as.matrix(dat_res)))
    d[is.na(d)] = max(d, na.rm = TRUE) + 0.1*max(d, na.rm = TRUE)
    row_hclust = hclust(as.dist(d))
    
    d = as.matrix(dist(as.matrix(t(dat_res))))
    d[is.na(d)] = max(d, na.rm = TRUE) + 0.1*max(d, na.rm = TRUE)
    col_hclust = hclust(as.dist(d))
    
    cell_assign_sub = cell_assign[rownames(dat_res),]
    row_colors = eight_wells_color_code[as.vector(cell_assign_sub$Well)]
    
    library("gplots")
    Colors <- c("grey", "red")
    pdf(paste0("~/share/", clone_i, "_specific_tape_heatmap.pdf"), 10, 5)
    p = heatmap.2(as.matrix(dat_res), 
                  col=Colors, 
                  RowSideColors = row_colors,
                  scale="none", 
                  Rowv = as.dendrogram(row_hclust), 
                  Colv = as.dendrogram(col_hclust),
                  key=T, 
                  density.info="none", 
                  trace="none", 
                  cexRow = 1, 
                  cexCol = 0.2,
                  margins = c(5,5))
    dev.off()
    
    write.table(colnames(dat_res)[p$colInd], paste0("~/share/", clone_i, "_specific_tape_heatmap.colnames.txt"), 
                row.names=F, col.names=F, sep="\t", quote=F)
    
}  
    
saveRDS(result, paste0(work_path, "/tree_sanjay/major_clades_specific_tapes.rds"))

res = NULL
for(i in names(result)){
    x = result[[i]]; x$clone = i
    res = rbind(res, x)
}
res %>% group_by(clone, tape, edit) %>% tally() %>% ungroup() %>% group_by(clone) %>% tally()
res %>% group_by(clone, clade) %>% tally() %>% arrange(n)


###########################################################################
### Step-8: Counting the number of edited sites per cell in each gastruloid

clone_i = "clone05"
dat_1 = read.csv(paste0(work_path, "/tree_sanjay/loci_", clone_i, ".csv"), header=T, row.names=1)
df_1 = melt(as.matrix(dat_1)) %>% filter(value != "None", value != "ETY") %>%
    group_by(Var1) %>% tally()

clone_i = "clone25"
dat_2 = read.csv(paste0(work_path, "/tree_sanjay/loci_", clone_i, ".csv"), header=T, row.names=1)
df_2 = melt(as.matrix(dat_2)) %>% filter(value != "None", value != "ETY") %>%
    group_by(Var1) %>% tally()

clone_i = "clone32"
dat_3 = read.csv(paste0(work_path, "/tree_sanjay/loci_", clone_i, ".csv"), header=T, row.names=1)
df_3 = melt(as.matrix(dat_3)) %>% filter(value != "None", value != "ETY") %>%
    group_by(Var1) %>% tally()

df = rbind(df_1, df_2, df_3)
colnames(df) = c("Cell", "num")

cell_assign = readRDS(paste0(work_path, "/cell_to_gastruloid/cell_assign.rds"))
colnames(cell_assign) = c("Cell", "Well", "celltype")
rownames(cell_assign) = as.vector(cell_assign$Cell)

df = df %>% left_join(cell_assign, by = "Cell")
df$Well = factor(df$Well, levels = c("Well14", "Well17", "Well25", "Well01",
                                     "Well03", "Well21", "Well16", "Well28"))

stats <- df %>%
    group_by(Well) %>%
    summarise(
        mean = mean(num, na.rm = TRUE),
        sd = sd(num, na.rm = TRUE))

p = df %>% ggplot() + geom_histogram(aes(num, fill = Well), bins = 100) +
    facet_wrap(~Well, nrow = 2) +
    labs(x = "The # of edited sites per cell", y = "Count") +
    theme_classic(base_size = 12) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    scale_fill_manual(values = c(gastruloid_celltype_color_code, eight_wells_color_code)) +
    geom_text(data = stats, aes(x = max(df$num), y = 90, # Adjust position as needed
                                label = paste0("Mean: ", round(mean, 2), "\nSD: ", round(sd, 2))), hjust = 1, size = 4)
ggsave(paste0("~/share/Num_edited_sites_by_Well.pdf"), p, height = 6, width = 12)


df_sub = df %>% filter(Well %in% c("Well14", "Well01", "Well03", "Well21"),
                       celltype %in% c("Somites", "Spinal cord", "Mesodermal progenitors", "NMPs", "Epiblast", "Transitional cells"))
df_sub$Well = factor(df_sub$Well, levels = c("Well14", "Well01", "Well03", "Well21"))
df_sub$celltype = factor(df_sub$celltype, levels = c("Somites", "Spinal cord", "Mesodermal progenitors", "NMPs", "Epiblast", "Transitional cells"))

stats <- df_sub %>%
    group_by(Well, celltype) %>%
    summarise(
        mean = mean(num, na.rm = TRUE),
        sd = sd(num, na.rm = TRUE))

p = df_sub %>% ggplot() + geom_histogram(aes(num, fill = Well), bins = 100) +
    facet_grid(Well ~ celltype) +
    labs(x = "The # of edited sites per cell", y = "Count") +
    theme_classic(base_size = 12) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    scale_fill_manual(values = c(gastruloid_celltype_color_code, eight_wells_color_code)) +
    geom_text(data = stats, aes(x = max(df$num), y = 40, # Adjust position as needed
                                label = paste0("Mean: ", round(mean, 2), "\nSD: ", round(sd, 2))), hjust = 1, size = 4)
ggsave(paste0("~/share/Num_edited_sites_by_Well_celltype.pdf"), p, height = 8, width = 12)





