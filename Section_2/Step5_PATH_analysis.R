
#################################################
### Performing PATH analysis on the combined tree
### Chengxiang Qiu
### Feb-20, 2025


source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

work_path = "Your_work_path"

library(PATH)
library(qlcMatrix)

tree = readRDS(paste0(work_path, "/tree_sanjay/combined_tree.rds"))

cell_assign = readRDS(paste0(work_path, "/cell_to_gastruloid/cell_assign.rds"))
colnames(cell_assign) = c("Cell", "Well", "celltype")
rownames(cell_assign) = as.vector(cell_assign$Cell)
dat_label = cell_assign[as.vector(tree$tip.label),]


dat_label_x = dcast(dat_label[,c("Cell", "celltype")] %>% group_by(Cell, celltype) %>% tally(), Cell~celltype, fill = 0)
rownames(dat_label_x) = dat_label_x[,1]
dat_label_x = dat_label_x[,-1]
dat_label_x = as.matrix(dat_label_x)
dat_label_x = dat_label_x[as.vector(tree$tip.label),]
### cell x celltype matrix, values are 0 or 1

W = one_node.tree.dist(tree)
# W = inv.tree.dist(tree)
# W = exp.tree.dist(tree)
# W = 1 - dat_dist; W = W/sum(W)
PATH_res = xcor(scale(dat_label_x), W)


CX_annot_orders <- c("Epiblast","Transitional cells","NMPs","Mesodermal progenitors",'Somites',
                     'Spinal cord','Definitive endoderm',"Endothelial cells","Notochord","Extraembryonic endoderm")

df_1 = data.frame(melt(PATH_res$Morans.I, value.name = "I"))
df_2 = data.frame(melt(PATH_res$Z.score, value.name = "Z"))
df = df_1 %>% left_join(df_2, by = c("Var1", "Var2"))
df$Var1 = factor(df$Var1, levels = CX_annot_orders)
df$Var2 = factor(df$Var2, levels = CX_annot_orders)


maxz <- max(abs(df$Z))
df %>%
    filter(Var2 == Var1) %>%
    ggplot(aes(x=Var1, y=Z, fill=Var1)) +
    geom_bar(stat="identity") +
    #  ylim(c(-maxz,maxz)) +
    theme_classic() +
    #scale_fill_brewer(palette = "Set2", type = "div") +
    scale_fill_manual(values = gastruloid_celltype_color_code) +
    labs(fill="Cell type") +
    ylab("Phylogenetic auto-correlation\n(z score)") + 
    xlab("Cell state") +
    #  geom_hline(yintercept = qnorm(0.05, lower.tail = F), col="black", lty=2) +
    ggtitle("Cell-type heritability", 
            subtitle = "z score")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('~/share/mGASv3_PATH_autoC.pdf', width = 6, height = 4)


df %>%
    ggplot(aes(x=Var1, y=Var2, fill=Z)) +
    geom_tile(col="white") +
    scale_fill_distiller(palette = 5, type = "div",limits=c(-maxz,maxz)) +
    #scale_fill_distiller(palette = 5, type = "div",limits=c(-35,35)) +
    theme_classic() +
    scale_y_discrete(limits=rev) +
    labs(fill="Phylogenetic\ncorrelation\nz score") +
    xlab("Cell type") + ylab("Cell type") +
    theme(aspect.ratio = 1) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('~/share/mGASv3_PATH_crossC.pdf', width = 6, height = 4.5)


pd_sub = dat_label %>% group_by(celltype) %>% tally()
pd_sub$celltype = factor(pd_sub$celltype, levels = CX_annot_orders)
p = pd_sub %>%
    ggplot(aes(x=celltype, y=n, fill=celltype)) +
    geom_bar(stat="identity") +
    labs(x="cell type", y="cell number", title="cell number") +
    theme_classic(base_size = 10) +
    scale_fill_manual(values=gastruloid_celltype_color_code) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black",angle = 45, hjust = 1), axis.text.y = element_text(color="black")) +
    ggsave("~/share/mGASv3_PATH_cell_num.pdf",
           height  = 4, 
           width = 6)


#############################
### Infer cell state dynamics

celltype = data.frame(celltype = CX_annot_orders,
                      states = 1:length(CX_annot_orders))
dat_label = dat_label %>% left_join(celltype, by = "celltype")
rownames(dat_label) = as.vector(dat_label$Cell)
dat_label = dat_label[as.vector(tree$tip.label),]
tree$states = as.vector(dat_label$states)

Pinf = PATH.inference(tree = tree, 
                       cell_states = "states", 
                       nstates = nrow(celltype))
trans_prop = Pinf$Pt
rownames(trans_prop) = colnames(trans_prop) = as.vector(celltype$celltype)

saveRDS(trans_prop, paste0(work_path, "/tree_sanjay/trans_prop_celltypes.rds"))

trans_prop %>% as.matrix() %>% melt() %>% filter(Var1 != Var2) %>%
    ggplot(aes(x=Var1, y=Var2, fill=value)) +
    geom_tile(col="white") +
    scale_fill_distiller(palette = 5, type = "div",limits=c(0,0.3)) +
    #scale_fill_distiller(palette = 5, type = "div",limits=c(-35,35)) +
    theme_classic() +
    scale_y_discrete(limits=rev) +
    labs(fill="Cell-type\ntransition\nprobolity") +
    xlab("Cell type") + ylab("Cell type") +
    theme(aspect.ratio = 1) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('~/share/mGASv3_PATH_trans_prop.pdf', width = 6, height = 4.5)


dat_trans = melt(trans_prop)
exclude_celltype = c("Definitive endoderm","Endothelial cells","Notochord","Extraembryonic endoderm")
dat_trans = subset(dat_trans, Var1 != Var2)
dat_trans = subset(dat_trans, !Var1 %in% exclude_celltype)
dat_trans = subset(dat_trans, !Var2 %in% exclude_celltype)
dat_trans = dat_trans[order(dat_trans$value, decreasing = T),]
dat_trans = dat_trans[dat_trans$value > 0.08,]
length(unique(c(dat_trans$Var1, dat_trans$Var2)))

res = NULL; res_list = NULL
for(i in 1:nrow(dat_trans)){
    if(!paste(dat_trans$Var1[i], ":", dat_trans$Var2[i]) %in% res_list){
        res = rbind(res, dat_trans[i,])
        res_list = c(res_list, paste(dat_trans$Var1[i], ":", dat_trans$Var2[i]))
        res_list = c(res_list, paste(dat_trans$Var2[i], ":", dat_trans$Var1[i]))
    }
}

write.table(res, "~/share/dat_trans.txt", row.names=F, quote=F, sep="\t")

write.table(dat_trans[dat_trans$value > 0.08,], "~/share/dat_trans.txt", row.names=F, quote=F, sep="\t")


##################################################################################
### comparing the distance of lineage and transcriptome for individual gastruloids

source("~/work/scripts/utils.R")
work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape/data_mGASv3"

library(ape)
library(castor)
library(distances)

combined_matrix = readRDS(paste0(work_path, "/tree_sanjay/combined_matrix.rds"))
### 7816 cells in total 

cell_assign = readRDS(paste0(work_path, "/cell_to_gastruloid/cell_assign.rds"))
colnames(cell_assign) = c("Cell", "Well", "celltype")
rownames(cell_assign) = as.vector(cell_assign$Cell)
cell_label = cell_assign[rownames(combined_matrix),]

obj = readRDS(paste0(work_path, "/obj_processed.rds"))
Cell = unlist(lapply(colnames(obj), function(x) strsplit(x,"[_]")[[1]][3])) 
x = rep(1, ncol(obj))
x[obj$experiment_id == "mGASv3_mGas2"] = 2
Cell = paste0(Cell, x)
pca_coor = Embeddings(obj, reduction = "pca")
rownames(pca_coor) = Cell

well_list = unique(cell_label$Well)

res = list()

for(well_i in well_list){
    print(well_i)
    
    cell_label_i = cell_label[cell_label$Well == well_i,]
    pca_coor_i = pca_coor[as.vector(cell_label_i$Cell),]
    mat_dist_i = combined_matrix[as.vector(cell_label_i$Cell), as.vector(cell_label_i$Cell)]
    
    tree_i = as.phylo(hclust(as.dist(mat_dist_i), "average"))
    dat_i = castor::get_all_pairwise_distances(tree_i, 
                                               only_clades = 1:length(tree_i$tip.label), 
                                               as_edge_counts = TRUE)
    rownames(dat_i) = colnames(dat_i) = 1:length(tree_i$tip.label)
    dat_i = melt(as.matrix(dat_i))
    dat_i = dat_i[dat_i$Var1 < dat_i$Var2,]
    dat_i$cell_1 = tree_i$tip.label[dat_i$Var1]
    dat_i$cell_2 = tree_i$tip.label[dat_i$Var2]
    dat_i$lineage_dist = dat_i$value
    dat_i$Var1 = dat_i$Var2 = dat_i$value = NULL
    
    ### Eucliean distances
    euclidean_dist = as.matrix(distances(as.matrix(pca_coor_i)))
    rownames(euclidean_dist) = colnames(euclidean_dist) = 1:nrow(pca_coor_i)
    euclidean_dist = melt(as.matrix(euclidean_dist))
    euclidean_dist = euclidean_dist[euclidean_dist$Var1 < euclidean_dist$Var2,]
    euclidean_dist$cell_1 = rownames(pca_coor_i)[euclidean_dist$Var1]
    euclidean_dist$cell_2 = rownames(pca_coor_i)[euclidean_dist$Var2]
    euclidean_dist$euclidean_dist = euclidean_dist$value
    euclidean_dist$Var1 = euclidean_dist$Var2 = euclidean_dist$value = NULL
    
    print(sum(dat_i$cell_1 != euclidean_dist$cell_1))
    print(sum(dat_i$cell_2 != euclidean_dist$cell_2))
    
    dat_i$euclidean_dist = euclidean_dist$euclidean_dist
    
    res[[well_i]] = dat_i
    
}


df = NULL

for(well_i in well_list){
    df_i = res[[well_i]]
    x = as.vector(df_i$lineage_dist)
    x[df_i$lineage_dist >= 10] = 10
    df_i$lineage_dist = as.vector(x)
    df_i = df_i %>% group_by(lineage_dist) %>% 
        summarize(mean_dist = mean(euclidean_dist), median_dist = median(euclidean_dist)) %>%
        mutate(well = well_i)
    df = rbind(df, df_i)
}

df$well = factor(df$well, levels = paste0("Well", c("14","17","25","01","03","21","16","28")))
df$lineage_dist = factor(df$lineage_dist, levels = 2:10)

p = df %>% ggplot(aes(x=lineage_dist, y=mean_dist, fill=well)) + 
    geom_bar(stat='identity') + 
    scale_fill_manual(values=eight_wells_color_code) +
    facet_grid(.~well) + labs(x='Cell-cell distance in the lineage tree',y='Cell-cell Elucidean distance in the transcriptome space') +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
    

pdf("~/share/dist_lineage_transcriptome.pdf", 10, 5)
print(p)
dev.off()


