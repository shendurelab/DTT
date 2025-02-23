
#####################################################################
### Performing the self-proliferation score and differentiation score
### Chengxiang Qiu
### Feb-20, 2025

###############################################################################
### Performing the self-proliferation score and differentiation score from Choi

source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

work_path = "Your_work_path"

tree = readRDS(paste0(work_path, "/tree_sanjay/combined_tree.rds"))

cell_assign = readRDS(paste0(work_path, "/cell_to_gastruloid/cell_assign.rds"))
colnames(cell_assign) = c("Cell", "Well", "celltype")
rownames(cell_assign) = as.vector(cell_assign$Cell)
dat_label = cell_assign[as.vector(tree$tip.label),]

Cell.Annotation = dat_label

distance.matrix.MNN <- castor::get_all_pairwise_distances(tree, only_clades = 1:length(tree$tip.label), as_edge_counts = TRUE)
rownames(distance.matrix.MNN) <- tree$tip.label
colnames(distance.matrix.MNN) <- tree$tip.label

distance.matrix.MNN <- as.data.frame(distance.matrix.MNN)
distance.matrix.MNN$A <- rownames(distance.matrix.MNN)
distance.matrix.MNN <- pivot_longer(distance.matrix.MNN, cols = -A, names_to = "B", values_to = "tree_distance")

Tree.Neighbors <- select(Cell.Annotation, c('Cell','celltype'))
colnames(Tree.Neighbors) <- c('A','TypeA')
distance.matrix.MNN <- left_join(distance.matrix.MNN, Tree.Neighbors, by = 'A')

colnames(Tree.Neighbors) <- c('B','TypeB')
distance.matrix.MNN <- left_join(distance.matrix.MNN, Tree.Neighbors, by = 'B')

CX_annot_orders <- c("Epiblast","Transitional cells","NMPs","Mesodermal progenitors",'Somites',
                     'Spinal cord','Definitive endoderm',"Endothelial cells","Notochord","Extraembryonic endoderm")

distance.matrix.MNN$TypeA <- factor(distance.matrix.MNN$TypeA, levels = CX_annot_orders)
distance.matrix.MNN$TypeB <- factor(distance.matrix.MNN$TypeB, levels = CX_annot_orders)
distance.matrix.MNN <- drop_na(distance.matrix.MNN)
distance.matrix.MNN$TypeAB <- paste0(distance.matrix.MNN$TypeA,'--',distance.matrix.MNN$TypeB)

##########################################################
# Getting Tree-Mutual Nearest Neighbor here by selecting cell-pairs with tree-distance == 2
##########################################################s

Tree.Mutual.Neighbors <- filter(distance.matrix.MNN, tree_distance == 2) %>%
  filter(TypeA == TypeB)

#distance.matrix.MNN <- pivot_longer(distance.matrix.MNN, cols = -A, names_to = "B", values_to = "clonal_distance") %>%
#  filter(A != B) %>%
#  group_by(A) %>%
#  arrange(clonal_distance) %>%
#  slice_head(n = 1) %>%
#  ungroup()


Tree.Mutual.Neighbors.Summary <- as.data.frame(table(Tree.Mutual.Neighbors$TypeA))
Tree.Raw.Freq.Summary <- as.data.frame(table(Cell.Annotation$celltype))
Tree.Mutual.Neighbors.Summary <- left_join(Tree.Mutual.Neighbors.Summary, Tree.Raw.Freq.Summary, by = 'Var1') %>%
  filter(Freq.x > 2) %>%
  filter(Freq.y > 0)
colnames(Tree.Mutual.Neighbors.Summary) <- c('celltype','MNN','Raw')
Tree.Mutual.Neighbors.Summary$Freq <- Tree.Mutual.Neighbors.Summary$MNN / Tree.Mutual.Neighbors.Summary$Raw
Tree.Mutual.Neighbors.Summary$LogFreq <- log10(Tree.Mutual.Neighbors.Summary$Freq)

Tree.Mutual.Neighbors.Summary$celltype <- factor(Tree.Mutual.Neighbors.Summary$celltype, levels = CX_annot_orders)
ggplot(Tree.Mutual.Neighbors.Summary, aes(x=celltype, y=Freq, fill = celltype)) + 
  geom_bar(stat="identity") +
  theme_classic() +
  scale_fill_manual(values = gastruloid_celltype_color_code) +
  ylab("Shared cell-type frequency") + 
  xlab("Cell state") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('~/share/mGASv3_MNN_self_241010.pdf', width = 4.5, height = 3)

### the plot displays the cell-type frequency table of the pairs with distance == 2, and they are from the same cell type
### dividing by
### the cell-type frequency table of all cells (as background)



Tree.Mutual.Neighbors.Others <- filter(distance.matrix.MNN, tree_distance == 2) %>%
    #filter(TypeA != TypeB) %>%
    filter(TypeA %in% CX_annot_orders) %>%
    filter(TypeB %in% CX_annot_orders)

Tree.Mutual.Neighbors.Others.Summary <- as.data.frame(table(Tree.Mutual.Neighbors.Others$TypeAB))
Tree.TypeAB.raw <- as.data.frame(table(distance.matrix.MNN$TypeAB))
colnames(Tree.TypeAB.raw) <- c('Var1','Raw')
Tree.Mutual.Neighbors.Others.Summary <- left_join(Tree.Mutual.Neighbors.Others.Summary, Tree.TypeAB.raw, by = join_by(Var1))
Tree.Mutual.Neighbors.Others.Summary <- Tree.Mutual.Neighbors.Others.Summary %>%
    filter(Freq > 2) %>%
    filter(Raw > 0)

Tree.Mutual.Neighbors.Others.Summary$Norm <- Tree.Mutual.Neighbors.Others.Summary$Freq / Tree.Mutual.Neighbors.Others.Summary$Raw

ggplot(Tree.Mutual.Neighbors.Others.Summary, aes(x=Var1, y=Norm, fill = Var1)) + 
    geom_bar(stat="identity") +
    theme_classic() +
    #  scale_fill_manual(values = CX_annot_colors) +
    ylab("Freq. mutual nearest neighbor") + 
    xlab("Cell state") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave('~/share/mGASv3_MNN_Other_241010.pdf', width = 40, height = 8)


Tree.Mutual.Neighbors.Others.Summary <- separate(Tree.Mutual.Neighbors.Others.Summary, col = Var1, into = c('TypeA','TypeB'), sep = '--')
Tree.Mutual.Neighbors.Others.Summary$LogNorm <- log10(Tree.Mutual.Neighbors.Others.Summary$Norm)

Tree.Mutual.Neighbors.Others.Summary$TypeA <- factor(Tree.Mutual.Neighbors.Others.Summary$TypeA, levels = CX_annot_orders)
Tree.Mutual.Neighbors.Others.Summary$TypeB <- factor(Tree.Mutual.Neighbors.Others.Summary$TypeB, levels = CX_annot_orders)


minLogNorm <- min(Tree.Mutual.Neighbors.Others.Summary$LogNorm)
maxLogNorm <- max(Tree.Mutual.Neighbors.Others.Summary$LogNorm)

Tree.Mutual.Neighbors.Others.Summary %>%
  ggplot(aes(x=TypeA, y=TypeB, fill=LogNorm)) +
  geom_tile(col="white") +
  scale_fill_distiller(palette = 5, type = "div", limits=c(minLogNorm,maxLogNorm)) +
  theme_classic() +
  scale_y_discrete(limits=rev) +
  labs(fill="Log10(Norm. Freq.)") +
  xlab("Cell type") + ylab("Cell type") +
  theme(aspect.ratio = 1) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('~/share/mGASv3_MNN_other_241010.pdf', width = 6, height = 4.5)





Tree.Nearest.Neighbors.Self <- filter(distance.matrix.MNN, TypeA == TypeB) %>%
  filter(A != B) %>%
  group_by(A) %>%
  arrange(tree_distance) %>%
  slice_head(n = 1) %>%
  ungroup()

Tree.Nearest.Neighbors.Self.Summary <- Tree.Nearest.Neighbors.Self %>%
  group_by(TypeA) %>%
  summarise(
    mean_tree_dist = mean(tree_distance),
    sd_tree_dist = sd(tree_distance),
    max_tree_dist = max(tree_distance),
    min_tree_dist = min(tree_distance)
  ) %>%
  ungroup()


Tree.Nearest.Neighbors.Self.Summary$TypeA <- factor(Tree.Nearest.Neighbors.Self.Summary$TypeA, levels = c('Spinal cord',"NMPs","Mesodermal progenitors",'Somites','Definitive endoderm','Notochord'))
ggplot(Tree.Nearest.Neighbors.Self.Summary, aes(x=TypeA, y=mean_tree_dist, fill = TypeA)) + 
  geom_bar(stat="identity") +
  theme_classic() +
  scale_fill_manual(values = CX_annot_colors) +
  ylab("Mean distance (tree)") + 
  xlab("Cell state") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('Plot/mGASv3L1_SelfNearN_240330.pdf', width = 4.5, height = 3)




Tree.Nearest.Neighbors.Other <- filter(distance.matrix.MNN, A != B) %>%
  group_by(A) %>%
  arrange(tree_distance) %>%
  slice_head(n = 1) %>%
  ungroup()

Tree.Nearest.Neighbors.Other.Summary <- Tree.Nearest.Neighbors.Other %>%
  group_by(TypeAB) %>%
  summarise(
    mean_tree_dist = mean(tree_distance),
    sd_tree_dist = sd(tree_distance),
    max_tree_dist = max(tree_distance),
    min_tree_dist = min(tree_distance)
  ) %>%
  ungroup()


#Tree.Nearest.Neighbors.Self.Summary$TypeA <- factor(Tree.Nearest.Neighbors.Self.Summary$TypeA, levels = c('Spinal cord',"NMPs","Mesodermal progenitors",'Somites','Definitive endoderm','Notochord'))
ggplot(Tree.Nearest.Neighbors.Other.Summary, aes(x=TypeAB, y=mean_tree_dist, fill = TypeAB)) + 
  geom_bar(stat="identity") +
  theme_classic() +
  #  scale_fill_manual(values = CX_annot_colors) +
  ylab("Mean distance (tree)") + 
  xlab("Cell state") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('Plot/mGASv3L1_SelfNearN_240330.pdf', width = 4.5, height = 3)


ggplot(Tree.Nearest.Neighbors.Other.Summary, aes(x=TypeAB, y=mean_tree_dist, fill = TypeAB)) + 
  geom_bar(stat="identity") +
  theme_classic() +
  #  scale_fill_manual(values = CX_annot_colors) +
  ylab("Mean distance (tree)") + 
  xlab("Cell state") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





Tree.Nearest.Neighbors.Self.Summary <- as.data.frame(table(Tree.Nearest.Neighbors.Self$TypeA) / table(Cell.Annotation$celltype))
Tree.Nearest.Neighbors.Self.Summary$Var1 <- factor(Tree.Mutual.Neighbors.Summary$Var1, levels = c('Spinal cord',"NMPs","Mesodermal progenitors",'Somites','Definitive endoderm','Notochord'))








Tree.Mutual.Neighbors2 <- as.data.frame(table(Tree.Mutual.Neighbors2$TypeAB) / 
                                          table(filter(distance.matrix.MNN, TypeA != TypeB)$TypeAB))
CellTypePair <- c('Neural_tube_2-Neural_tube_1','Neural_tube_2-NMP','Neural_tube_2-pPSM','Neural_tube_2-aPSM','Neural_tube_2-Somites',
                  'Neural_tube_1-NMP','Neural_tube_1-pPSM','Neural_tube_1-aPSM','Neural_tube_1-Somites',
                  'NMP-pPSM','NMP-aPSM','NMP-Somites',
                  'pPSM-aPSM','pPSM-Somites','aPSM-Somites')
NN.DiffCellAnn.group1 <- filter(NN.DiffCellAnn.group1, Var1 %in% CellTypePair)
NN.DiffCellAnn.group1$Var1 <- factor(NN.DiffCellAnn.group1$Var1, levels = CellTypePair)

ggplot(NN.DiffCellAnn.group1, aes(x=Var1, y=Freq, fill = Var1)) + 
  geom_bar(stat="identity") +
  theme_classic() +
  #  scale_fill_manual(values = Sam_annot_colors) +
  ylab("Different cell-type frequency") + 
  xlab("Cell state") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
ggsave('Plot/mGASv3_Group1_NNdiff.pdf', width = 6, height = 4)




