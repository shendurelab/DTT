
########################################################################
### Bootstrapping analysis to evaluate the robust of tree reconstruction
### Chengxiang Qiu
### Feb-20, 2025

### Creating cell-cell distance by bootstrapping
### python Perform_cell_cell_distance_bootstrapping.py

source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

work_path = "Your_work_path"

cell_assign = readRDS(paste0(work_path, "/cell_to_gastruloid/cell_assign.rds"))
colnames(cell_assign) = c("Cell", "Well", "celltype")
rownames(cell_assign) = as.vector(cell_assign$Cell)

########################################################
### Step-1: Creating clone-level tree for each bootstrap

args = commandArgs(trailingOnly=TRUE)
kk = as.numeric(args[1]) ### 1:100

clone_list = c("clone05", "clone32", "clone25")

for(clone_i in clone_list){
    print(clone_i)
    
    tree_orig = readRDS(paste0(work_path, "/tree_sanjay/", clone_i, "_tree.rds"))
    cell_include = as.vector(tree_orig$tip.label)
    
    dat_1 = read.csv(paste0(work_path, "/tree_sanjay/bootstrapping/cell_cell_distance/cell_cell_orig_distance_", clone_i, "_", kk, ".csv"), header=F)
    dat_2 = read.csv(paste0(work_path, "/tree_sanjay/bootstrapping/cell_cell_distance/cell_cell_comp_", clone_i, "_", kk, ".csv"), header=F)
    matrix_clone = as.matrix(dat_1/dat_2)
    matrix_clone[is.nan(matrix_clone)] = 0
    matrix_clone_cellname = read.csv(paste0(work_path, "/tree_sanjay/cell_cell_", clone_i, ".cell_list.csv"), header=F)
    rownames(matrix_clone) = colnames(matrix_clone) = as.vector(matrix_clone_cellname$V1)
    
    dat = matrix_clone[cell_include, cell_include]
    
    tree = as.phylo(hclust(as.dist(dat), "average"))
    saveRDS(tree, paste0(work_path, "/tree_sanjay/bootstrapping/bootstrap_tree/", clone_i, "_tree_", kk, ".rds"))
}


###################################################################
### Step-2: Counting how many clades are repeated by 100 bootstraps

### Here, instead of searching perfect matching clades (Felsenstein's bootstrap proportions; FBPs), 
### I am going to apply a different method which is called "transfer bootstrap" and report 
### transfer bootstrap expectation (TBE), Ref: https://pubmed.ncbi.nlm.nih.gov/29670290/

clone_list = c("clone05", "clone32", "clone25")

for(clone_i in clone_list){
    print(clone_i)
    
    tree_orig = readRDS(paste0(work_path, "/tree_sanjay/", clone_i, "_tree.rds"))
    dat_orig = extract_ancestor_nodes(tree_orig)
    
    tree = readRDS(paste0(work_path, "/tree_sanjay/bootstrapping/bootstrap_tree/", clone_i, "_tree_", kk, ".rds"))
    dat_anc = extract_ancestor_nodes(tree)
    dat_tip = diag(1, length(tree$tip.label), length(tree$tip.label)); rownames(dat_tip) = tree$tip.label
    dat_com = cbind(dat_tip, dat_anc); colnames(dat_com) = paste0("node_", 1:ncol(dat_com))
    dat = dat_com
    
    print(sum(tree$tip.label == tree_orig$tip.label))
    
    ham_dist = function(x, y){
        a = length(union(x, y))
        b = length(intersect(x, y))
        return(a - b)
    }
    
    res = NULL
    for(i in 1:ncol(dat_orig)){
        print(i)
        score = NULL
        a1 = c(1:nrow(dat_orig))[dat_orig[,i] == 1]
        a0 = c(1:nrow(dat_orig))[dat_orig[,i] == 0]
        for(j in 1:ncol(dat)){
            b1 = c(1:nrow(dat))[dat[,j] == 1]
            b0 = c(1:nrow(dat))[dat[,j] == 0]
            score_tmp = min(ham_dist(a1, b1), ham_dist(a1, b0))
            score = c(score, score_tmp)
            if(score_tmp == 0){
                break
            }
        }
        res = rbind(res, data.frame(node = i,
                                    score = min(score)))
    }
    
    saveRDS(res, paste0(work_path, "/tree_sanjay/bootstrapping/summary_result/", clone_i, "_TBEs_", kk, ".rds"))
    
}



####################################################
### faster way to calculate TBEs (dat_orig data dat)

### Method-1

res_x = NULL
for(i in 1:ncol(dat_orig)){
    print(i)
    score = nrow(dat_orig)
    a = dat_orig[,i]
    for(j in 1:ncol(dat)){
        b = dat[,j]
        score_tmp = min(sum(a == b), sum(a != b))
        if(score_tmp == 0){
            score = 0
            break
        }
        if(score_tmp < score){
            score = score_tmp
        }
    }
    res_x = rbind(res_x, data.frame(node = i,
                                    score = score))
}

### Method-2

dat_orig_ = 1 - dat_orig
dat_ = 1 - dat

dat_1 = t(dat_orig) %*% dat + t(dat_orig_) %*% dat_
dat_1_min = apply(dat_1, 1, min)

dat_2 = t(dat_orig_) %*% dat + t(dat_orig) %*% dat_
dat_2_min = apply(dat_2, 1, min)

dat_3 = as.matrix(cbind(dat_1_min, dat_2_min))
res_y = data.frame(node = c(1:ncol(dat_orig)), score = apply(dat_3, 1, min))


########################################################################
### Step-3: plot the distribution of bootstrapping values for each clone 
### as a function of distance from root, as well as the bootstrapping values plotted onto the tree

clone_i = "clone25"

tree_orig = readRDS(paste0(work_path, "/tree_sanjay/", clone_i, "_tree.rds"))
N = tree_orig$Nnode
dat_orig = extract_ancestor_nodes(tree_orig)
N_cell = nrow(dat_orig)

N_p = NULL
for(i in 1:ncol(dat_orig)){
    N_p = c(N_p, min(sum(dat_orig[,i] == 1), sum(dat_orig[,i] == 0)))
}

res = NULL
for(kk in 1:100){
    res_kk = readRDS(paste0(work_path, "/tree_sanjay/bootstrapping/summary_result/", clone_i, "_TBEs_", kk, ".rds"))
    res = cbind(res, as.vector(res_kk$score))
}

score = apply(as.matrix(res), 1, mean)

x2 = NULL
for(i in 1:nrow(res)){
    if(score[i] == 0){
        x2 = rbind(x2, data.frame(node = i, TBE = 1 - score[i]))
    } else {
        x2 = rbind(x2, data.frame(node = i, TBE = 1 - score[i]/(N_p[i] - 1)))
    }
}

print(sum(x2$TBE >= 0.7)/nrow(x2))
### clone05: 47.80%
### clone25: 76.75%
### clone32: 59.85%
 
node_depths = data.frame(node = 1:nrow(x2), 
                         depth = node.depth.edgelength(tree_orig)[(length(tree_orig$tip.label)+1):(length(tree_orig$tip.label)+tree_orig$Nnode)])
x2 = x2 %>% left_join(node_depths, by = "node")

p = x2 %>% ggplot() + geom_point(aes(x = TBE, y = depth)) +
    geom_vline(xintercept = 0.7, color = "red", size = 1.5) + 
    labs(x = "Transfer bootstrap expectation (TBE)", y = "Distance from root") +
    theme_classic(base_size = 12) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
ggsave(paste0("~/share/", clone_i, "_scatter_TBE_depth.pdf"), p, width = 5, height = 5)

p = x2 %>% ggplot() + geom_histogram(aes(TBE), bins = 50) +
    geom_vline(xintercept = 0.7, linetype = "longdash", color = "red", size = 1.5) + 
    labs(x = "Transfer bootstrap expectation (TBE)", y = "Count") +
    theme_classic(base_size = 12) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
ggsave(paste0("~/share/", clone_i, "_hist_TBE.pdf"), p, width = 5, height = 5)




#############################################################
### Step-4: Creating gastruloid-level tree for each bootstrap

cell_assign = readRDS(paste0(work_path, "/cell_to_gastruloid/cell_assign.rds"))
colnames(cell_assign) = c("Cell", "Well", "celltype")
rownames(cell_assign) = as.vector(cell_assign$Cell)

args = commandArgs(trailingOnly=TRUE)
kk = as.numeric(args[1]) ### 1:100

clone_list = c("clone05", "clone25", "clone32")

ham_dist = function(x, y){
    a = length(union(x, y))
    b = length(intersect(x, y))
    return(a - b)
}

for(clone_i in clone_list){
    print(clone_i)
    tree_tmp = readRDS(paste0(work_path, "/tree_sanjay/", clone_i, "_tree.rds"))
    well_list = unique(cell_assign$Well[cell_assign$Cell %in% tree_tmp$tip.label])
    
    dat_1 = read.csv(paste0(work_path, "/tree_sanjay/bootstrapping/cell_cell_distance/cell_cell_orig_distance_", clone_i, "_", kk, ".csv"), header=F)
    dat_2 = read.csv(paste0(work_path, "/tree_sanjay/bootstrapping/cell_cell_distance/cell_cell_comp_", clone_i, "_", kk, ".csv"), header=F)
    matrix_clone = as.matrix(dat_1/dat_2)
    matrix_clone[is.nan(matrix_clone)] = 0
    matrix_clone_cellname = read.csv(paste0(work_path, "/tree_sanjay/cell_cell_", clone_i, ".cell_list.csv"), header=F)
    rownames(matrix_clone) = colnames(matrix_clone) = as.vector(matrix_clone_cellname$V1)
    
    for(well_i in well_list){
        print(well_i)
        tree_orig = readRDS(paste0(work_path, "/tree_sanjay/", clone_i, "_", well_i, "_tree.rds"))
        cell_include = as.vector(tree_orig$tip.label)
        dat = matrix_clone[cell_include, cell_include]
        tree = as.phylo(hclust(as.dist(dat), "average"))
        saveRDS(tree, paste0(work_path, "/tree_sanjay/bootstrapping/bootstrap_tree/", clone_i, "_", well_i, "_tree_", kk, ".rds"))
    
        ### Counting how many clades are repeated by 100 bootstraps
        dat_orig = extract_ancestor_nodes(tree_orig)
        dat_anc = extract_ancestor_nodes(tree)
        dat_tip = diag(1, length(tree$tip.label), length(tree$tip.label)); rownames(dat_tip) = tree$tip.label
        dat_com = cbind(dat_tip, dat_anc); colnames(dat_com) = paste0("node_", 1:ncol(dat_com))
        dat = dat_com
        
        if(sum(tree$tip.label != tree_orig$tip.label) != 0){
            break
        }
        
        res = NULL
        for(i in 1:ncol(dat_orig)){
            print(i)
            score = NULL
            a1 = c(1:nrow(dat_orig))[dat_orig[,i] == 1]
            a0 = c(1:nrow(dat_orig))[dat_orig[,i] == 0]
            for(j in 1:ncol(dat)){
                b1 = c(1:nrow(dat))[dat[,j] == 1]
                b0 = c(1:nrow(dat))[dat[,j] == 0]
                score_tmp = min(ham_dist(a1, b1), ham_dist(a1, b0))
                score = c(score, score_tmp)
                if(score_tmp == 0){
                    break
                }
            }
            res = rbind(res, data.frame(node = i,
                                        score = min(score)))
        }
        saveRDS(res, paste0(work_path, "/tree_sanjay/bootstrapping/summary_result/", clone_i, "_", well_i, "_TBEs_", kk, ".rds"))
    }
}


### summary the result, for each node, save the TBE

for(clone_i in clone_list){
    print(clone_i)
    tree_tmp = readRDS(paste0(work_path, "/tree_sanjay/", clone_i, "_tree.rds"))
    cell_assign = readRDS(paste0(work_path, "/cell_to_gastruloid/cell_assign.rds"))
    colnames(cell_assign) = c("Cell", "Well", "celltype")
    rownames(cell_assign) = as.vector(cell_assign$Cell)
    well_list = unique(cell_assign$Well[cell_assign$Cell %in% tree_tmp$tip.label])

    for(well_i in well_list){
        print(well_i)
        tree_orig = readRDS(paste0(work_path, "/tree_sanjay/", clone_i, "_", well_i, "_tree.rds"))
    
        N = tree_orig$Nnode
        dat_orig = extract_ancestor_nodes(tree_orig)
        N_cell = nrow(dat_orig)
        
        N_p = NULL
        for(i in 1:ncol(dat_orig)){
            N_p = c(N_p, min(sum(dat_orig[,i] == 1), sum(dat_orig[,i] == 0)))
        }
        
        res = NULL
        for(kk in 1:100){
            res_kk = readRDS(paste0(work_path, "/tree_sanjay/bootstrapping/summary_result/", clone_i, "_", well_i, "_TBEs_", kk, ".rds"))
            res = cbind(res, as.vector(res_kk$score))
        }
        
        score = apply(as.matrix(res), 1, mean)
        
        x2 = NULL
        for(i in 1:nrow(res)){
            if(score[i] == 0){
                x2 = rbind(x2, data.frame(node = i, TBE = 1 - score[i]))
            } else {
                x2 = rbind(x2, data.frame(node = i, TBE = 1 - score[i]/(N_p[i] - 1)))
            }
        }
        
        print(paste0(well_i, ", how many clades >70% TBE: ", round(100*sum(x2$TBE >= 0.7)/nrow(x2), 2), "%"))
        saveRDS(x2, paste0(work_path, "/tree_sanjay/bootstrapping/", clone_i, "_", well_i, "_TBEs.rds"))
    }
}



