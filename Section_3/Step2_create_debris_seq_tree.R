
##############################################################
### Calculating well-well distance and creating the phylo-tree
### Chengxiang Qiu
### Feb-20, 2025

### we first create the Well-DTT matrix from Debris-seq data using the below Python script
### python DebrisSeqTAPE_process.py

### The output can be downloaded directly from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/sam_tape
### Step1File_DebrisSeqTAPE_mGASv5_uncollapsed.csv

################################################################
### Step-1, identifying dominant edits for each DTT in each Well

source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

work_path = "Your_work_path"

dat_orig = read.csv(paste0(work_path, "/Step1File_DebrisSeqTAPE_mGASv5_uncollapsed.csv"), header=F)

### no filtering, using all the 76 DTTs
dat = dat_orig
colnames(dat) = c("Well", "TargetBC", "Site1", "Site2", "Site3", "Site4", "Site5", "Site6", "nRead")

### Replace all ETYs with None for simplicity
tmp = as.vector(dat$Site1)
dat$Site1[tmp == "ETY"] = "None"

### As each well and each targetBC has a different # of reads, let's first normalize 
### with respect to that, such that each well-targetBC combination has the same number of reads.
dat_sum = dat %>% group_by(Well, TargetBC) %>% summarize(nRead_sum = sum(nRead))

dat = dat %>% left_join(dat_sum, by = c("Well", "TargetBC")) %>%
    mutate(nRead_norm = nRead/nRead_sum)

bps <- expand.grid(bp1 = c("A", "G", "C", "T"),
                   bp2 = c("A", "G", "C", "T"),
                   bp3 = c("A", "G", "C", "T"))
bps_list <- c(apply(bps, 1, paste0, collapse = ""), "None")

well_list = unique(dat$Well)

tape_list = unique(dat$TargetBC)

result = NULL

for(i in well_list){
    print(i)
    for(j in tape_list){
        for(k in paste0("Site", 1:6)){
            dat_x = dat[(dat$Well == i & dat$TargetBC == j), c(k, "nRead_norm")]
            colnames(dat_x) = c("edit", "nRead_norm")
            dat_x = dat_x %>% group_by(edit) %>% summarise(total = sum(nRead_norm)) %>%
                arrange(desc(total)) %>% filter(edit != "None", total >= 0.2)
            if(nrow(dat_x) > 0){
                dat_x$Well = i
                dat_x$TargetBC = j
                dat_x$Site = k
                dat_x = dat_x[,c("Well","TargetBC","Site","edit","total")]
                result = rbind(result, dat_x)
            }
        }
    }
}

saveRDS(result, paste0(work_path, "/data_mGASv5/tree/well_tape_dominant_edits.rds"))


########################################
### Step-2: Calculating jaccard distance

res = list()
for(i in well_list){
    res_i = result[result$Well == i,]
    res[[i]] = paste0(res_i$TargetBC, "-", res_i$Site, "-", res_i$edit)
}

jac_sim = matrix(0, length(well_list), length(well_list))
for(i in 1:length(well_list)){
    for(j in 1:length(well_list)){
        x = res[[well_list[i]]]
        y = res[[well_list[j]]]
        jac_sim[i,j] = sum(x %in% y)/length(unique(c(x,y)))
    }
}
rownames(jac_sim) = colnames(jac_sim) = well_list

dis_mat = 1 - jac_sim

tree = as.phylo(hclust(as.dist(1 - jac_sim), "average"))
saveRDS(tree, paste0(work_path, "/data_mGASv5/tree/tree_well_tape_dominant_edits.rds"))

p_rec <- ggtree(tree) +
    geom_tiplab(aes(label = label), size = 3) +
    theme(legend.position="none")

ggsave(paste0("~/share/tree_8_DTTs.pdf"), p_rec,
       height  = 10, 
       width = 5)


#######################################################
### Step-3: Making a cell x tape table for visulization

dat_new = data.frame(Well = result$Well,
                     edit = paste0(result$TargetBC, "-", result$Site, "-", result$edit))

x = dat_new %>% group_by(edit) %>% tally() %>% as.data.frame()
x$TargetBC = unlist(lapply(as.vector(x$edit), function(x) strsplit(x,"[-]")[[1]][1])) 
x$site = unlist(lapply(as.vector(x$edit), function(x) strsplit(x,"[-]")[[1]][2])) 
x = x %>% arrange(desc(n))
x$x_axis = 1:nrow(x)

well_list = p_rec$data %>%
    dplyr::filter(isTip) %>%     # Filter for tip nodes
    dplyr::arrange(y) %>%        # Arrange by y-axis position (plot order)
    dplyr::pull(label)
y = data.frame(Well = well_list, y_axis = 1:length(well_list))

dat_new = dat_new %>% left_join(y, by = "Well") %>%
    left_join(x, by = "edit")
### 108 wells x 7,146 TapeBC-Site-Edit combinations


p = dat_new %>% 
    ggplot() + 
    geom_tile(aes(x = x_axis, y = y_axis), color = "red") + 
    theme_void() +
    scale_x_continuous(
        breaks = as.vector(x$x_axis),
        labels = as.vector(x$edit)) + 
    scale_y_continuous(
        breaks = as.vector(y$y_axis),
        labels = as.vector(y$Well)) +
    theme(legend.position="none")

ggsave(paste0("~/share/mGASv5_debris_well_tape.png"), p, 
       dpi = 300,
       height  = 10, 
       width = 50, limitsize = FALSE)

write.table(x[,1], paste0("~/share/mGASv5_debris_well_tape.colnames.txt"), row.names=F, col.names=F, sep="\t", quote=F)

y = y[order(y$y_axis, decreasing=T),]
write.table(y[,1], paste0("~/share/mGASv5_debris_well_tape.rownames.txt"), row.names=F, col.names=F, sep="\t", quote=F)


#########################################################################
### Step-4: Applying bootstrap to create 100 trees, validating robustness

tree_orig = readRDS(paste0(work_path, "/data_mGASv5/tree/tree_well_tape_dominant_edits.rds"))
dat_orig = extract_ancestor_nodes(tree_orig)

result = readRDS(paste0(work_path, "/data_mGASv5/tree/well_tape_dominant_edits.rds"))
well_list = unique(result$Well)

bootstrap_time = 100

TBE = NULL

for(iter in 1:bootstrap_time){
    
    print(iter)
    
    tape_list = unique(result$TargetBC)
    
    index = sample(1:length(tape_list), length(tape_list), replace = T)
    
    result_bs = NULL
    for(i in 1:length(index)){
        result_i = result[result$TargetBC == tape_list[index[i]],]
        result_i$TargetBC = paste0(result_i$TargetBC, "_", i)
        result_bs = rbind(result_bs, result_i)
    }
    
    res = list()
    for(i in well_list){
        res_i = result_bs[result_bs$Well == i,]
        res[[i]] = paste0(res_i$TargetBC, "-", res_i$Site, "-", res_i$edit)
    }
    
    jac_sim = matrix(0, length(well_list), length(well_list))
    for(i in 1:length(well_list)){
        for(j in 1:length(well_list)){
            x = res[[well_list[i]]]
            y = res[[well_list[j]]]
            jac_sim[i,j] = sum(x %in% y)/length(unique(c(x,y)))
        }
    }
    rownames(jac_sim) = colnames(jac_sim) = well_list
    
    dis_mat = 1 - jac_sim
    
    tree = as.phylo(hclust(as.dist(1 - jac_sim), "average"))
    
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
    
    TBE = cbind(TBE, res$score)
}

res = as.matrix(TBE)

N = tree_orig$Nnode
N_cell = nrow(dat_orig)
N_p = NULL
for(i in 1:ncol(dat_orig)){
    N_p = c(N_p, min(sum(dat_orig[,i] == 1), sum(dat_orig[,i] == 0)))
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
### 99.1%

saveRDS(x2, paste0(work_path, "/data_mGASv5/tree/tree_well_tape_dominant_edits.TBE.rds"))


### Plot TBE on the tree
x2 = readRDS(paste0(work_path, "/data_mGASv5/tree/tree_well_tape_dominant_edits.TBE.rds"))
tree_orig = readRDS(paste0(work_path, "/data_mGASv5/tree/tree_well_tape_dominant_edits.rds"))
p_rec <- ggtree(tree_orig)

tmp = p_rec$data
x2$node = x2$node + 108
tmp = tmp %>% left_join(x2, by = "node")
p_rec$data$TBE = round(as.vector(tmp$TBE),2)

p_rec <- p_rec +
    geom_label2(aes(subset = node %in% c(108:215), label = TBE), size = 2) +
    theme(legend.position = "none")

ggsave(paste0("~/share/tree_8_DTTs_label_TBE.pdf"), p_rec,
       height  = 10, 
       width = 3)

