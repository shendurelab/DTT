
######################################################
### For selected branches, performing chi-squared test
### Chengxiang Qiu
### Feb-20, 2025

################################################################################
### I suggest starting by pulling the “top” set of cell divisions for each gastruloid, 
### and plot the # of cells that finally were recovdered 
### from each of these; for each, and the % (pie chart?) of each cell type

x_color_code = c("Definitive.endoderm" = "#5d71db",
                 "Endothelial.cells" = "#904bb8",
                 "Epiblast" = "#52c05a",
                 "Mesodermal.progenitors" = "#dd5732",
                 "NMPs" = "#5fc08c",
                 "Notochord" = "#a14c78",
                 "Extraembryonic.endoderm" = "#37835d",
                 "Somites" = "#dc87ba",
                 "Spinal.cord" = "#a4b46c",
                 "Transitional.cells" = "#ac4c55",
                 "Mesodermal" = "#dc87ba",
                 "Neuroectodermal" = "#a4b46c",
                 "PSC" = "#52c05a",
                 "ExE" = "#37835d",
                 "Other" = "grey",
                 "None" = "grey")


source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

work_path = "Your_work_path"

clone_list = c("clone05", "clone05", "clone05", "clone25", "clone25", "clone25", "clone32", "clone32")
well_list = c("Well14", "Well17", "Well25", "Well01", "Well03", "Well21", "Well16", "Well28")

celltype_org_list = list()

### Meso vs. Neuro
celltype_org = NULL
celltype_org = rbind(celltype_org, data.frame(celltype_org = "None", celltype = "Epiblast"))
celltype_org = rbind(celltype_org, data.frame(celltype_org = "Mesodermal", celltype = "Somites"))
celltype_org = rbind(celltype_org, data.frame(celltype_org = "Neuroectodermal", celltype = "Spinal cord"))
celltype_org = rbind(celltype_org, data.frame(celltype_org = "Neuroectodermal", celltype = "NMPs"))
celltype_org = rbind(celltype_org, data.frame(celltype_org = "Mesodermal", celltype = "Mesodermal progenitors"))
celltype_org = rbind(celltype_org, data.frame(celltype_org = "None", celltype = "Transitional cells"))
celltype_org = rbind(celltype_org, data.frame(celltype_org = "None", celltype = "Notochord"))
celltype_org = rbind(celltype_org, data.frame(celltype_org = "None", celltype = "Definitive endoderm"))
celltype_org = rbind(celltype_org, data.frame(celltype_org = "None", celltype = "Endothelial cells"))
celltype_org = rbind(celltype_org, data.frame(celltype_org = "None", celltype = "Extraembryonic endoderm"))
celltype_org_list[["meso_neuro"]] = celltype_org

### PSC vs. Other
celltype_org = NULL
celltype_org = rbind(celltype_org, data.frame(celltype_org = "PSC", celltype = "Epiblast"))
celltype_org = rbind(celltype_org, data.frame(celltype_org = "Other", celltype = "Somites"))
celltype_org = rbind(celltype_org, data.frame(celltype_org = "Other", celltype = "Spinal cord"))
celltype_org = rbind(celltype_org, data.frame(celltype_org = "Other", celltype = "NMPs"))
celltype_org = rbind(celltype_org, data.frame(celltype_org = "Other", celltype = "Mesodermal progenitors"))
celltype_org = rbind(celltype_org, data.frame(celltype_org = "PSC", celltype = "Transitional cells"))
celltype_org = rbind(celltype_org, data.frame(celltype_org = "Other", celltype = "Notochord"))
celltype_org = rbind(celltype_org, data.frame(celltype_org = "Other", celltype = "Definitive endoderm"))
celltype_org = rbind(celltype_org, data.frame(celltype_org = "Other", celltype = "Endothelial cells"))
celltype_org = rbind(celltype_org, data.frame(celltype_org = "Other", celltype = "Extraembryonic endoderm"))
celltype_org_list[["psc_other"]] = celltype_org

### ExE endoderm vs. Other
celltype_org = NULL
celltype_org = rbind(celltype_org, data.frame(celltype_org = "Other", celltype = "Epiblast"))
celltype_org = rbind(celltype_org, data.frame(celltype_org = "Other", celltype = "Somites"))
celltype_org = rbind(celltype_org, data.frame(celltype_org = "Other", celltype = "Spinal cord"))
celltype_org = rbind(celltype_org, data.frame(celltype_org = "Other", celltype = "NMPs"))
celltype_org = rbind(celltype_org, data.frame(celltype_org = "Other", celltype = "Mesodermal progenitors"))
celltype_org = rbind(celltype_org, data.frame(celltype_org = "Other", celltype = "Transitional cells"))
celltype_org = rbind(celltype_org, data.frame(celltype_org = "Other", celltype = "Notochord"))
celltype_org = rbind(celltype_org, data.frame(celltype_org = "Other", celltype = "Definitive endoderm"))
celltype_org = rbind(celltype_org, data.frame(celltype_org = "Other", celltype = "Endothelial cells"))
celltype_org = rbind(celltype_org, data.frame(celltype_org = "ExE", celltype = "Extraembryonic endoderm"))
celltype_org_list[["exe_other"]] = celltype_org


for(kk in 1:8){
    
    well_i = well_list[kk]
    print(well_i)
    
    tree = readRDS(paste0(work_path, "/tree_sanjay/", clone_list[kk], "_", well_i, "_tree.rds"))
    
    dat_anc = extract_ancestor_nodes(tree)
    dat_tip = diag(1, length(tree$tip.label), length(tree$tip.label)); rownames(dat_tip) = tree$tip.label
    dat = cbind(dat_tip, dat_anc); colnames(dat) = paste0("node_", 1:ncol(dat))
    
    cell_assign = readRDS(paste0(work_path, "/cell_to_gastruloid/cell_assign.rds"))
    colnames(cell_assign) = c("Cell", "Well", "celltype")
    rownames(cell_assign) = as.vector(cell_assign$Cell)
    cell_assign = cell_assign[row.names(dat),]
    
    for(compare_i in c("meso_neuro", "psc_other", "exe_other")){
        print(compare_i)
        celltype_org = celltype_org_list[[compare_i]]
        celltype_list = unique(as.vector(celltype_org$celltype_org))
        celltype_list = celltype_list[celltype_list != "None"]
        cell_assign_mat = matrix(0, length(celltype_list), nrow(cell_assign))
        for(i in 1:length(celltype_list)){
            cell_assign_mat[i,cell_assign$celltype %in% as.vector(celltype_org$celltype[celltype_org$celltype_org == celltype_list[i]])] = 1 
        }
        rownames(cell_assign_mat) = celltype_list
        
        x = cell_assign_mat %*% as.matrix(dat)
        x_norm = t(t(x)/apply(x, 2, sum))
        
        node_depths = data.frame(node = colnames(dat), 
                                 depth = node.depth.edgelength(tree))
        
        ### chi squared test
        cs_test = run_chi_squared_test(tree, x)
        cs_test$log10fdr = round(-log10(cs_test$fdr), 2)
        cs_test$node = as.numeric(as.vector(gsub("node_","",as.vector(cs_test$node))))
        cs_test_sig = as.vector(cs_test$node[cs_test$fdr < 0.05])
        
        saveRDS(cs_test, paste0(work_path, "/tree_sanjay/fisher_test/", compare_i, "_", well_i, "_fisher_test_res.rds"))
        
        if(length(cs_test_sig) != 0){
            ### pie chart
            edge = tree$edge
            edge_sub = edge[edge[,1] %in% cs_test_sig,]
            pie_include = unique(c(edge_sub))
            x_sub = data.frame(t(x_norm[,pie_include]))
            x_sub$node = as.numeric(as.vector(gsub("node_","",pie_include)))
            pies <- nodepie(x_sub, cols = 1:(ncol(x_sub)-1))
            pies <- lapply(pies, function(g) g+scale_fill_manual(values = x_color_code))
            
            p = ggtree(tree)
            
            tmp = p$data
            tmp = tmp %>% left_join(cs_test[,c("node", "log10fdr")], by = "node")
            p$data$log10fdr = paste0(tmp$node - length(tree$tip.label), "/", as.vector(tmp$log10fdr))
            
            p_x = p +
                geom_inset(pies, width = 0.03, height = 0.03) +
                geom_fruit(data=cell_assign[,c("Cell","celltype","Well")], geom=geom_tile,
                           mapping=aes(y=Cell, fill=celltype),
                           offset = 0.05, pwidth = 0.03, size = 10) +
                scale_fill_manual(values = c(gastruloid_celltype_color_code, eight_wells_color_code)) +
                geom_label2(aes(subset=(node %in% cs_test_sig), label=log10fdr), fill='grey', hjust = 1.3) +
                theme(legend.position="none")
            
            ggsave(paste0("~/share/", compare_i, "_", well_i, "_pie_chart.pdf"), p_x, height  = 10, width = 10)
        }
    }
}



##################################
### Plot specific events in Fig. 4

how_many_division = function(tree, node_i){
    edge = tree$edge
    root = length(tree$tip.label) + 1
    cnt = 1
    current_node = node_i
    while(current_node != root){
        current_node = edge[edge[,2] == current_node,1]
        cnt = cnt + 1
    }
    return(cnt)
}

candidate = NULL
candidate = rbind(candidate, data.frame(clone = "clone32", well = "Well16", node = 3, category = "exe_other"))
candidate = rbind(candidate, data.frame(clone = "clone32", well = "Well28", node = 6, category = "exe_other"))
candidate = rbind(candidate, data.frame(clone = "clone25", well = "Well01", node = 23, category = "psc_other"))
candidate = rbind(candidate, data.frame(clone = "clone32", well = "Well16", node = 1, category = "psc_other"))
candidate = rbind(candidate, data.frame(clone = "clone32", well = "Well28", node = 18, category = "psc_other"))
candidate = rbind(candidate, data.frame(clone = "clone5", well = "Well17", node = 1, category = "meso_neuro"))
candidate = rbind(candidate, data.frame(clone = "clone05", well = "Well14", node = 121, category = "meso_neuro"))
candidate = rbind(candidate, data.frame(clone = "clone05", well = "Well25", node = 4, category = "meso_neuro"))

result = list()

for(candidate_i in 1:nrow(candidate)){
    print(candidate_i)
    clone_i = candidate$clone[candidate_i]; well_i = candidate$well[candidate_i] 
    node_i = candidate$node[candidate_i]; category_i = candidate$category[candidate_i]

    tree = readRDS(paste0(work_path, "/tree_sanjay/", clone_i, "_", well_i, "_tree.rds"))
    
    fisher_test = readRDS(paste0(work_path, "/tree_sanjay/fisher_test/", category_i, "_", well_i, "_fisher_test_res.rds"))
    node_i = node_i + length(tree$tip.label)
    fdr_i = signif(fisher_test$fdr[fisher_test$node == node_i], 3)
 
    dat_anc = extract_ancestor_nodes(tree)
    dat_tip = diag(1, length(tree$tip.label), length(tree$tip.label)); rownames(dat_tip) = tree$tip.label
    dat = cbind(dat_tip, dat_anc); colnames(dat) = paste0("node_", 1:ncol(dat))
    
    cell_assign = readRDS(paste0(work_path, "/cell_to_gastruloid/cell_assign.rds"))
    colnames(cell_assign) = c("Cell", "Well", "celltype")
    rownames(cell_assign) = as.vector(cell_assign$Cell)
    cell_assign = cell_assign[row.names(dat),]
    
    celltype_org = celltype_org_list[[category_i]]
    celltype_list = unique(as.vector(celltype_org$celltype_org))
    celltype_list = celltype_list[celltype_list != "None"]
    cell_assign_mat = matrix(0, length(celltype_list), nrow(cell_assign))
    for(i in 1:length(celltype_list)){
        cell_assign_mat[i,cell_assign$celltype %in% as.vector(celltype_org$celltype[celltype_org$celltype_org == celltype_list[i]])] = 1 
    }
    rownames(cell_assign_mat) = celltype_list
    
    edge = tree$edge
    edge_sub = edge[edge[,1] == node_i,]
    x = t(cell_assign_mat %*% as.matrix(dat))
    
    result[[paste0(well_i, "_", node_i-length(tree$tip.label))]] = x[unique(c(edge_sub)),]
    
    p_cir = ggtree(tree, layout="fan", size=0.15, open.angle=5) +
        geom_fruit(data=cell_assign[,c("Cell","celltype","Well")], geom=geom_tile,
                   mapping=aes(y=Cell, fill=celltype), offset = 0.05, pwidth = 0.05) +
        scale_fill_manual(values = c(gastruloid_celltype_color_code, eight_wells_color_code)) +
        geom_label2(aes(subset=(node == node_i), label=node), fill='red', size = 1) +
        geom_label2(aes(subset=(node %in% c(1:nrow(dat))[dat[,node_i] == 1]), label=node), fill='red', size = 0.3) 
    ggsave(paste0("~/share/", well_i, "_subtree_", node_i-length(tree$tip.label), "_cir.pdf"), p_cir, height = 10, width = 10)
    
    #### Create the subtree by keeping only the selected tips
    selected_tips <- rownames(dat)[dat[,node_i] == 1]
    sub_tree <- drop.tip(tree, setdiff(tree$tip.label, selected_tips))
    
    cell_assign = readRDS(paste0(work_path, "/cell_to_gastruloid/cell_assign.rds"))
    colnames(cell_assign) = c("Cell", "Well", "celltype")
    rownames(cell_assign) = as.vector(cell_assign$Cell)
    cell_assign = cell_assign[sub_tree$tip.label,]
    
    p_rec = ggtree(sub_tree) +
        geom_fruit(data=cell_assign[,c("Cell","celltype","Well")], geom=geom_tile,
                   mapping=aes(y=Cell, fill=celltype), offset = 0.05, pwidth = 0.05) +
        scale_fill_manual(values = c(gastruloid_celltype_color_code, eight_wells_color_code)) +
        theme(legend.position="none")
    ggsave(paste0("~/share/", well_i, "_subtree_", node_i-length(tree$tip.label), "_rec.pdf"), p_rec, height = 10, width = 10)
    
}







