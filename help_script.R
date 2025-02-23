
print("Loading Monocle3 and Seurat")
suppressMessages(library(monocle3))
suppressMessages(library(Seurat))

print("Loading packages for regular data analysis, e.g. dplyr")
suppressMessages(library(Matrix))
suppressMessages(library(dplyr))
suppressMessages(library(reshape2))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))

print("Loading packages for plotting, e.g. ggplot2")
suppressMessages(library(ggplot2))
suppressMessages(library(plotly))
suppressMessages(library(htmlwidgets))
suppressMessages(library(gridExtra))
suppressMessages(library(viridis))
suppressMessages(library(gplots))

print("Loading packages for phylo tree analysis")
suppressMessages(library(ggtreeExtra))
suppressMessages(library(ggtree))
suppressMessages(library(PATH))
suppressMessages(library(ape))
suppressMessages(library(qlcMatrix))
suppressMessages(library(castor))
suppressMessages(library(distances))
suppressMessages(library(DescTools))

### Support data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/sam_tape
### Please let me know if you cannot find them.
### CX Qiu (cxqiu@uw.edu)

##############################################
### cell lineage tracing in gastruloid project

gastruloid_celltype_color_code = c("Cardiac mesoderm" = "#dc9436",
                                   "Definitive endoderm" = "#5d71db",
                                   "Early neurons" = "#85b937",
                                   "Endothelial cells" = "#904bb8",
                                   "PSC-like cells" = "#52c05a",
                                   "ESCs (2-cell state)" = "#cf439f",
                                   "Fibroblasts" = "#cc79dc",
                                   "Floor plate" = "#c2b73d",
                                   "Node-like cells" = "#7e72b7",
                                   "Hindbrain" = "#8e8927",
                                   "Anterior mesendoderm" = "#5a9bd5",
                                   "Mesodermal progenitors" = "#dd5732",
                                   "MHB" = "#3fc1bf",
                                   "Motor neurons" = "#dd4168",
                                   "NMPs" = "#5fc08c",
                                   "Notochord" = "#a14c78",
                                   "Extraembryonic endoderm" = "#37835d",
                                   "Somites" = "#dc87ba",
                                   "Spinal cord" = "#a4b46c",
                                   "Transitional cells" = "#ac4c55")

eight_wells_color_code = c("Well01" = "#0000FF",
                           "Well03" = "#FFA500",
                           "Well14" = "#FF0000",
                           "Well16" = "#800000", 
                           "Well17" = "#006400",
                           "Well21" = "#800080",
                           "Well25" = "#696969",
                           "Well28" = "#0d0887")


###############################
### Editing pattern color plate

bps <- expand.grid(bp1 = c("A", "G", "C", "T"),
                   bp2 = c("A", "G", "C", "T"),
                   bp3 = c("A", "G", "C", "T"))
bps_list <- apply(bps, 1, paste0, collapse = "")
edit_color_plate = c("#7fc64f", "#9645bf", "#46ca63", "#bd39aa", "#39a12d", "#4959cf", "#a8bc34", "#a074ea",
                     "#719b2a", "#dd73de", "#52c281", "#d7419a", "#559a4a", "#6951b5", "#d0b940", "#567fef",
                     "#dd9d2d", "#5168bb", "#da821e", "#4b8fcd", "#e66526", "#54b1e0", "#bc371c", "#3dc2cc",
                     "#eb5348", "#58c8a9", "#e33f78", "#3b772c", "#7b51a7", "#a9962f", "#bb86d9", "#64771e",
                     "#a54a97", "#8ab16b", "#ab2e6f", "#428f5c", "#ea72ae", "#2a673e", "#c83748", "#3a997f",
                     "#b94668", "#2b775d", "#f07a61", "#5565a1", "#d07036", "#9498df", "#74691a", "#814d93",
                     "#c1bb75", "#864e78", "#d69b54", "#b36d9e", "#5e6930", "#e295c4", "#a16c27", "#924465",
                     "#988e4e", "#e07e87", "#8e5e32", "#ef9676", "#9e4b50", "#d6a075", "#a14826", "#c16858",
                     "grey80", "grey80", "white")
names(edit_color_plate) = c(bps_list, "None", "ETY", "Unobserved")


###########################################
### published datasets used for integration

jax_celltype_color_code = c("Epithelium" = "#9c476b",
                            "Mesoderm" = "#5dbb4d",
                            "Endothelium" = "#7459c6",
                            "Cardiomyocytes" = "#a8b439",
                            "Megakaryocytes" = "#c845a1",
                            "CNS_neurons" = "#46873f",
                            "Hepatocytes" = "#bc76d5",
                            "Olfactory_sensory_neurons" = "#d59a35",
                            "Ependymal_cells" = "#6c80c9",
                            "Lung_and_airway" = "#d05130",
                            "Neuroectoderm_and_glia" = "#4eacd7",
                            "Intestine" = "#73732a",
                            "Neural_crest_PNS_glia" = "#5dc598",
                            "Primitive_erythroid" = "#d681b2",
                            "Neural_crest_PNS_neurons" = "#9cb26a",
                            "Muscle_cells" = "#d3756f",
                            "White_blood_cells" = "#379884",
                            "Eye_and_other" = "#9f5f2b",
                            "Definitive_erythroid" = "#d69f6a")

pijuan_celltype_color_code = c("Epiblast" = "#5dbb4d",
                               "ExE ectoderm" = "#5dc54f",
                               "ExVE" = "#d59a35",
                               "ExE endoderm" = "#93b730",
                               "Nascent mesoderm" = "#d05130",
                               "PGC" = "#9c476b",
                               "Rostral neurectoderm" = "#7459c6",
                               "Def. endoderm" = "#bc76d5",
                               "Mesenchyme" = "#4eacd7",
                               "Blood progenitors 2" = "#85bf71",
                               "Gut" = "#4866d3",
                               "Erythroid3" = "#cabd39",
                               "Notochord" = "#46873f",
                               "Caudal Mesoderm" =              "#dc9535",
                               "Paraxial mesoderm" =            "#446bae",
                               "Somitic mesoderm" =             "#d4552a",
                               "Allantois" =                    "#4cacdb",
                               "Forebrain/Midbrain/Hindbrain" = "#d4424c",
                               "Cardiomyocytes" =               "#73732a",
                               "Neural crest" =                 "#b4397f",
                               "Primitive Streak" =             "#347e47",
                               "EmVE" =                         "#da84cf",
                               "Parietal endoderm" = "#6c80c9",
                               "Visceral endoderm" = "#b282d5",
                               "Anterior Primitive Streak" = "#aa9732",
                               "Mixed mesoderm" = "#894e98",
                               "Surface ectoderm" = "#c845a1",
                               "Haematoendothelial progenitors" = "#88a0e5",
                               "Blood progenitors 1" = "#d681b2",
                               "ExE mesoderm" = "#67b88c",
                               "Caudal epiblast" = "#459939",
                               "Erythroid2" = "#2b7f63",
                               "Intermediate mesoderm" = "#dd7f9e",
                               "Pharyngeal mesoderm" = "#547739",
                               "Caudal neurectoderm" = "#ea8e72",
                               "Erythroid1" = "#9cb26a",
                               "Endothelium" = "#5dc598",
                               "Spinal cord" = "#bead6c",
                               "NMP" = "#b98554")

Liberali_celltype_color_code = c("Naive pluripotency" = "#9c476b",
                                 "Exiting naive pluripotency" = "#9cb26a",
                                 "Zscan4+ Artefact" = "#7459c6",
                                 "Epiblast" = "#a8b439",
                                 "Pre-somitic mesoderm" = "#c845a1",
                                 "Somite differentiation front" = "#46873f",
                                 "Ectopic pluripotency" = "#9f5f2b",
                                 "Caudal mesoderm" = "#d59a35",
                                 "Gut" = "#6c80c9",
                                 "Paraxial mesoderm" = "#d05130",
                                 "Hemogenic endothelium" = "#5dc598",
                                 "Somite" = "#73732a",
                                 "Neuromesodermal progenitors" = "#4eacd7",
                                 "Cd63+ ectoderm-like artefact" = "#d681b2",
                                 "Epiblast/primitive streak" = "#5dbb4d",
                                 "Primitive streak" = "#d3756f",
                                 "Caudal epiblast/primitive streak" = "#379884",
                                 "Anterior primitive streak/Def. endoderm" = "#bc76d5",
                                 "Caudal epiblast" = "#d69f6a")

TLS_celltype_color_code = c("PCGLC" = "#9c476b",
                            "NeuralTube2" = "#9cb26a",
                            "NeuralTube1" = "#7459c6",
                            "NMPs" = "#a8b439",
                            "pPSM" = "#c845a1",
                            "aPSM" = "#46873f",
                            "Somite-1" = "#9f5f2b",
                            "Somite0" = "#d59a35",
                            "Somite" = "#6c80c9",
                            "SomiteDermo" = "#d05130",
                            "SomiteSclero" = "#5dc598",
                            "Endothelial" = "#73732a",
                            "Endoderm" = "#4eacd7",
                            "Unknown" = "#d681b2")

Rosen_celltype_color_code = c("Primitive Streak"          = "#dc3c6e",
                              "Early Nascent Mesoderm"    = "#5cc151",
                              "Anterior Primitive Streak" = "#b15ecf",
                              "Early Posterior PSM"       = "#9bb837",
                              "Mature Endoderm"           = "#6a61cd",
                              "Caudal Epiblast"           = "#c9a93a",
                              "PGCs"                      = "#6576c1",
                              "Early Neurectoderm"        = "#df852f",
                              "Epiblast"                  = "#5fa1d8",
                              "Caudal Neurectoderm"       = "#d54637",
                              "Late Posterior PSM"        = "#42c3c2",
                              "Early Spinal Cord"         = "#d149a2",
                              "Late Neurectoderm"         = "#4f923a",
                              "Notochord"                 = "#d28fd0",
                              "Late Nascent Mesoderm"     = "#66c388",
                              "Cardiopharyngeal Mesoderm" = "#974f8a",
                              "NMPs"                      = "#747c29",
                              "Mature Somites"            = "#df7a89",
                              "Somites"                   = "#3f926d",
                              "Early Anterior PSM"        = "#a44657",
                              "Head Mesoderm"             = "#476b31",
                              "Late Spinal Cord"          = "#a8532d",
                              "Late Anterior PSM"         = "#aeb06a",
                              "Endothelium"               = "#e1956d",
                              "Neurons"                   = "#957030")

Van_celltype_color_code = c("Cardiac"               = "#be73a7",
                            "Paraxial MD"           = "#5eb956",
                            "Differentiated somite" = "#c152b8",
                            "Somite"                = "#b0b33e",
                            "Differentiation front" = "#7b66cd",
                            "PSM"                   = "#d88c38",
                            "NMPs"                  = "#678ccc",
                            "Spinal cord"           = "#ce4c33",
                            "Mesenchyme"            = "#50b9a3",
                            "Endothelium"           = "#cc4270",
                            "Allantois"             = "#558140",
                            "PGC-like/ExE EcD"      = "#c56f61",
                            "Endoderm"              = "#9c8240")


#####################################################
### Function: doing regular analysis using Seurat ###
#####################################################

doClusterSeurat <- function(obj, 
                            nfeatures = 2500, 
                            resolution = 1, 
                            k.filter = 200, 
                            doClustering = TRUE, 
                            n_dim = 30, 
                            min.dist = 0.3, 
                            n.components = 2){
    
    if(length(table(obj$group))!=1){
        
        obj.list <- SplitObject(object = obj, split.by = "group")
        
        for (i in 1:length(x = obj.list)) {
            obj.list[[i]] <- NormalizeData(object = obj.list[[i]], verbose = FALSE)
            obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]], 
                                                  selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
        }
        
        reference.list <- obj.list[names(table(obj$group))]
        obj.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:n_dim, k.filter = k.filter)
        
        obj.integrated <- IntegrateData(anchorset = obj.anchors, dims = 1:n_dim)
        
        DefaultAssay(object = obj.integrated) <- "integrated"
        
        obj <- obj.integrated 
        
    } else {
        
        obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
        obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures)
        
    }
    
    obj <- ScaleData(object = obj, verbose = FALSE)
    obj <- RunPCA(object = obj, npcs = n_dim, verbose = FALSE)
    
    if (doClustering == TRUE){
        obj <- FindNeighbors(object = obj, dims = 1:n_dim, reduction = "pca")
        obj <- FindClusters(object = obj, resolution = resolution)
    }
    
    obj <- RunUMAP(object = obj, reduction = "pca", dims = 1:n_dim, min.dist = min.dist, n.components = n.components)
    
    return(obj)
}

#####################################################################
### Function: create ancestor node x tips matrix from an ape tree ###
#####################################################################

extract_ancestor_nodes = function(phylo_tree){
    # Get the number of nodes in the tree (ancestor nodes and terminal nodes)
    # in this example, n = 1034 terminal nodes, and n = 1033 ancestor nodes
    ancestor_nodes_num = phylo_tree$Nnode
    terminal_nodes_num = length(phylo_tree$tip.label)
    
    ancestor_nodes = paste0("node_", 1:ancestor_nodes_num)
    terminal_nodes = phylo_tree$tip.label
    
    # In the phylo_tree$edge, 1,...,terminal_nodes_num, 
    # and then (terminal_nodes_num+1),...,(terminal_nodes_num+ancestor_nodes_num)
    # in this example, 1-1034 are ternimal_nodes, 1035-2067 are ancestor nodes
    edge = phylo_tree$edge
    
    # Change names in the edge, to avoid potential mistakes due to index and numeric value
    node = c(terminal_nodes, ancestor_nodes)
    edge_label = data.frame(from = node[edge[,1]],
                            to = node[edge[,2]])
    
    # Using igraph is faster
    library(igraph)
    g = graph_from_data_frame(edge_label, directed = TRUE)
    g_rev = reverse_edges(g)
    
    ancestor_matrix <- matrix(0, nrow = length(terminal_nodes), ncol = length(ancestor_nodes), 
                              dimnames = list(terminal_nodes, ancestor_nodes))
    
    for (term in terminal_nodes) {
        ancestors <- neighborhood(g_rev, order = vcount(g), nodes = term, mode = "out")[[1]]
        ancestor_names <- V(g)$name[ancestors]
        ancestor_matrix[term, ancestor_names[ancestor_names %in% ancestor_nodes]] <- 1
    }
    
    return(ancestor_matrix)
}

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

################################################################################
### Function: Performing chi squared test for ancestor nodes in a phylo tree ###
################################################################################

run_chi_squared_test = function(tree, x){
    
    edge = tree$edge
    x_t = t(x)
    sub_node = c((length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode))
    
    ### performing chi-squared test
    res_obs = NULL
    for(i in 1:length(sub_node)){
        edge_sub = edge[edge[,1] == sub_node[i],]
        contingency_table = matrix(c(x_t[edge_sub[1,2], ], x_t[edge_sub[2,2], ]), nrow = 2, byrow = TRUE)
        contingency_table = contingency_table[,apply(contingency_table, 2, sum) != 0,drop=FALSE]
        
        if(sum(contingency_table[1,]) > 1 & sum(contingency_table[2,]) > 1) {
            if(ncol(contingency_table) == 2){
                fit = fisher.test(contingency_table)
                res_obs = rbind(res_obs, data.frame(node = paste0("node_", sub_node[i]),
                                                    pval = fit$p.value))
            }
            if(ncol(contingency_table) > 2){
                fit = chisq.test(contingency_table)
                res_obs = rbind(res_obs, data.frame(node = paste0("node_", sub_node[i]),
                                                    pval = chi_sq_test$p.value))
            }
        }
    }
    rownames(res_obs) = NULL
    
    res_obs = res_obs[!is.nan(res_obs$pval),]
    res_obs$fdr = p.adjust(res_obs$pval, method = "fdr")
    
    return(res_obs)
}
