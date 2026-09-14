
#################################################################
### Identifying differential expressed genes in the founder cells

###########################
### single cell level

work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape"
source("~/work/scripts/utils.R")
library(Signac)
library(DESeq2)

mouse_gene = read.table("~/work/tome/code/mouse.v23.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

obj = readRDS(paste0(work_path, "/revision/data_10X_multiome/obj_combine_processed.rds"))
mat = GetAssayData(obj, assay = "RNA", layer = "counts")

cell_assignment = readRDS(paste0(work_path, "/revision/data_10X_multiome/cell_well_cluster_assignment.rds"))
dat_1 = data.frame(well_cluster = paste0("cluster_", c(8,16,18,20)), category = "PSC")
dat_2 = data.frame(well_cluster = paste0("cluster_", c(3,4,5,21,22)), category = "Neuroectodermal")
dat_3 = data.frame(well_cluster = paste0("cluster_", c(9,14,15,19)), category = "Mesodermal")
dat_4 = data.frame(well_cluster = paste0("cluster_", c(10,11,17,2,23,24,25,26,27,28,29,6,7)), category = "Mixed")
pd = rbind(dat_1, dat_2, dat_3, dat_4)
cell_assignment = cell_assignment %>% left_join(pd, by = "well_cluster") %>% filter(!is.na(category))
rownames(cell_assignment) = cell_assignment$cell_id

mat = mat[,as.vector(cell_assignment$cell_id)]
mat = mat[Matrix::rowSums(mat) >= 10,]
obj = CreateSeuratObject(mat, meta.data = cell_assignment[,c("well_cluster", "category")])
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
Idents(obj) = as.vector(obj$category)

res = NULL
for(category in c("PSC", "Neuroectodermal", "Mesodermal")){
    res_x = FindMarkers(obj,
                  ident.1 = category,
                  ident.2 = "Mixed",
                  logfc.threshold = 0,
                  min.pct = 0,
                  return.thresh = 1,
                  only.pos = FALSE)
    res_x$gene = rownames(res_x)
    res_x$category = category
    res = rbind(res, res_x)
}
res = res %>% left_join(mouse_gene %>% select(gene = gene_ID, gene_short_name), by = "gene")
print(res[res$p_val_adj < 0.05,])


res_perm_list = list()
for(i in 1:10){
    print(paste0("perm: ", i))
    pd_perm = pd
    set.seed(i)
    pd_perm$category = sample(as.vector(pd$category))
    obj_perm = obj
    obj_perm_x = data.frame(obj_perm[[]]) %>% select(-category) %>% left_join(pd_perm, by = "well_cluster")
    Idents(obj_perm) = sample(as.vector(obj_perm_x$category))
    res_x = NULL
    for(category in c("PSC", "Neuroectodermal", "Mesodermal")){
        res_xx = FindMarkers(obj_perm,
                      ident.1 = category,
                      ident.2 = "Mixed",
                      logfc.threshold = 0,
                      min.pct = 0,
                      return.thresh = 1,
                      only.pos = FALSE)
        res_xx$gene = rownames(res_xx)
        res_xx$category = category
        res_x = rbind(res_x, res_xx)
    }
    res_perm_list[[i]] = res_x %>% left_join(mouse_gene %>% select(gene = gene_ID, gene_short_name), by = "gene")
}
res_perm = do.call(rbind, res_perm_list)
print(res_perm[res_perm$p_val_adj < 0.05,])


compute_perm_fdr <- function(p_real, p_perm) {
    p_real <- p_real[!is.na(p_real)]
    p_perm <- p_perm[!is.na(p_perm)]
    n_real <- length(p_real)
    n_perm <- length(p_perm)

    # for each real p-value (as a threshold), count how many real and perm pass
    ord <- order(p_real)
    p_sorted <- p_real[ord]

    # number of real discoveries at threshold = rank
    R <- seq_along(p_sorted)
    # expected false discoveries: scale perm count to real sample size
    V <- sapply(p_sorted, function(t) sum(p_perm <= t)) * (n_real / n_perm)

    fdr_raw <- V / R
    # enforce monotonicity (FDR shouldn't decrease as threshold loosens)
    fdr <- rev(cummin(rev(fdr_raw)))
    fdr <- pmin(fdr, 1)

    out <- data.frame(p_val = p_real, fdr = NA_real_)
    out$fdr[ord] <- fdr
    out
}

res$fdr_perm <- compute_perm_fdr(res$p_val, res_perm$p_val)$fdr

x = res[res$fdr_perm < 0.05,]
write.csv(x, '~/share/DAG_perm_fdr.csv', row.names=F, quote=F)

saveRDS(res, paste0(work_path, "/revision/data_10X_multiome/DAG/single_cell_result.rds"))
saveRDS(res_perm, paste0(work_path, "/revision/data_10X_multiome/DAG/single_cell_result_perm.rds"))
### 16,605 genes

### making QQ-plot

compute_qq <- function(pvals, group_name, gene, category) {
    n <- length(pvals)
    ord <- order(pvals)
    data.frame(
        expected = -log10(ppoints(n)),
        observed = -log10(pvals[ord]),
        group    = group_name,
        gene     = gene[ord],
        category = category[ord]
    )
}


res_real = readRDS(paste0(work_path, "/revision/data_10X_multiome/DAG/single_cell_result.rds"))
res_perm = readRDS(paste0(work_path, "/revision/data_10X_multiome/DAG/single_cell_result_perm.rds"))

res_real = res_real[!is.na(res_real$p_val),]
res_perm = res_perm[!is.na(res_perm$p_val),]

qq_df <- rbind(
    compute_qq(res_real$p_val, "Real", res_real$gene_short_name, res_real$category),
    compute_qq(res_perm$p_val, "Permutation", res_perm$gene_short_name, res_perm$category)
)

x1 = min(res_real$p_val[res_real$fdr_perm > 0.05])
x2 = max(res_real$p_val[res_real$fdr_perm < 0.05])
x = -log10((x1+x2)/2)

y = NULL
y = rbind(y, data.frame(gene = "Tet1", category = "PSC"))
y = rbind(y, data.frame(gene = "Cald1", category = "PSC"))
y = rbind(y, data.frame(gene = "Meg3", category = "PSC"))
y = rbind(y, data.frame(gene = "Rian", category = "PSC"))
y = rbind(y, data.frame(gene = "Mirg", category = "PSC"))
y = rbind(y, data.frame(gene = "Peg10", category = "PSC"))
y = rbind(y, data.frame(gene = "Plagl1", category = "PSC"))
y = rbind(y, data.frame(gene = "Nedd4", category = "Neuroectodermal"))
y = rbind(y, data.frame(gene = "Meis2", category = "Mesodermal"))


library(ggrepel)

label_df <- merge(y, qq_df, by = c("gene", "category")) %>% filter(group == "Real")


p = ggplot(qq_df, aes(x = expected, y = observed, color = group)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("Real" = "#E64B35", "Permutation" = "#4DBBD5")) +
    labs(x = expression(-log[10]("Expected"~italic(P))),
         y = expression(-log[10]("Observed"~italic(P))),
         color = NULL) +
    theme_classic(base_size = 12) +
    geom_hline(yintercept = x) +
    theme(legend.position = c(0.2, 0.85))

p = p + geom_text_repel(
  data = label_df,
  aes(label = gene),
  size = 6,
  color = "black",
  family = "Arial",
  fontface = "bold.italic",
  max.overlaps = 20,
  show.legend = FALSE,
  segment.color = "black",
  segment.size = 0.4,
  min.segment.length = 0,
  box.padding = 0.5,
  point.padding = 0.3,
  nudge_x = 1,
  nudge_y = -1
)

ggsave("~/share/QQ_plot_DAG.png",p,width=5,height=5,dpi=300)


save_path = "/net/shendure/vol10/www/content/members/cxqiu/public/nobackup/sam_tape/mESCs_DAG_DAP"

write.csv(res_real, paste0(save_path, "/single_cell_DAG_real.csv"), quote=F, row.names=F)
write.csv(res_perm, paste0(save_path, "/single_cell_DAG_perm.csv"), quote=F, row.names=F)






