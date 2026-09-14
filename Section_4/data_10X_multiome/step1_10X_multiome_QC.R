
###################################
### Revision-1: 

### Experiment 1 - 10X Multiome run from 9/15:
###    scATACseq -  sequenced 4 libraries on NextSeq 2000 P2 100 cycle kit 
###    scRNAseq -  sequenced 4 libraries on NextSeq 2000 P2 100 cycle kit 
###    Aidan has demultiplexed the runs and handed off to CX for sequencing assessment.

### Notes from Aidan:

### We did 4 lanes of 10X -- here is the path to the unfiltered data from cellranger for each of the 4 lanes:
### Lane 1: /net/shendure/vol10/projects/akeith/nobackup/analyses/SR_revisions/ES_multiome/count/arc_L1/outs/raw_feature_bc_matrix/
### Lane 2: /net/shendure/vol10/projects/akeith/nobackup/analyses/SR_revisions/ES_multiome/count/arc_L2/outs/raw_feature_bc_matrix/
### Lane 3: /net/shendure/vol10/projects/akeith/nobackup/analyses/SR_revisions/ES_multiome/count/arc_L3/outs/raw_feature_bc_matrix/
### Lane 4: /net/shendure/vol10/projects/akeith/nobackup/analyses/SR_revisions/ES_multiome/count/arc_L4/outs/raw_feature_bc_matrix/
### 
### (there are also .h5 files if that is easier at <...>/outs/raw_feature_bc_matrix.h5).
### 
### I also ran cellranger aggregate to combine all of the lanes into a single output. Unfiltered data for the merged lanes can be found in the directory: /net/shendure/vol10/projects/akeith/nobackup/analyses/SR_revisions/ES_multiome/count/arc_merged/outs/
###
### Here is the script I used to run cellranger:

######################################
# cellranger-arc_count-aggr.sh 
######################################
# Description: script to process
# 4 lanes of 10X Multiome that were
# generated from barcoded mESCs for
# Sam and Choi's revision experiments
######################################

# paths to working directory and reference
work_dir=/net/shendure/vol10/projects/akeith/nobackup/analyses/SR_revisions/ES_multiome/
ref=/net/shendure/vol10/projects/akeith/nobackup/refs/refdata-cellranger-arc-mm10-2020-A

# load modules
module load cellranger-arc/2.0.2
module load bcl2fastq/2.20

# change into working directory
cd ${work_dir}/count

######################################
# Step 1
# Running cellranger-arc count
######################################
echo Running cellranger-arc count for each lane...
# run cellranger-arc count on each of the 4 lanes individually
for suffix in $(cat suffixes.tsv); do echo Running cellranger count for lane $suffix; \
cellranger-arc count --id=arc_${suffix} \
--reference=$ref \
--libraries=libraries_${suffix}.csv \
--localcores=16;
echo Finished running cellranger count on $suffix; echo ; done
echo

######################################
# Step 2 
# Running cellranger-arc aggregate 
######################################
echo Running cellranger-arc aggregate to combine lanes...
cellranger-arc aggr  --id=arc_merged \
--csv=aggr_libraries.csv \
--normalize=depth \
--reference=$ref \
--localcores=16 \
--nosecondary
echo Finished running cellranger-arc aggregate
    


###########################################
# Step 3
# Analyzing the multiome data using Signac
###########################################

source("~/work/scripts/utils.R")
work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape"
data_path = "/net/shendure/vol10/projects/akeith/nobackup/analyses/SR_revisions/ES_multiome/count/arc_merged/outs"

mouse_gene = read.table("~/work/tome/code/mouse.v23.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

############
### RNA only

counts = Read10X(paste0(data_path, "/raw_feature_bc_matrix"), gene.column = 1)
count_matrix = counts[["Gene Expression"]]

count_matrix_nonzero = count_matrix
count_matrix_nonzero@x[count_matrix_nonzero@x >= 1] = 1

pd = data.frame(cell_id = colnames(count_matrix),
                UMI_count = colSums(count_matrix),
                gene_count = colSums(count_matrix_nonzero),
                lane = unlist(lapply(colnames(count_matrix), function(x) strsplit(x,"[-]")[[1]][2])))
rownames(pd) = as.vector(pd$cell_id)
pd$lane = paste0("lane_", pd$lane)
pd$experiment_id = "10X_multiome"
pd$log2_UMI_count = log2(pd$UMI_count + 1)
pd$log2_gene_count = log2(pd$gene_count + 1)

p1 = ggplot() + geom_histogram(data = pd[pd$UMI_count > 50,], aes(x = log2_UMI_count), bins = 100) +
    geom_vline(xintercept = log2(500+1)) 
p2 = ggplot() + geom_histogram(data = pd[pd$UMI_count > 50,], aes(x = log2_gene_count), bins = 100) +
    geom_vline(xintercept = log2(250+1)) 
ggsave("~/share/hist_umi_raw.pdf", p1+p2, width = 10, height = 5)

pd$log2_UMI_count = pd$log2_gene_count = NULL

mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "X", "Y", "M")),]

count_matrix_filter = count_matrix[rownames(count_matrix) %in% as.vector(mouse_gene_sub$gene_ID),
                                   pd$UMI_count >= 500 & pd$gene_count >= 250]
pd_filter = pd[pd$UMI_count >= 500 & pd$gene_count >= 250,]
pd_filter$log2_umi = log2(pd_filter$UMI_count)

MT_gene = as.vector(mouse_gene[grep("^mt-",mouse_gene$gene_short_name),]$gene_ID)
Rpl_gene = as.vector(mouse_gene[grep("^Rpl",mouse_gene$gene_short_name),]$gene_ID)
Mrpl_gene = as.vector(mouse_gene[grep("^Mrpl",mouse_gene$gene_short_name),]$gene_ID)
Rps_gene = as.vector(mouse_gene[grep("^Rps",mouse_gene$gene_short_name),]$gene_ID)
Mrps_gene = as.vector(mouse_gene[grep("^Mrps",mouse_gene$gene_short_name),]$gene_ID)
RIBO_gene = c(Rpl_gene, Mrpl_gene, Rps_gene, Mrps_gene)

pd_filter$MT_pct = 100 * Matrix::colSums(count_matrix_filter[rownames(count_matrix_filter) %in% MT_gene,])/Matrix::colSums(count_matrix_filter)
pd_filter$RIBO_pct = 100 * Matrix::colSums(count_matrix_filter[rownames(count_matrix_filter) %in% RIBO_gene,])/Matrix::colSums(count_matrix_filter)

Matrix::writeMM(t(count_matrix_filter), paste0(work_path, "/revision/data_10X_multiome/gene_count.mtx"))
write.csv(pd_filter, paste0(work_path, "/revision/data_10X_multiome/df_cell.csv"))
write.csv(mouse_gene[rownames(count_matrix_filter),], paste0(work_path, "/revision/data_10X_multiome/df_gene.csv"))

### running scrublet to identify doublets
### python Scrublet.py

scrublet_score = read.csv(paste0(work_path, "/revision/data_10X_multiome/doublet_scores_observed_cells.csv"), header=F)
pd_filter$doublet_score = as.vector(scrublet_score$V1)

simulated_scrublet_score = read.csv(paste0(work_path, "/revision/data_10X_multiome/doublet_scores_simulated_doublets.csv"), header=F)
df = data.frame(score = c(as.vector(scrublet_score$V1), as.vector(simulated_scrublet_score$V1)),
                group = rep(c("observed", "simulated"), times = c(nrow(scrublet_score), nrow(simulated_scrublet_score))))
p1 = ggplot(df, aes(score, fill = group)) + geom_histogram(alpha=0.6, binwidth = 0.02)
p2 = ggplot() + geom_histogram(data = pd_filter, aes(x = log2_umi), bins = 100) 
p3 = ggplot() + geom_histogram(data = pd_filter, aes(x = MT_pct), bins = 100) +
    geom_vline(xintercept = 2) 
ggsave("~/share/hist_umi.pdf", p1+p2+p3, width = 12, height = 5)

obj = CreateSeuratObject(count_matrix_filter, meta.data = pd_filter)
saveRDS(obj, paste0(work_path, "/revision/data_10X_multiome/obj_RNA_raw.rds"))

### filtering cells 
pd_filter = pd_filter[pd_filter$log2_umi >= 0 &
                      pd_filter$UMI_count <= quantile(pd_filter$UMI_count, 0.995) &
                      pd_filter$doublet_score <= 0.2 &
                      pd_filter$MT_pct <= 20 &
                      pd_filter$MT_pct >= 2 &
                      pd_filter$RIBO_pct <= 40,]

count_matrix_filter = count_matrix_filter[,as.vector(pd_filter$cell_id)]

obj = CreateSeuratObject(count_matrix_filter, meta.data = pd_filter)
saveRDS(obj, paste0(work_path, "/revision/data_10X_multiome/obj_RNA.rds"))

## rm gene_count.mtx
## rm df_cell.csv
## rm df_gene.csv
## rm doublet_scores_observed_cells.csv 
## rm doublet_scores_simulated_doublets.csv



#############
### ATAC only

######################################
### manually calculating FRiP and FRiT

bedtools intersect -a /net/shendure/vol10/projects/akeith/nobackup/analyses/SR_revisions/ES_multiome/count/arc_merged/outs/atac_fragments.tsv.gz \
-b /net/shendure/vol2/projects/cxqiu/JAX_atac/Novaseq/analysis_mouse/precheck/tss_file_1Kb.sorted.bed \
-u > /net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape/revision/data_10X_multiome/atac_fragments_TSS_1Kb.tsv
bgzip -f /net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape/revision/data_10X_multiome/atac_fragments_TSS_1Kb.tsv
tabix -p bed /net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape/revision/data_10X_multiome/atac_fragments_TSS_1Kb.tsv.gz

bedtools intersect -a /net/shendure/vol10/projects/akeith/nobackup/analyses/SR_revisions/ES_multiome/count/arc_merged/outs/atac_fragments.tsv.gz \
-b /net/shendure/vol10/projects/akeith/nobackup/analyses/SR_revisions/ES_multiome/count/arc_merged/outs/atac_peaks.bed \
-u > /net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape/revision/data_10X_multiome/atac_fragments_peaks.tsv
bgzip -f /net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape/revision/data_10X_multiome/atac_fragments_peaks.tsv
tabix -p bed /net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape/revision/data_10X_multiome/atac_fragments_peaks.tsv.gz


################################
#### filtering cells in ATAC-seq

library(Signac)

counts = Read10X(paste0(data_path, "/raw_feature_bc_matrix"), gene.column = 1)
counts = counts[["Peaks"]]

frag_info = CountFragments(paste0(data_path, "/atac_fragments.tsv.gz"))
FRIP = CountFragments(paste0(work_path, "/revision/data_10X_multiome/atac_fragments_peaks.tsv.gz"))
FRIP = FRIP[,c(1,2,5)]; colnames(FRIP) = c("CB", "Freq_Peaks", "Reads_Peaks")
FRIT = CountFragments(paste0(work_path, "/revision/data_10X_multiome/atac_fragments_TSS_1Kb.tsv.gz"))
FRIT = FRIT[,c(1,2,5)]; colnames(FRIT) = c("CB", "Freq_TSSs", "Reads_TSSs")
frag_info = frag_info %>% left_join(FRIP, by = "CB") %>% left_join(FRIT, by = "CB") %>% as.data.frame()
saveRDS(frag_info, paste0(work_path, "/revision/data_10X_multiome/ATAC_frag_info.rds"))

frag_info$pct_reads_in_peaks = 100*frag_info$Reads_Peaks/frag_info$reads_count
frag_info$pct_reads_in_tss = 100*frag_info$Reads_TSSs/frag_info$reads_count
frag_info$cell_id = frag_info$CB
frag_info$pct_reads_in_peaks[is.na(frag_info$pct_reads_in_peaks)] = 0
frag_info$pct_reads_in_tss[is.na(frag_info$pct_reads_in_tss)] = 0


chrom_assay = CreateChromatinAssay(
    counts = counts,
    sep = c(":","-"),
    genome = "mm10",
    fragments = paste0(data_path, "/atac_fragments.tsv.gz"),
    min.cells = 10,
    min.features = 200
)

metadata = data.frame(cell_id = colnames(chrom_assay),
                      lane = unlist(lapply(colnames(chrom_assay), function(x) strsplit(x,"[-]")[[1]][2])))
metadata$lane = paste0("lane_", metadata$lane)
metadata = metadata %>% left_join(frag_info[,c("cell_id","pct_reads_in_peaks","pct_reads_in_tss")]) %>% as.data.frame()
obj = CreateSeuratObject(counts=chrom_assay, assay='peaks', meta.data=metadata)    

peaks.keep = seqnames(granges(obj)) %in% standardChromosomes(granges(obj))
obj = obj[as.vector(peaks.keep), ]

annotation = readRDS("/net/shendure/vol10/projects/cxqiu/nobackup/genome/atac_data/mm10/annotations.RDS")
Annotation(obj) = annotation

obj = NucleosomeSignal(object = obj)
obj = TSSEnrichment(object = obj)
obj$blacklist_ratio = FractionCountsInRegion(
    object = obj, 
    assay = 'peaks',
    regions = blacklist_mm10
)

saveRDS(obj, paste0(work_path, "/revision/data_10X_multiome/obj_ATAC_raw.rds"))


#################################
### run AMULET to detect doublets

cell_list = data.frame(obj[[]])[,c("cell_id","lane")]
cell_list$cell_id = unlist(lapply(cell_list$cell_id, function(x) strsplit(x,"[-]")[[1]][1]))
cell_list$cell_id = paste0(cell_list$cell_id, "-1")
write.table(cell_list, paste0(work_path, "/revision/data_10X_multiome/amulet_analysis/ATAC_cell_list_raw.txt"), row.names=F, col.names=F, quote=F, sep="\t")

# AMULET.sh

amu_res = NULL
for(i in 1:4){
    amu_i = read.table(paste0(work_path, "/revision/data_10X_multiome/amulet_analysis/arc_L", i, "/MultipletCellIds_01.txt"))
    amu_i = unlist(lapply(as.vector(amu_i$V1), function(x) strsplit(x,"[-]")[[1]][1]))
    amu_res = c(amu_res, paste0(amu_i, "-", i))
}
saveRDS(amu_res, paste0(work_path, "/revision/data_10X_multiome/amulet_analysis/amulet_res.rds"))

############
### checking quality

p = DensityScatter(obj, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
ggsave("~/share/ATAC_QC_1.pdf", p, width = 6, height = 4)

p = VlnPlot(
    object = obj,
    features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks', 'pct_reads_in_tss'),
    pt.size = 0,
    ncol = 3
)
ggsave("~/share/ATAC_QC_2.pdf", p, width = 6, height = 6)

pd = data.frame(obj[[]])
pd$log2_nCount_peaks = log2(pd$nCount_peaks)
p1 = ggplot() + 
    geom_histogram(data = pd, aes(x = log2_nCount_peaks), bins = 100) + 
    theme_classic(base_size = 5)
p2 = ggplot() + 
    geom_histogram(data = pd, aes(x = pct_reads_in_peaks), bins = 100) + 
    theme_classic(base_size = 5) + 
    scale_x_continuous(breaks = seq(floor(min(pd$pct_reads_in_peaks)),ceiling(max(pd$pct_reads_in_peaks)), by = 1))
p3 = ggplot() + 
    geom_histogram(data = pd, aes(x = pct_reads_in_tss), bins = 100) + 
    theme_classic(base_size = 5) + 
    scale_x_continuous(breaks = seq(floor(min(pd$pct_reads_in_tss)),ceiling(max(pd$pct_reads_in_tss)), by = 1))
p4 = ggplot() + 
    geom_histogram(data = pd[pd$nucleosome_signal < 2,], aes(x = nucleosome_signal), bins = 100) + 
    theme_classic(base_size = 5)
p5 = ggplot() + 
    geom_histogram(data = pd[pd$TSS.enrichment < 20,], aes(x = TSS.enrichment), bins = 100) + 
    theme_classic(base_size = 5)
p6 = ggplot() + 
    geom_histogram(data = pd[pd$blacklist_ratio < 0.1,], aes(x = blacklist_ratio), bins = 100) + 
    theme_classic(base_size = 5)
ggsave("~/share/ATAC_QC_3.pdf", (p1 | p2 | p3) / (p4 | p5 | p6), width = 15, height = 6)

### filtering cells
obj$keep = pd$log2_nCount_peaks > mean(pd$log2_nCount_peaks) - sd(pd$log2_nCount_peaks) &
       pd$log2_nCount_peaks < mean(pd$log2_nCount_peaks) + 2 * sd(pd$log2_nCount_peaks) &
       pd$pct_reads_in_peaks > mean(pd$pct_reads_in_peaks) - sd(pd$pct_reads_in_peaks) &
       pd$pct_reads_in_tss > mean(pd$pct_reads_in_tss) - sd(pd$pct_reads_in_tss) &
       pd$TSS.enrichment > mean(pd$TSS.enrichment) - sd(pd$TSS.enrichment) &
       pd$nucleosome_signal < 1.5 &
       pd$blacklist_ratio < 0.05 &
       !pd$cell_id %in% amu_res

obj = subset(obj, keep)
saveRDS(obj, paste0(work_path, "/revision/data_10X_multiome/obj_ATAC.rds"))



###########################################
# Step 4
# Combining RNA and ATAC
###########################################

library(Signac)
source("~/work/scripts/utils.R")
work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape"
data_path = "/net/shendure/vol10/projects/akeith/nobackup/analyses/SR_revisions/ES_multiome/count/arc_merged/outs"

mouse_gene = read.table("~/work/tome/code/mouse.v23.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

annotation = readRDS("/net/shendure/vol10/projects/cxqiu/nobackup/genome/atac_data/mm10/annotations.RDS")

obj_rna = readRDS(paste0(work_path, "/revision/data_10X_multiome/obj_RNA.rds"))
obj_atac = readRDS(paste0(work_path, "/revision/data_10X_multiome/obj_ATAC.rds"))

count_rna = GetAssayData(obj_rna, slot = "counts")
pd_rna = data.frame(obj_rna[[]])
count_atac = GetAssayData(obj_atac, slot = "counts")
pd_atac = data.frame(obj_atac[[]])

cell_overlap = intersect(rownames(pd_rna), rownames(pd_atac))
pd = cbind(pd_rna[cell_overlap,], 
           pd_atac[cell_overlap, !colnames(pd_atac) %in% c("orig.ident","cell_id","lane","keep")])

obj = CreateSeuratObject(count_rna[,cell_overlap],
                         assay = "RNA",
                         meta.data = pd)

obj[["ATAC"]] = CreateChromatinAssay(
    counts = count_atac[,cell_overlap],
    sep = c("-","-"),
    fragments = paste0(data_path, "/atac_fragments.tsv.gz"),
    annotation = annotation
)

saveRDS(obj, paste0(work_path, "/revision/data_10X_multiome/obj_combine.rds"))

### n = 63,889 cells, 32,208 genes, and 205,458 peaks
### median UMI counts of RNA-seq is 1,169; median gene counts is 745
### median read counts of ATAC-seq is 2,667; median read counts in peaks of ATAC-seq is 2,314;
### PCR duplicate rate for RNA is ~25%; for ATAC is ~10%





#######################################
### How many mito% per cell in ATAC-seq

cd /net/shendure/vol10/projects/akeith/nobackup/analyses/SR_revisions/ES_multiome/count/arc_merged/outs
tabix atac_fragments.tsv.gz chr1 > /net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tape/revision/data_10X_multiome/atac_fragments.chr1.tsv



