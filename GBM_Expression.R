####### TCGA Expression Validation ########

library(GenomicDataCommons)
library(magrittr)
GenomicDataCommons::status()
ge_manifest = files() %>% 
  filter( ~ cases.project.project_id == 'TCGA-GBM' &
            type == 'gene_expression' &
            analysis.workflow_type == 'HTSeq - Counts') %>%
  manifest()

destdir = "/Users/amitra2/Documents/CompBioRice17/NIfTI_20140122/Expression/"
fnames = lapply(ge_manifest$id[1:174],gdcdata,
                destination_dir=destdir,overwrite=TRUE,
                progress=TRUE)

###Finding number of files in directory and names
uiud = read_xlsx("/Users/amitra2/Downloads/TCGA Clinical Data_ID.xlsx", col_names = TRUE)
uiud = data.frame(uiud$uuid)
colnames(uiud)[1] = "id"
temp = list.files(pattern = "*.gz") 
###Concatenating all files together
myfiles = lapply(temp, gunzip)

temp = list.files(pattern = "*.htseq.counts")
temp1 = gsub("\\..*","",temp)
temp1 = gsub("\"","", temp1)
include <- function (theList, toMatch){
  matches <- unique (grep(paste(toMatch), 
                          theList, value=TRUE))
  return(matches)
}
match_samples = include(temp1,uiud)
temp1 = lapply(list.files(pattern = temp))
               
###Seperate solution from above; creating multiple dataframes for storing different files
for (i in 1:length(temp)) assign(temp[i], gunzip(temp[i]))
myfiles = lapply(temp, read.delim)

dat1 = 'fffc1088-c5a6-46a0-b050-860184f6ded2.htseq.counts'
datfile = read.delim(dat1, sep = "\t")

for (i in 3 : (length(temp)+1))
  datfile[ , i] <- c(myfiles[[i-1]][[2]])

samples_in_gbm = intersect(temp1,uiud$id)
indx = grep(paste(samples_in_gbm,collapse="|"), temp1)
samples_expr = datfile[,indx]

samples_expr$ENSG00000000003.13 = gsub("\\..*","",samples_expr$ENSG00000000003.13)
gene_id_to_gene_name = function(x) return(gene_id_name[gene_id_name$ensembl_gene_id %in% x,])

temp_genes = gene_id_to_gene_name(samples_expr$ENSG00000000003.13)
temp_genes = temp_genes[match(samples_expr$ENSG00000000003.13, temp_genes$ensembl_gene_id),]
rownames(samples_expr) = make.names(as.character(temp_genes$hgnc_symbol), unique = TRUE)
samples_expr = samples_expr[,-1]
samples_expr = subset(samples_expr, rowMeans(samples_expr) > 0.05)

# create random test condition
condition = as.factor(c(sample(c(0,1), replace=TRUE, size=9)))

library(DESeq2)
library("BiocParallel")
library("vsn")
register(SnowParam(3))

###Creating DeSeqData Object
gbm_ddsdata <- DESeqDataSetFromMatrix(samples_expr, DataFrame(condition), ~ condition) 

gbm_dds <- DESeq(gbm_ddsdata, parallel = TRUE)

#Find differential gne expression
gbm_res <- results(gbm_dds, parallel = TRUE)
summary(gbm_res)
gbm_res.05 <- results(gbm_dds, alpha=.05)
summary(gbm_res.05)
#res_trim = results(dds_trim, parallel = TRUE)

topGene <- rownames(gbm_res.05)[which.min(gbm_res$padj)]
plotCounts(gbm_dds, gene=topGene)

plotMA(res, ylim=c(-5,5))
plotMA(res, ylim=c(-5,5))
with(gbm_res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

gbm_resOrdered <- gbm_res[order(gbm_res$pvalue),]
gbm_resOrderedDF <- as.data.frame(gbm_resOrdered)[1:100, ]
write.csv(gbm_resOrderedDF, file = "GBM_DEG.csv")


# Normalizes gene expression results from DESeq2 using Log 2 Norm

gbm_nt <- normTransform(gbm_dds) # defaults to log2(x+1)
gbm_log2.norm.counts <- assay(gbm_nt)

log2.norm.counts.mad <- apply(gbm_log2.norm.counts,1,mad)

log2.norm.counts3 <- gbm_log2.norm.counts
log2.norm.counts3 <- data.frame(log2.norm.counts3)

log2.norm.counts3$MAD <- log2.norm.counts.mad
log2.norm.counts3 <- log2.norm.counts3[!is.na(log2.norm.counts3$MAD),]
log2.norm.counts3 <- log2.norm.counts3[rev(order(log2.norm.counts3$MAD)),]
log2.norm.counts3.Top5000 <- log2.norm.counts3[1:5000,]
log2.norm.counts3.Top5000 <- log2.norm.counts3.Top5000[,colnames(log2.norm.counts3.Top5000)[colnames(log2.norm.counts3.Top5000)!="MAD"]]
colnames(log2.norm.counts3.Top5000) <- colnames(gbm_log2.norm.counts)


pdf(file="/Users/amitra2/Documents/CompBioRice17/NIfTI_20140122/GBM_Top5000_DESeq2_heatmap_K7.pdf",width=12,height=13)

row_dend = hclust(dist(log2.norm.counts3.Top5000,method = "euclidean")) # row clustering
col_dend = hclust(as.dist(1-cor(log2.norm.counts3.Top5000, method="pearson"))) # column clustering
set.seed(12) # for k-means clustering 12
hm_top5000 <- Heatmap(log2.norm.counts3.Top5000, name = "Log2 normalized count",
                      clustering_distance_rows = "euclidean",
                      clustering_distance_columns = "pearson",
                      clustering_method_rows = "complete",
                      clustering_method_columns = "complete",
                      
                      cluster_columns = color_branches(col_dend, k = 5),
                      
                      km = 7, 
                      
                      
                      show_row_names = FALSE,
                      row_dend_side="left",
                      column_dend_side="top",
                      column_dend_height = unit(2, "cm"),
                      row_dend_width = unit(3, "cm")
) 

hm_top5000

dev.off()


row_order_list = row_order(hm_top5000)
col_order_list = column_order(hm_top5000)
#colnames(log2.norm.counts3.Top5000)[col_order_list$]

cl1 <- gene_id_to_gene_name(rownames(log2.norm.counts3.Top5000)[row_order_list[[1]]])$hgnc_symbol # cluster 1
cl2 <- gene_id_to_gene_name(rownames(log2.norm.counts3.Top5000)[row_order_list[[2]]])$hgnc_symbol # cluster 2
cl3 <- gene_id_to_gene_name(rownames(log2.norm.counts3.Top5000)[row_order_list[[3]]])$hgnc_symbol # cluster 3
cl4 <- rownames(log2.norm.counts3.Top5000)[row_order_list[[4]]] # cluster 4
cl5 <- rownames(log2.norm.counts3.Top5000)[row_order_list[[5]]]# cluster 5
cl6 <- rownames(log2.norm.counts3.Top5000)[row_order_list[[6]]]# cluster 5
cl7 <- rownames(log2.norm.counts3.Top5000)[row_order_list[[7]]] # cluster 7
cl4 = gene_id_to_gene_name(cl4)
cl5 = gene_id_to_gene_name(cl5)

max.cl <- max(length(cl5),length(cl4),length(cl7))

df <- data.frame(
                 cluster5=c(cl5,rep("",max.cl-length(cl5))),
                 cluster4=c(cl4,rep("",max.cl-length(cl4))),
                 cluster7=c(cl7,rep("",max.cl-length(cl7)))
)

write.table(df,"/Users/amitra2/Documents/CompBioRice17/NIfTI_20140122/DESeq2_heatmap_genes_cluster_1567.txt",quote=F,row.names=F,col.names=TRUE,sep="\t")
write.table(cl4, "/Users/amitra2/Documents/CompBioRice17/NIfTI_20140122/Cluster4.txt",quote=F,row.names=F,col.names=TRUE,sep="\t")
write.table(cl5, "/Users/amitra2/Documents/CompBioRice17/NIfTI_20140122/Cluster5.txt",quote=F,row.names=F,col.names=TRUE,sep="\t")

temp = read.delim("/Users/amitra2/Documents/CompBioRice17/NIfTI_20140122/HetPathways", header=F,check.names=F,as.is=TRUE, sep = "\t")


#### Survival #####

clinical_info = read_xlsx("/Users/amitra2/Downloads/TCGA Clinical Data_ID.xlsx", col_names = TRUE)
gbm_survival = clinical_info[c(3,18,20,22)]
colnames(gbm_survival)[1] = "PatientID"
colnames(gbm_survival)[2] = "VitalStatus"
colnames(gbm_survival)[3] = "Time"
colnames(gbm_survival)[4] = "KarnofskyScore"

gbm_survival$Time = as.numeric(gbm_survival$Time)
gbm_survival$VitalStatus = ifelse(gbm_survival$VitalStatus == "Dead", 1, 0)

gbm_fit = survfit(Surv(Time,VitalStatus) ~ KarnofskyScore,
              data = gbm_survival)
# Visualize with survminer
ggsurvplot(gbm_fit, data = gbm_survival, risk.table = TRUE)

ggsurvplot(
  gbm_fit,                     # survfit object with calculated statistics.
  data = gbm_survival,  # data used to fit survival curves. 
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  xlim = c(0,1600),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 300,     # break X axis in time intervals by 500.
  ggtheme = theme_minimal(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
)



###### Methylation Validation ######### Not enough methylation matched files to validate

# uiud = read_xlsx("/Users/amitra2/Downloads/TCGA Clinical Data_ID.xlsx", col_names = TRUE)
# uiud = data.frame(uiud[c(1,4)])
# 
# library(GenomicDataCommons)
# library(magrittr)
# GenomicDataCommons::status()
# ge_manifest = files() %>% 
#   filter( ~ cases.project.project_id == 'TCGA-GBM' &
#             type == 'methylation_beta_value') %>%
#   manifest()
# 
# samples_in_gbm = as.data.frame(samples_in_gbm)
# colnames(samples_in_gbm)[1] = "id"
# samples_in_gbm$id = as.character(samples_in_gbm$id)
# 
# 
# destdir = "/Users/amitra2/Documents/CompBioRice17/NIfTI_20140122/Methylation/"
# fnames = lapply(ge_manifest$id[1:420],gdcdata,
#                 destination_dir=destdir,overwrite=TRUE,
#                 progress=TRUE)
# 
# indx2 = grep(paste(samples_in_gbm,collapse="|"), ge_manifest$id)


##### Mutect Calls #####


####################
# MuTect Filtering 
####################

mutect = read.table("/Users/amitra2/Downloads/gdc_download_20171202_043235/d8b4cd20-9cdf-49fe-9677-72ad67a4dc07/TCGA.GBM.mutect.d8b4cd20-9cdf-49fe-9677-72ad67a4dc07.DR-7.0.somatic.maf",
                    sep = "\t", header = TRUE)


tCount = 30
tVaf = 0.01 # tVaf = 0.05 (for WES)

######################$
"NRAS" %in% unique(mutect.filtered.exonic$gene)
"NRAS" %in% unique(mutect.filtered$gene)
"NRAS" %in% unique(mutect$gene)


mutect[mutect$Hugo_Symbol %in% driver_genes,c("gene","snp129","snp138")] # check if any mutations in driver genes were filtered out as SNPs -> one sample with BRAF SNP (rs71645984), clinical significance NA

mutect.filtered <- mutect[mutect$esp6500siv2_all < 0.01 & mutect$x1kg2015aug_max < 0.01 &  mutect$snp138 %in% c(".",NA) | mutect$snp138 == "rs11554290",]



######################$

mutect <- mutect[(as.numeric(mutect[, "t_ref_count"]) + as.numeric(mutect[,"t_alt_count"])) > tCount, , drop = FALSE]
#mutect.filtered <- mutect.filtered[(as.numeric(mutect.filtered[, "n_ref_count"]) + as.numeric(mutect.filtered[,"n_alt_count"])) > nCount, , drop = FALSE]
mutect.filtered <- mutect[as.numeric(mutect[, "tumor_f"]) > tVaf, , drop = FALSE]
#normal_f <- as.numeric(mutect.filtered[, "n_alt_count"])/(as.numeric(mutect.filtered[,"n_ref_count"]) + as.numeric(mutect.filtered[, "n_alt_count"]))
mutect.filtered <- mutect.filtered[normal_f < nVaf, , drop = FALSE]
## no segdup filtering ##
#mutect.filtered <- mutect.filtered[is.na(mutect.filtered$segdup), ]
mutect.filtered <- mutect.filtered[mutect.filtered$covered=="COVERED",]


mutect.filtered.exonic <- mutect.filtered[!(mutect.filtered$exonicfunc %in% c(".","unknown")) | 
                                            mutect.filtered$func == "splicing" #| # include splice site mutations as exonic mutations for now
                                          #mutect.filtered$gene == "TERT" & mutect.filtered$func == "upstream" # include TERT promoter mutation
                                          ,]
# for functional mutation filtering

rownames(mutect.filtered.exonic) <- NULL

mutect.filtered.nonsyn.exonic <- mutect.filtered[mutect.filtered$exonicfunc %in% c("nonsynonymous SNV","stopgain"),]
rownames(mutect.filtered.nonsyn.exonic) <- NULL

mutect.filtered.syn.exonic <- mutect.filtered[mutect.filtered$exonicfunc %in% c("synonymous SNV"),]
rownames(mutect.filtered.syn.exonic) <- NULL

### number of exonic mutations per sample

n_nonsyn_per_sample <- unlist(lapply(unique(mutect.filtered.nonsyn.exonic$sample_name), 
                                     function(x) nrow(mutect.filtered.nonsyn.exonic[mutect.filtered.nonsyn.exonic$sample_name==x,]) ))

summary(n_nonsyn_per_sample)

n_syn_per_sample <- unlist(lapply(unique(mutect.filtered.syn.exonic$sample_name), 
                                  function(x) nrow(mutect.filtered.syn.exonic[mutect.filtered.syn.exonic$sample_name==x,]) ))

summary(n_syn_per_sample)

n_exonic_per_sample <- n_nonsyn_per_sample + n_syn_per_sample


#### Mutation plot ####

###### TCGA Biolinks ######

source("https://bioconductor.org/biocLite.R")
biocLite("TCGAbiolinks")

gbm.mutect.maf <- GDCquery_Maf('GBM', pipelines = "mutect2")

mut <- GDCquery_Maf("GBM", pipelines = "mutect2")
clin <- GDCquery_clinic("TCGA-GBM","clinical")
clin <- clin[,c("bcr_patient_barcode","disease","gender","tumor_stage","race","vital_status")]

TCGAvisualize_oncoprint(mut = mut, genes = mut$Hugo_Symbol[1:20],
                        filename = "oncoprint.pdf",
                        annotation = clin,
                        color=c("background"="#CCCCCC","DEL"="purple","INS"="yellow","SNP"="brown"),
                        rows.font.size=10,
                        width = 10,
                        heatmap.legend.side = "right",
                        dist.col = 0,
                        label.font.size = 10)

library(devtools)
library("GenVisR")
set.seed(426)

gbmMAF = read.table("/Users/amitra2/Downloads/gdc_download_20171202_043235/d8b4cd20-9cdf-49fe-9677-72ad67a4dc07/TCGA.GBM.mutect.d8b4cd20-9cdf-49fe-9677-72ad67a4dc07.DR-7.0.somatic.maf",
                    sep = "\t", header = TRUE)

gbmMAF$Variant_Classification = gsub("\t", " ", gbmMAF$Variant_Classification)
gbmMAF$Variant_Classification = sub(" .*", "", gbmMAF$Variant_Classification)

gbmMAF = subset(gbmMAF, subset =  Variant_Classification %in% c("Missense_Mutation","Silent","RNA","Nonsense_mutation", 
                                                              "Frame_Shift_Del","Splice_Site", "In_Frame_Ins", "In_Frame_Del", "Frame_Shift_Ins", "Nonstop_Mutation"))
gbmMAF$Variant_Classification = as.factor(gbmMAF$Variant_Classification)
# Plot with the MAF file type specified (default) The mainRecurCutoff
# parameter is described in the next section
Waterfall(gbmMAF, fileType = "MAF", mainRecurCutoff = 0.05)

# Plot the mutation landscape
aterfall(gbmMAF, fileType="MAF", mainRecurCutoff = 0.10)
# Alter the coverage space to whole genome space
waterfall(gbmMAF, fileType="MAF", mainRecurCutoff = 0.05, maxGenes = 25, coverageSpace = 3.2e+09)

genes_to_plot = c("TP53","PTEN","RB1","IDH1","JAK","STAT","PDGFRA","TNF1")
waterfall(gbmMAF, fileType="MAF", plotGenes = genes_to_plot)

# Create a data frame specifying the mutation burden for each sample
tumor_sample = unique(gbmMAF$Tumor_Sample_Barcode)
mutation_burden = sample(1:387, length(tumor_sample), replace = TRUE)
mutation_rate = data.frame(sample = tumor_sample, mut_burden = mutation_burden)

# Alter the coverage space to whole genome space
waterfall(gbmMAF, mutBurden = mutation_rate, mainRecurCutoff = 0.05, maxGenes = 25)
