# 2024-5-14
# Perform gene set enrichment analysis using Reactome for each dataset and create heatmap


################################################################################################
######################## 1. Load in necessary data and packages
BiocManager::install("ReactomePA", force = TRUE)

library(clusterProfiler)
library(enrichplot)
library(ggplot2)
require(DOSE)
library(ReactomePA)
library(UpSetR)
library(ComplexHeatmap)
library(pheatmap)

# SET THE DESIRED ORGANISM HERE
organism = library('org.Mm.eg.db', character.only = TRUE)


setwd("C:/Users/livis/Documents/Benayoun_Lab/11-5_rerun/deseq")
#read in DEseq all gene statistics
BMDM_Cohort1_ <- read.csv("2024-11-05_Benayoun_BMDM_Cohort1_SEX_DIM_all_genes_statistics.txt"  , header = T, sep = "\t", row.names = NULL)
BMDM_Cohort2_ <- read.csv("2024-11-05_Benayoun_BMDM_Cohort2_SEX_DIM_all_genes_statistics.txt"  , header = T, sep = "\t", row.names = NULL)
Peritoneal_Cohort1_ <- read.csv("2024-11-05_Benayoun_Peritoneal_Cohort1_SEX_DIM_all_genes_statistics.txt", header = T, sep = "\t", row.names = NULL)
Peritoneal_Cohort2_ <- read.csv("2024-11-05_Benayoun_Peritoneal_Cohort2_SEX_DIM_all_genes_statistics.txt"  , header = T, sep = "\t", row.names = NULL)
GSE41879_Peritoneal_ <- read.csv("2024-11-05_GSE41879_Peritoneal_SEX_DIM_all_genes_statistics.txt"  , header = T, sep = "\t", row.names = NULL)
GSE99622_Microglia_whole_brain_P14_ <- read.csv("2024-11-05_GSE99622_Microglia_whole_brain_P14_SEX_DIM_all_genes_statistics.txt"  , header = T, sep = "\t", row.names = NULL)
GSE109099_Alveolar_ <- read.csv("2024-11-05_GSE109099_Alveolar_SEX_DIM_all_genes_statistics.txt", header = T, sep = "\t", row.names = NULL)
GSE124829_Microglia_whole_brain_ImmGen_8w_ <- read.csv("2024-11-05_GSE124829_Microglia_whole_brain_ImmGen_8w_SEX_DIM_all_genes_statistics.txt"  , header = T, sep = "\t", row.names = NULL)
GSE124829_Peritoneal_ImmGen_8w_ <- read.csv("2024-11-05_GSE124829_Peritoneal_ImmGen_6w_SEX_DIM_all_genes_statistics.txt"  , header = T, sep = "\t", row.names = NULL)
GSE153299_OCP_MCSF_RANKL_ <- read.csv("2024-11-05_GSE153299_OCP_MCSF_RANKL_SEX_DIM_all_genes_statistics.txt"  , header = T, sep = "\t", row.names = NULL)
GSE153299_OCP_MCSF_ <- read.csv("2024-11-05_GSE153299_OCP_MCSF_SEX_DIM_all_genes_statistics.txt"  , header = T, sep = "\t", row.names = NULL)
GSE174207_Alveolar_CC_ <- read.csv("2024-11-05_GSE174207_Alveolar_CC_SEX_DIM_all_genes_statistics.txt"  , header = T, sep = "\t", row.names = NULL)
GSE149014_Peritoneal_ <- read.csv("2024-11-05_GSE149014_Peritoneal_SEX_DIM_all_genes_statistics.txt"  , header = T, sep = "\t", row.names = NULL)
PRJNA383777_Microglia_whole_brain_ <- read.csv("2024-11-05_PRJNA383777_Microglia_whole_brain_SEX_DIM_all_genes_statistics.txt"  , header = T, sep = "\t", row.names = NULL)
PRJNA408225_Microglia_frontal_lobe_ <- read.csv("2024-11-05_PRJNA408225_Microglia_frontal_lobe_SEX_DIM_all_genes_statistics.txt"  , header = T, sep = "\t", row.names = NULL)
GSE156799_Alveolar_ <- read.csv("2024-11-05_GSE156799_Alveolar_SEX_DIM_all_genes_statistics.txt"  , header = T, sep = "\t", row.names = NULL)
GSE109099_Exudate_ <- read.csv("2024-11-05_GSE109099_Exudate_SEX_DIM_all_genes_statistics.txt"  , header = T, sep = "\t", row.names = NULL)
PRJNA408225_Microglia_Hippocampus_ <- read.csv("2024-11-05_PRJNA408225_microglia_Hippocampus_SEX_DIM_all_genes_statistics.txt"  , header = T, sep = "\t", row.names = NULL)

#put into list
my_deseq_list <- list("BMDM_Cohort1" = BMDM_Cohort1_,
                      "BMDM_Cohort2" = BMDM_Cohort2_,
                      "Peritoneal_Cohort1" = Peritoneal_Cohort1_,
                      "Peritoneal_Cohort2" = Peritoneal_Cohort2_,
                      "GSE41879_Peritoneal" = GSE41879_Peritoneal_,
                      "GSE99622_Microglia_whole_brain_P14" = GSE99622_Microglia_whole_brain_P14_,
                      "GSE109099_Alveolar" = GSE109099_Alveolar_,
                      "GSE124829_Microglia_whole_brain_ImmGen" = GSE124829_Microglia_whole_brain_ImmGen_8w_,
                      "GSE124829_Peritoneal_ImmGen" = GSE124829_Peritoneal_ImmGen_8w_,
                      "GSE153299_OCP_MCSF_RANKL" = GSE153299_OCP_MCSF_RANKL_,
                      "GSE153299_OCP_MCSF" = GSE153299_OCP_MCSF_,
                      "GSE174207_Alveolar_CC" = GSE174207_Alveolar_CC_,
                      "GSE149014_Peritoneal" = GSE149014_Peritoneal_,
                      "PRJNA383777_Microglia_whole_brain" = PRJNA383777_Microglia_whole_brain_,
                      "PRJNA408225_Microglia_frontal_lobe" = PRJNA408225_Microglia_frontal_lobe_,
                      "PRJNA408225_Microglia_Hippocampus" = PRJNA408225_Microglia_Hippocampus_,
                      "GSE156799_Alveolar" = GSE156799_Alveolar_,
                      "GSE109099_Exudate" = GSE109099_Exudate_)


setwd("C:/Users/livis/Documents/Benayoun_Lab/11-5_rerun/gsea-reactome")

################################################################################################
######################## 2. Perform gene set enrichment analysis for each dataset in list

#this runs gsea from reactome for all datasets
for(i in 1:18){
  df = my_deseq_list[[i]]
  my.outprefix <- paste(Sys.Date(), names(my_deseq_list[i]) , sep='_')
  # we want the t-statistic
  original_gene_list <- df$stat
  head(original_gene_list)
  # name the vector
  names(original_gene_list) <- df$row.names
  
  # omit any NA values 
  gene_list<-na.omit(original_gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  head(gene_list)
  
  #convert to entrezID
  entrezID.genes <- bitr(names(gene_list) , fromType="SYMBOL", toType="ENTREZID", 
                         OrgDb="org.Mm.eg.db")
  names(gene_list) <- entrezID.genes$ENTREZID
  head(gene_list)
  
  #GSEA
  reactome_results <- gsePathway(geneList=gene_list, 
                                 minGSSize = 3, 
                                 maxGSSize = 800, 
                                 pvalueCutoff = 0.05, 
                                 organism = "mouse" ,
                                 pAdjustMethod = "none")
  
  x <- as.data.frame(reactome_results)
  
  write.table(x, file = paste(my.outprefix,"GSEA_REACTOME.txt", sep ="_"), sep = "\t" , row.names = F, quote=F)
  
}

# 50 warnings about package "stats" may not be available or low p-values - ok to proceed.

################################################################################################
######################## 3. Compile GSEA results into UpSet plot and Heatmap

#read in reactome files
ra_BMDM_Cohort1_ <- read.csv("2024-11-20_BMDM_Cohort1_GSEA_REACTOME.txt"  , header = T, sep = "\t", row.names = NULL)
ra_BMDM_Cohort2_ <- read.csv("2024-11-20_BMDM_Cohort2_GSEA_REACTOME.txt"  , header = T, sep = "\t", row.names = NULL)
ra_Peritoneal_Cohort1_ <- read.csv("2024-11-20_Peritoneal_Cohort1_GSEA_REACTOME.txt", header = T, sep = "\t", row.names = NULL)
ra_Peritoneal_Cohort2_ <- read.csv("2024-11-20_Peritoneal_Cohort2_GSEA_REACTOME.txt"  , header = T, sep = "\t", row.names = NULL)
ra_GSE41879_Peritoneal_ <- read.csv("2024-11-20_GSE41879_Peritoneal_GSEA_REACTOME.txt"  , header = T, sep = "\t", row.names = NULL)
ra_GSE99622_Microglia_whole_brain_P14_ <- read.csv("2024-11-20_GSE99622_Microglia_whole_brain_P14_GSEA_REACTOME.txt"  , header = T, sep = "\t", row.names = NULL)
ra_GSE109099_Alveolar_ <- read.csv("2024-11-20_GSE109099_Alveolar_GSEA_REACTOME.txt", header = T, sep = "\t", row.names = NULL)
ra_GSE124829_Microglia_whole_brain_ImmGen_8w_ <- read.csv("2024-11-20_GSE124829_Microglia_whole_brain_ImmGen_GSEA_REACTOME.txt"  , header = T, sep = "\t", row.names = NULL)
ra_GSE124829_Peritoneal_ImmGen_8w_ <- read.csv("2024-11-20_GSE124829_Peritoneal_ImmGen_GSEA_REACTOME.txt"  , header = T, sep = "\t", row.names = NULL)
ra_GSE153299_OCP_MCSF_RANKL_ <- read.csv("2024-11-20_GSE153299_OCP_MCSF_RANKL_GSEA_REACTOME.txt"  , header = T, sep = "\t", row.names = NULL)
ra_GSE153299_OCP_MCSF_ <- read.csv("2024-11-20_GSE153299_OCP_MCSF_GSEA_REACTOME.txt"  , header = T, sep = "\t", row.names = NULL)
ra_GSE174207_Alveolar_CC_ <- read.csv("2024-11-20_GSE174207_Alveolar_CC_GSEA_REACTOME.txt"  , header = T, sep = "\t", row.names = NULL)
ra_GSE149014_Peritoneal_ <- read.csv("2024-11-20_GSE149014_Peritoneal_GSEA_REACTOME.txt"  , header = T, sep = "\t", row.names = NULL)
ra_PRJNA383777_Microglia_whole_brain_ <- read.csv("2024-11-20_PRJNA383777_Microglia_whole_brain_GSEA_REACTOME.txt"  , header = T, sep = "\t", row.names = NULL)
ra_PRJNA408225_Microglia_frontal_lobe_ <- read.csv("2024-11-20_PRJNA408225_Microglia_frontal_lobe_GSEA_REACTOME.txt"  , header = T, sep = "\t", row.names = NULL)
ra_GSE156799_Alveolar_ <- read.csv("2024-11-20_GSE156799_Alveolar_GSEA_REACTOME.txt"  , header = T, sep = "\t", row.names = NULL)
ra_GSE109099_Exudate_ <- read.csv("2024-11-20_GSE109099_Exudate_GSEA_REACTOME.txt"  , header = T, sep = "\t", row.names = NULL)
ra_PRJNA408225_Microglia_Hippocampus_           <- read.csv("2024-11-20_PRJNA408225_Microglia_Hippocampus_GSEA_REACTOME.txt"  , header = T, sep = "\t", row.names = NULL)

#make upset plot


my.exp.pathways <- list(
                        "Peritoneal_Cohort1" = ra_Peritoneal_Cohort1_$Description ,
                        "Peritoneal_Cohort2" = ra_Peritoneal_Cohort2_$Description ,  
                        "GSE124829_Peritoneal_ImmGen_8w" = ra_GSE124829_Peritoneal_ImmGen_8w_$Description ,
                        "GSE149014_Peritoneal" = ra_GSE149014_Peritoneal_$Description , 
                        "GSE41879_Peritoneal"   = ra_GSE41879_Peritoneal_$Description ,
                        "GSE99622_Microglia_whole_brain_P14" = ra_GSE99622_Microglia_whole_brain_P14_$Description ,
                        "GSE124829_Microglia_whole_brain_ImmGen_8w" = ra_GSE124829_Microglia_whole_brain_ImmGen_8w_$Description ,
                        "PRJNA383777_Microglia_whole_brain" = ra_PRJNA383777_Microglia_whole_brain_$Description ,
                        "PRJNA408225_Microglia_frontal_lobe" = ra_PRJNA408225_Microglia_frontal_lobe_$Description,
                        "PRJNA408225_Microglia_Hippocampus" = ra_PRJNA408225_Microglia_Hippocampus_$Description, 
                        "GSE109099_Alveolar" = ra_GSE109099_Alveolar_$Description ,
                        "GSE174207_Alveolar_CC" = ra_GSE174207_Alveolar_CC_$Description,
                        "GSE156799_Alveolar" = ra_GSE156799_Alveolar_$Description,
                        "GSE153299_OCP_MCSF_RANKL" = ra_GSE153299_OCP_MCSF_RANKL_$Description ,
                        "GSE153299_OCP_MCSF" = ra_GSE153299_OCP_MCSF_$Description ,
                        "BMDM_Cohort1" = ra_BMDM_Cohort1_$Description ,
                        "BMDM_Cohort2" = ra_BMDM_Cohort2_$Description,
                        "GSE109099_Exudate" = ra_GSE109099_Exudate_$Description )


my.pathways                        <- unique(unlist(my.exp.pathways))
my.exp.pathways.mat                 <- data.frame(matrix(0,length(my.pathways),length(my.exp.pathways)))
colnames(my.exp.pathways.mat)       <- names(my.exp.pathways)
rownames(my.exp.pathways.mat)       <- my.pathways



#populate upset matrix
for (j in 1:length(my.exp.pathways)) {
  my.exp.pathways.mat[my.pathways %in% my.exp.pathways[[j]],j] <- 1
}

dim(my.exp.pathways.mat)
# [1] 985   18

##write upset results to txt file
write.table(my.exp.pathways.mat,  file = paste("all_genes_gsea_reactome_upset_list", Sys.Date(), sep = "_", ".txt"), quote = F, sep = "\t")


## pdf FDR (by degree) or FDR reversed (by frequency)
pdf((paste("all_genes_gsea_reactome_upset", Sys.Date(), sep = "_", ".pdf" )), width = 25, height = 25) # degree

## upset by degree or frequency
upset(my.exp.pathways.mat, nsets = 18, nintersects = 100, order.by = "degree")

dev.off()



#sort upset list
sorted.matrix <- my.exp.pathways.mat
sorted.matrix$total <- rowSums(sorted.matrix)
sorted.matrix <- sorted.matrix[order(-sorted.matrix$total),]


#save sorted upset list
write.table(sorted.matrix,  file = paste("all_genes_gsea_reactome_upset_list_SORTED", Sys.Date(), sep = "_", ".txt"), quote = F, sep = "\t")


#####################################
#make heatmap
top100_pathways <- sorted.matrix[(1:100),]
top100_pathways <- top100_pathways[,-19]

my.toppathways <- as.data.frame(top100_pathways)


#populate new dataframe with enrichment score from reactome gsea output
for(i in 1:18){
  #start by creating the string for the name of the dataframe (use whatever you entered them as above)
  dataset.name <- colnames(my.toppathways[i])
  dataset.name <- paste("ra", dataset.name, sep="_")
  dataset.name <- paste(dataset.name, "_", sep = "")
  df <- get(dataset.name)
  #populate the dataframe
  for(j in 1:100){
    if(top100_pathways[j,i] == 0){
      next
    }
    row.name <- rownames(my.toppathways[j,])
    my.row <- which(df==row.name, arr.ind=TRUE)
    #may need to change this depending on which column contains the enrichment score
    my.toppathways[j,i] <- df[my.row[1], 5]
  }
}


my.topmatrix <- data.matrix(my.toppathways)

my.topmatrix <- my.topmatrix[(1:50),]

my.heatmap.out <- paste(Sys.Date(),"reactome_GSEA_heatmap_t100.pdf", sep = "_")
pdf(my.heatmap.out, onefile = F, width=15, height = 18)

my.heatmap.title <- "Reactome GSEA Heatmap"
pheatmap(my.topmatrix,
         cluster_cols = F,
         cluster_rows = F,
         display_numbers = T,
         fontsize_number = 3,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = T, scale="none",
         border_color = NA,angle_col = '45',
         main = my.heatmap.title, cellwidth = 10, cellheight = 10)

dev.off()

############################
#top30 heatmap
my.topmatrix <- data.matrix(my.toppathways[ 1:30,])

my.heatmap.out <- paste(Sys.Date(),"reactome_GSEA_heatmap_t30.pdf", sep = "_")
pdf(my.heatmap.out, onefile = F, width=15, height = 18)
my.heatmap.title <- "Reactome GSEA Heatmap"
pheatmap(my.topmatrix,
         cluster_cols = F,
         cluster_rows = T,
         display_numbers = F,
         fontsize_number = 5,
         colorRampPalette(rev(c("#FF1493","#FF69B4","#FFB6C1","white","#ADD8E6", "#87CEFA", "#00BFFF")))(50),
         border_color = NA,
         angle_col = '45',
         show_rownames = T, scale="none",
         main = my.heatmap.title, cellwidth = 15, cellheight = 15)

dev.off()

#######################
sink(paste0(Sys.Date(),"SMAC_RCT_Enrichment_analysis_sessionInfo.txt"))
sessionInfo()
sink()

