# 2024-5-14
# Perform gene set enrichment analysis for each dataset and create heatmap


################################################################################################
######################## 1. Load in necessary data and packages

library(clusterProfiler)
library(enrichplot)
library(ggplot2)
require(DOSE)

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


setwd("C:/Users/livis/Documents/Benayoun_Lab/11-5_rerun/gsea")

################################################################################################
######################## 2. Perform gene set enrichment analysis for each dataset in list


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
  
  #GSEA
  gse <- gseGO(geneList=gene_list, 
               ont ="ALL", 
               keyType = "SYMBOL", 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = org.Mm.eg.db, 
               pAdjustMethod = "none")
  
  x <- as.data.frame(gse)
  
  write.table(x, file = paste(my.outprefix,"GSEA_ALL_sva.txt", sep ="_"), sep = "\t" , row.names = F, quote=F)
  
  
}

# 50 warnings about package "stats" may not be available or low p-values - ok to proceed

################################################################################################
######################## 3. Compile GSEA results into UpSet plot and Heatmap


setwd("C:/Users/livis/Documents/Benayoun_Lab/11-5_rerun/gsea")

##########
#read in gsea files
gsea_BMDM_Cohort1_                                <- read.csv("2024-11-20_BMDM_Cohort1_GSEA_ALL_sva.txt"  , header = T, sep = "\t", row.names = NULL)
gsea_BMDM_Cohort2_                                <- read.csv("2024-11-20_BMDM_Cohort2_GSEA_ALL_sva.txt"  , header = T, sep = "\t", row.names = NULL)
gsea_Peritoneal_Cohort1_                          <- read.csv("2024-11-20_Peritoneal_Cohort1_GSEA_ALL_sva.txt"  , header = T, sep = "\t", row.names = NULL)
gsea_Peritoneal_Cohort2_                          <- read.csv("2024-11-20_Peritoneal_Cohort2_GSEA_ALL_sva.txt"  , header = T, sep = "\t", row.names = NULL)
gsea_GSE41879_Peritoneal_                         <- read.csv("2024-11-20_GSE41879_Peritoneal_GSEA_ALL_sva.txt"  , header = T, sep = "\t", row.names = NULL)
gsea_GSE99622_Microglia_whole_brain_P14_          <- read.csv("2024-11-20_GSE99622_Microglia_whole_brain_P14_GSEA_ALL_sva.txt"  , header = T, sep = "\t", row.names = NULL)
gsea_GSE109099_Alveolar_                          <- read.csv("2024-11-20_GSE109099_Alveolar_GSEA_ALL_sva.txt"  , header = T, sep = "\t", row.names = NULL)
gsea_GSE124829_Microglia_whole_brain_ImmGen_8w_   <- read.csv("2024-11-20_GSE124829_Microglia_whole_brain_ImmGen_GSEA_ALL_sva.txt"  , header = T, sep = "\t", row.names = NULL)
gsea_GSE124829_Peritoneal_ImmGen_8w_              <- read.csv("2024-11-20_GSE124829_Peritoneal_ImmGen_GSEA_ALL_sva.txt"  , header = T, sep = "\t", row.names = NULL)
gsea_GSE153299_OCP_MCSF_RANKL_                    <- read.csv("2024-11-20_GSE153299_OCP_MCSF_RANKL_GSEA_ALL_sva.txt"  , header = T, sep = "\t", row.names = NULL)
gsea_GSE153299_OCP_MCSF_                          <- read.csv("2024-11-20_GSE153299_OCP_MCSF_GSEA_ALL_sva.txt"  , header = T, sep = "\t", row.names = NULL)
gsea_GSE174207_Alveolar_CC_                       <- read.csv("2024-11-20_GSE174207_Alveolar_CC_GSEA_ALL_sva.txt"  , header = T, sep = "\t", row.names = NULL)
gsea_GSE149014_Peritoneal_                        <- read.csv("2024-11-20_GSE149014_Peritoneal_GSEA_ALL_sva.txt"  , header = T, sep = "\t", row.names = NULL)
gsea_PRJNA383777_Microglia_whole_brain_           <- read.csv("2024-11-20_PRJNA383777_Microglia_whole_brain_GSEA_ALL_sva.txt"  , header = T, sep = "\t", row.names = NULL)
gsea_PRJNA408225_Microglia_frontal_lobe_          <- read.csv("2024-11-20_PRJNA408225_Microglia_frontal_lobe_GSEA_ALL_sva.txt"  , header = T, sep = "\t", row.names = NULL)
gsea_GSE156799_Alveolar_                          <- read.csv("2024-11-20_GSE156799_Alveolar_GSEA_ALL_sva.txt"  , header = T, sep = "\t", row.names = NULL)
gsea_GSE109099_Exudate_                           <- read.csv("2024-11-20_GSE109099_Exudate_GSEA_ALL_sva.txt"  , header = T, sep = "\t", row.names = NULL)
gsea_PRJNA408225_Microglia_Hippocampus_           <- read.csv("2024-11-20_PRJNA408225_Microglia_Hippocampus_GSEA_ALL_sva.txt"  , header = T, sep = "\t", row.names = NULL)


#make upset plot

my.exp.pathways <- list(
  "Peritoneal_Cohort1_" = gsea_Peritoneal_Cohort1_$Description ,
  "Peritoneal_Cohort2_" = gsea_Peritoneal_Cohort2_$Description ,  
  "GSE124829_Peritoneal_ImmGen_8w_" = gsea_GSE124829_Peritoneal_ImmGen_8w_$Description ,
  "GSE149014_Peritoneal_" = gsea_GSE149014_Peritoneal_$Description , 
  "GSE41879_Peritoneal_"   = gsea_GSE41879_Peritoneal_$Description ,
  "GSE99622_Microglia_whole_brain_P14_" = gsea_GSE99622_Microglia_whole_brain_P14_$Description ,
  "GSE124829_Microglia_whole_brain_ImmGen_8w_" = gsea_GSE124829_Microglia_whole_brain_ImmGen_8w_$Description ,
  "PRJNA383777_Microglia_whole_brain_" = gsea_PRJNA383777_Microglia_whole_brain_$Description ,
  "PRJNA408225_Microglia_frontal_lobe_" = gsea_PRJNA408225_Microglia_frontal_lobe_$Description,
  "PRJNA408225_Microglia_Hippocampus_" = gsea_PRJNA408225_Microglia_Hippocampus_$Description,
  "GSE109099_Alveolar_" = gsea_GSE109099_Alveolar_$Description ,
  "GSE174207_Alveolar_CC_" = gsea_GSE174207_Alveolar_CC_$Description,
  "GSE156799_Alveolar_" = gsea_GSE156799_Alveolar_$Description,
  "GSE153299_OCP_MCSF_RANKL_" = gsea_GSE153299_OCP_MCSF_RANKL_$Description ,
  "GSE153299_OCP_MCSF_" = gsea_GSE153299_OCP_MCSF_$Description ,
  "BMDM_Cohort1_" = gsea_BMDM_Cohort1_$Description ,
  "BMDM_Cohort2_" = gsea_BMDM_Cohort2_$Description ,
  "GSE109099_Exudate_" = gsea_GSE109099_Exudate_$Description 
)



my.pathways                        <- unique(unlist(my.exp.pathways))
my.exp.pathways.mat                 <- data.frame(matrix(0,length(my.pathways),length(my.exp.pathways)))
colnames(my.exp.pathways.mat)       <- names(my.exp.pathways)
rownames(my.exp.pathways.mat)       <- my.pathways




for (j in 1:length(my.exp.pathways)) {
  my.exp.pathways.mat[my.pathways %in% my.exp.pathways[[j]],j] <- 1
}

dim(my.exp.pathways.mat)
# [1] 9378   18


# update # of sets
sum(rowSums(my.exp.pathways.mat) == 18)
#0

#0 06042024
my.pathways.expressed <- my.pathways[rowSums(my.exp.pathways.mat) == 18]
print(my.pathways.expressed, quote = F)



##write upset results to txt file
write.table(my.exp.pathways.mat,  file = paste("all_genes_gsea_upset_list", Sys.Date(), sep = "_", ".txt"), quote = F, sep = "\t")




#######################
#gsea heatmap

sorted.matrix <- my.exp.pathways.mat
sorted.matrix$total <- rowSums(sorted.matrix)
sorted.matrix <- sorted.matrix[order(-sorted.matrix$total),]
top.pathways <- sorted.matrix[,-19]

#top.pathways <- read.csv("all_genes_gsea_upset_list_2024-11-06_.txt" , nrows=100, header = T, sep = "\t" )
#rownames(top.pathways) = top.pathways[,1]
#top.pathways <- top.pathways[,-1]
#top.pathways <- top.pathways[,-18]

my.toppathways <- top.pathways


#populate new dataframe with enrichment score
for(i in 1:18){
  dataset.name <- colnames(top.pathways[i])
  dataset.name <- paste("gsea", dataset.name, sep="_")
  df <- get(dataset.name)
  for(j in 1:100){
    if(top.pathways[j,i] == 0){
      next
    }
    row.name <- rownames(top.pathways[j,])
    my.row <- which(df==row.name, arr.ind=TRUE)
    my.toppathways[j,i] <- df[my.row[1], 6]
  }
}

colnames(my.toppathways) <- colnames(top.pathways)

my.topmatrix <- data.matrix(my.toppathways)
my.topmatrix <- my.topmatrix[1:100,]

library(pheatmap)

my.heatmap.out <- paste(Sys.Date(),"t100_GSEA_heatmap.pdf", sep = "_")
pdf(my.heatmap.out, onefile = F, width=7, height = 25)
my.heatmap.title <- "GSEA Heatmap"
pheatmap(my.topmatrix,
         cluster_cols = F,
         cluster_rows = T,
         display_numbers = F,
         fontsize_number = 3,
         colorRampPalette(rev(c("#00BFFF", "#87CEFA", "#ADD8E6", "white", "#FFB6C1", "#FF69B4", "#FF1493")))(50),
         border_color = NA,
         angle_col = '45',
         show_rownames = T, scale="none",
         main = my.heatmap.title, cellwidth = 10, cellheight = 10)

dev.off()
############################
#make top 30 heatmap
my.topmatrix <- data.matrix(my.toppathways[ 1:30,])


my.heatmap.out <- paste(Sys.Date(), "t30_GSEA_heatmap.pdf", sep = "_")
pdf(my.heatmap.out, onefile = F, width=10, height = 25)
my.heatmap.title <- "GSEA Heatmap"
pheatmap(my.topmatrix,
         cluster_cols = F,
         cluster_rows = T,
         display_numbers = F,
         colorRampPalette(rev(c("#FF1493","#FF69B4","#FFB6C1","white","#ADD8E6", "#87CEFA", "#ADD8E6","#00BFFF"  )))(50),
         fontsize_number = 5,
         border_color = NA,
         angle_col = '45',
         show_rownames = T, scale="none",
         main = my.heatmap.title, cellwidth = 15, cellheight = 15)

dev.off()

#######################
sink(paste0(Sys.Date(),"SMAC_Enrichment_analysis_sessionInfo.txt"))
sessionInfo()
sink()

