# 2024-10-24
# Analysis with ECM related gene lists based on Sex

################################################################################################
######################## 1. Load in necessary data and packages

library(pheatmap)
library(phenoTest)

######################## A. Prep ECM gene Sets ########################
setwd("C:/Users/livis/Documents/Benayoun_Lab/ecm")
my.1.ECM     <- read.table('./GeneSets/CoreMatrisome.txt', sep = "\t", header = T)
my.2.ECM     <- read.table('./GeneSets/Collagens.txt', sep = "\t", header = T)
my.3.ECM     <- read.table('./GeneSets/ECMaffiliated.txt', sep = "\t", header = T)
my.4.ECM     <- read.table('./GeneSets/ECMglycoproteins.txt', sep = "\t", header = T)
my.5.ECM     <- read.table('./GeneSets/ECMregulators.txt', sep = "\t", header = T)
my.6.ECM     <- read.table('./GeneSets/Proteoglycans.txt', sep = "\t", header = T)
my.7.ECM     <- read.table('./GeneSets/Matrisomeassociated.txt', sep = "\t", header = T)
my.8.ECM     <- read.table('./GeneSets/Secretedfactors.txt', sep = "\t", header = T)

my.ECM.curated.gs <- list("CoreMatrisome"       = my.1.ECM$GeneSymbol,
                          "Collagens"           = my.2.ECM$GeneSymbol,
                          "ECMaffiliated"       = my.3.ECM$GeneSymbol,
                          "ECMglycoproteins"    = my.4.ECM$GeneSymbol,
                          "ECMregulators"       = my.5.ECM$GeneSymbol,
                          "Proteoglycans"       = my.6.ECM$GeneSymbol,
                          "Matrisomeassociated" = my.7.ECM$GeneSymbol,
                          "Secretedfactors"     = my.8.ECM$GeneSymbol)

######################## A. Load DESeq Files ########################
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




################################################################################################
######################## 2. Run GSEA analysis
setwd("C:/Users/livis/Documents/Benayoun_Lab/ecm_11-5")
for(i in 1:length(my_deseq_list)){
  my.ECM.sex = my_deseq_list[[i]]
  my.outprefix <- paste(Sys.Date(), names(my_deseq_list[i]) , sep='_')
  
  my.ECM.sex$GeneName  <- my.ECM.sex$row.names
  
  ECM.sex.geneList         <- my.ECM.sex$stat
  names(ECM.sex.geneList)  <- my.ECM.sex$GeneName
  ECM.sex.geneList         <- sort(ECM.sex.geneList , decreasing = TRUE)
  
  
  set.seed(123456789)
  
  # run phenotest GSEA
  gsea.data <- gsea( x         =  ECM.sex.geneList ,
                     gsets     =  my.ECM.curated.gs,
                     logScale  =  FALSE             ,
                     B         =  8000              ,
                     minGenes  =  5                 ,
                     maxGenes  =  5000              ,
                     center = TRUE)
  my.summary <- data.frame(gsea.data$significance$summary)
  
  
  # write results to file
  my.outfile <- paste(my.outprefix, "sex_ECM_Gene_Lists_GSEA_Analysis_table.txt", sep = "_")
  write.table(my.summary, file = my.outfile, quote = F, sep = "\t")
  
  save(gsea.data , file = paste(my.outprefix,"sex_ECM_Gene_Lists_GSEA.RData", sep = "_") )
}

################################################################################################
######################## 3. Create heatmap

#read in ECM GSEA
BMDM_Cohort1_ <- read.csv("2024-12-05_BMDM_Cohort1_sex_ECM_Gene_Lists_GSEA_Analysis_table.txt"  , header = T, sep = "\t", row.names = 1)
BMDM_Cohort2_ <- read.csv("2024-12-05_BMDM_Cohort2_sex_ECM_Gene_Lists_GSEA_Analysis_table.txt"  , header = T, sep = "\t", row.names = 1)
Peritoneal_Cohort1_ <- read.csv("2024-12-05_Peritoneal_Cohort1_sex_ECM_Gene_Lists_GSEA_Analysis_table.txt", header = T, sep = "\t", row.names = 1)
Peritoneal_Cohort2_ <- read.csv("2024-12-05_Peritoneal_Cohort2_sex_ECM_Gene_Lists_GSEA_Analysis_table.txt"  , header = T, sep = "\t", row.names = 1)
GSE41879_Peritoneal_ <- read.csv("2024-12-05_GSE41879_Peritoneal_sex_ECM_Gene_Lists_GSEA_Analysis_table.txt"  , header = T, sep = "\t", row.names = 1)
GSE99622_Microglia_whole_brain_P14_ <- read.csv("2024-12-05_GSE99622_Microglia_whole_brain_P14_sex_ECM_Gene_Lists_GSEA_Analysis_table.txt"  , header = T, sep = "\t", row.names = 1)
GSE109099_Alveolar_ <- read.csv("2024-12-05_GSE109099_Alveolar_sex_ECM_Gene_Lists_GSEA_Analysis_table.txt", header = T, sep = "\t", row.names = 1)
GSE124829_Microglia_whole_brain_ImmGen_8w_ <- read.csv("2024-12-05_GSE124829_Microglia_whole_brain_ImmGen_sex_ECM_Gene_Lists_GSEA_Analysis_table.txt"  , header = T, sep = "\t", row.names = 1)
GSE124829_Peritoneal_ImmGen_8w_ <- read.csv("2024-12-05_GSE124829_Peritoneal_ImmGen_sex_ECM_Gene_Lists_GSEA_Analysis_table.txt"  , header = T, sep = "\t", row.names = 1)
GSE153299_OCP_MCSF_RANKL_ <- read.csv("2024-12-05_GSE153299_OCP_MCSF_RANKL_sex_ECM_Gene_Lists_GSEA_Analysis_table.txt"  , header = T, sep = "\t", row.names = 1)
GSE153299_OCP_MCSF_ <- read.csv("2024-12-05_GSE153299_OCP_MCSF_sex_ECM_Gene_Lists_GSEA_Analysis_table.txt"  , header = T, sep = "\t", row.names = 1)
GSE174207_Alveolar_CC_ <- read.csv("2024-12-05_GSE174207_Alveolar_CC_sex_ECM_Gene_Lists_GSEA_Analysis_table.txt"  , header = T, sep = "\t", row.names = 1)
GSE149014_Peritoneal_ <- read.csv("2024-12-05_GSE149014_Peritoneal_sex_ECM_Gene_Lists_GSEA_Analysis_table.txt"  , header = T, sep = "\t", row.names = 1)
PRJNA383777_Microglia_whole_brain_ <- read.csv("2024-12-05_PRJNA383777_Microglia_whole_brain_sex_ECM_Gene_Lists_GSEA_Analysis_table.txt"  , header = T, sep = "\t", row.names = 1)
PRJNA408225_Microglia_frontal_lobe_ <- read.csv("2024-12-05_PRJNA408225_Microglia_frontal_lobe_sex_ECM_Gene_Lists_GSEA_Analysis_table.txt"  , header = T, sep = "\t", row.names = 1)
GSE156799_Alveolar_ <- read.csv("2024-12-05_GSE156799_Alveolar_sex_ECM_Gene_Lists_GSEA_Analysis_table.txt"  , header = T, sep = "\t", row.names = 1)
GSE109099_Exudate_ <- read.csv("2024-12-05_GSE109099_Exudate_sex_ECM_Gene_Lists_GSEA_Analysis_table.txt"  , header = T, sep = "\t", row.names = 1)
PRJNA408225_Microglia_hippocampus_ <- read.csv("2024-12-05_PRJNA408225_Microglia_Hippocampus_sex_ECM_Gene_Lists_GSEA_Analysis_table.txt"  , header = T, sep = "\t", row.names = 1)

#put into list
data_list <- list("Peritoneal_Cohort1" = Peritoneal_Cohort1_,
                  "Peritoneal_Cohort2" = Peritoneal_Cohort2_,
                  "GSE41879_Peritoneal" = GSE41879_Peritoneal_,
                  "GSE124829_Peritoneal_ImmGen" = GSE124829_Peritoneal_ImmGen_8w_,
                  "GSE149014_Peritoneal" = GSE149014_Peritoneal_,
                  "GSE99622_Microglia_whole_brain_P14" = GSE99622_Microglia_whole_brain_P14_,
                  "GSE124829_Microglia_whole_brain_ImmGen" = GSE124829_Microglia_whole_brain_ImmGen_8w_,
                  "PRJNA383777_Microglia_whole_brain" = PRJNA383777_Microglia_whole_brain_,
                  "PRJNA408225_Microglia_frontal_lobe" = PRJNA408225_Microglia_frontal_lobe_,
                  "PRJNA408225_Microglia_hippocampus" = PRJNA408225_Microglia_hippocampus_,
                  "GSE109099_Alveolar" = GSE109099_Alveolar_,
                  "GSE174207_Alveolar_CC" = GSE174207_Alveolar_CC_,
                  "GSE156799_Alveolar" = GSE156799_Alveolar_,
                  "GSE153299_OCP_MCSF_RANKL" = GSE153299_OCP_MCSF_RANKL_,
                  "GSE153299_OCP_MCSF" = GSE153299_OCP_MCSF_,
                  "BMDM_Cohort1" = BMDM_Cohort1_,
                  "BMDM_Cohort2" = BMDM_Cohort2_,
                  "GSE109099_Exudate" = GSE109099_Exudate_)


# Assuming `data_list` is your list of datasets
# Combine the "nes" values into a matrix
nes_matrix <- do.call(cbind, lapply(data_list, function(x) x$nes))
rownames(nes_matrix) <- rownames(data_list[[1]])  # Use row names from the first dataset

# Combine the "fdr" values into a matrix for annotation
fdr_matrix <- do.call(cbind, lapply(data_list, function(x) x$fdr))
rownames(fdr_matrix) <- rownames(data_list[[1]])

# Create an annotation matrix with asterisks where fdr < 0.05
annotation_matrix <- ifelse(fdr_matrix < 0.05, "*", "")

# Define the color scale
heatmap_colors <- colorRampPalette(c("deepskyblue", "white", "deeppink"))(100)

# Create the heatmap
my.heatmap.title <- paste("ECM_GSEA_Results_Heatmap.pdf")
pdf(my.heatmap.title)
pheatmap(
  nes_matrix,
  color = heatmap_colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = annotation_matrix,  # Overlay the asterisks
  number_color = "black",  # Color of the asterisks
  fontsize_number = 10  # Adjust the size of the asterisks
)
dev.off()

#######################
sink(file = paste(Sys.Date(),"SMAC_GSEA_ECM_session_Info.txt", sep =""))
sessionInfo()
sink()
