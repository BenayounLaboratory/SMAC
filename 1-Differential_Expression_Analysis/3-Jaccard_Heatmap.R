# 2023-8-30
# Use Jaccard index to create heatmaps of up and down regulated genes


################################################################################################
######################## 1. Define Jaccard function

# logFC.data = my.logFC.data
# FDR.data <- my.FDR.data
# fdr.thres = 0.05

calc_jaccard <- function(logFC.data, FDR.data, fdr.thres = 0.05){
  
  my.up  <- logFC.data > 0
  my.sig <- FDR.data   < fdr.thres
  
  up.genes   <- data.frame(matrix(0,length(my.exp.genes),length(deseq.res.list.genes)))
  colnames(up.genes) <- names(deseq.res.list.genes)
  rownames(up.genes) <- my.exp.genes
  
  down.genes <- data.frame(matrix(0,length(my.exp.genes),length(deseq.res.list.genes)))
  colnames(down.genes) <- names(deseq.res.list.genes)
  rownames(down.genes) <- my.exp.genes
  
  for (i in 1:length(deseq.res.list.genes)) {
    for (j in 1:length(my.exp.genes)) {
      
      up.genes[j,i]   <- bitAnd(my.up[j,i], my.sig[j,i])
      down.genes[j,i] <- bitAnd(!my.up[j,i], my.sig[j,i])
      
    }
  }
  up.genes[is.na(up.genes)]     <- 0
  down.genes[is.na(down.genes)] <- 0
  
  my.jaccard.up <- matrix(0,ncol(up.genes),ncol(up.genes))
  rownames(my.jaccard.up) <- colnames(up.genes)
  colnames(my.jaccard.up) <- colnames(up.genes)
  
  for (i in 1:ncol(up.genes) ) {
    for (j in 1:ncol(up.genes) ) {
      
      my.inter <- sum(bitAnd(up.genes[,i],up.genes[,j]))
      my.union <- sum(bitOr(up.genes[,i],up.genes[,j]))
      
      my.jaccard.up[i,j] <- my.inter/my.union
    }
  }
  
  
  my.jaccard.dwn <- matrix(0,ncol(down.genes),ncol(down.genes))
  rownames(my.jaccard.dwn) <- colnames(down.genes)
  colnames(my.jaccard.dwn) <- colnames(down.genes)
  
  for (i in 1:ncol(down.genes) ) {
    for (j in 1:ncol(down.genes) ) {
      
      my.inter <- sum(bitAnd(down.genes[,i],down.genes[,j]))
      my.union <- sum(bitOr(down.genes[,i],down.genes[,j]))
      
      my.jaccard.dwn[i,j] <- my.inter/my.union
    }
  }
  
  return(list(my.jaccard.up, my.jaccard.dwn))
}

################################################################################################
######################## 2. Create files of top 1500 upregulated and top 1500 downregulated genes
######################## from deseq files via log fold change

library(stringr)
setwd("C:/Users/livis/Documents/Benayoun_Lab/11-5_rerun/deseq")

# Generate list of files from working directory of DESeq files(thus make sure everything in the directory is a deseq file)
files <- list.files(pattern = "\\.txt$")

for( i in 1:length(files)){
  #setwd("C:/Users/livis/Documents/Benayoun_Lab/11-5_rerun/deseq") 
  dataset <- read.csv(files[i],  header = T, sep = "\t")
  dataset_name <- str_extract(files[i], "(?<=\\d{4}-\\d{2}-\\d{2}_).*?(?=_SEX)") 
  
  # Order by log2
  dataset <- dataset[order(dataset$log2FoldChange),]
  
  # Grab top and bottom most changed
  top <- head(dataset,1500)
  bottom <- tail(dataset,1500)
  t3000 <- rbind(top,bottom)
  
  #setwd("C:/Users/livis/Documents/Benayoun_Lab/11-5_rerun/jaccard")
  filename = paste(dataset_name, "50_top3000.txt", sep = '_')
  write.table(t3000, file = filename, quote = F, sep = "\t")
}

################################################################################################
######################## 3. Create the heatmap plot

options(stringsAsFactors = F)
# Load packages
require(DOSE)             # DOSE_3.12.0     
library(ggplot2)          # ggplot2_3.3.5    
library(scales)           # scales_1.1.1 

library("ComplexHeatmap") # 
library("bitops")         # 
library(circlize)         #
library(Polychrome)

theme_set(theme_bw())   


# Load in top 3000 up/downregulated files
setwd("C:/Users/livis/Documents/Benayoun_Lab/11-5_rerun/jaccard")
jaccard_files <- list.files(pattern = "\\.txt$")

my.data.list <- list()
# Loop through each file in the "jaccard_50" directory
for (file in jaccard_files) {
  dataset_name <- str_extract(file, ".*(?=_50_top3000\\.txt)")
  dataset <- read.csv(file, header = TRUE, sep = "\t")
  my.data.list[[dataset_name]] <- dataset
}

# Reset order
new_order <- c(
  "GSE99622_Microglia_whole_brain_P14",
  "GSE124829_Microglia_whole_brain_ImmGen_8w",
  "PRJNA383777_Microglia_whole_brain",
  "PRJNA408225_Microglia_frontal_lobe",
  "PRJNA408225_Microglia_hippocampus",
  "Benayoun_Peritoneal_Cohort1",
  "Benayoun_Peritoneal_Cohort2",
  "GSE41879_Peritoneal",
  "GSE149014_Peritoneal",
  "GSE124829_Peritoneal_ImmGen_6w",
  "Benayoun_BMDM_Cohort1",
  "Benayoun_BMDM_Cohort2",
  "GSE153299_OCP_MCSF",
  "GSE153299_OCP_MCSF_RANKL",
  "GSE109099_Exudate",
  "GSE109099_Alveolar",
  "GSE156799_Alveolar",
  "GSE174207_Alveolar_CC"
)

my.data.list <- my.data.list[new_order]

# Create background
my.exp.genes <- list()
for (i in 1:length(my.data.list)) {
  my.exp.genes <- unique(c(my.exp.genes,row.names(my.data.list[[i]])))
}
length(my.exp.genes)
#[1] 15958

deseq.res.list.genes <- my.data.list

# Create matrices
my.logFC.data <- data.frame(matrix(NA,length(my.exp.genes),length(deseq.res.list.genes)))
colnames(my.logFC.data) <- names(deseq.res.list.genes)
rownames(my.logFC.data) <- my.exp.genes

my.FDR.data <- data.frame(matrix(NA,length(my.exp.genes),length(deseq.res.list.genes)))
colnames(my.FDR.data) <- names(deseq.res.list.genes)
rownames(my.FDR.data) <- my.exp.genes

# Populate matrices
for (i in 1:length(deseq.res.list.genes)) {
  for (j in 1:length(my.exp.genes)) {
    my.idx <- row.names(deseq.res.list.genes[[i]]) %in% my.exp.genes[j]
    
    if (sum(my.idx) > 0) {
      
      my.logFC.data[j,i] <- deseq.res.list.genes[[i]]$log2FoldChange[my.idx]
      my.FDR.data[j,i]   <- deseq.res.list.genes[[i]]$padj[my.idx]
      
    }
  }
}

# Get spearman rank correlation
my.cors <-cor(my.logFC.data, method = 'spearman', use = "complete.obs")

# Calculate Jaccard index
jacc.5 <- calc_jaccard(my.logFC.data, my.FDR.data, 0.05)

# Export to pdf
pdf(paste0(Sys.Date(),"_Jaccard_Index_FDR5.pdf"))

Heatmap(jacc.5[[1]], 
        col = colorRamp2(c(0,1), c("white","deeppink"), transparency = 0, space = "LAB"),
        border = T, rect_gp = gpar(col = "grey", lwd = 0.5)  ,
        column_title = "Jaccard Female-biased")

Heatmap(jacc.5[[2]], 
        col = colorRamp2(c(0,1), c("white","deepskyblue"), transparency = 0, space = "LAB"),
        border = T, rect_gp = gpar(col = "grey", lwd = 0.5) ,
        column_title = "Jaccard Male-biased")
dev.off()

#######################
sink(file = paste(Sys.Date(),"SMAC_Jaccard_RNAseq_analysis_session_Info.txt", sep =""))
sessionInfo()
sink()

