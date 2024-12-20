# 2023-9-18
# Use DESeq files to create a strip plot of differentially expressed genes


################################################################################################
######################## 1. Load in necessary data and packages
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
my_deseq_list <- list("Peritoneal_Cohort1" = Peritoneal_Cohort1_,
                      "Peritoneal_Cohort2" = Peritoneal_Cohort2_,
                      "GSE41879_Peritoneal" = GSE41879_Peritoneal_,
                      "GSE124829_Peritoneal_ImmGen" = GSE124829_Peritoneal_ImmGen_8w_,
                      "GSE149014_Peritoneal" = GSE149014_Peritoneal_,
                      "GSE99622_Microglia_whole_brain_P14" = GSE99622_Microglia_whole_brain_P14_,
                      "GSE124829_Microglia_whole_brain_ImmGen" = GSE124829_Microglia_whole_brain_ImmGen_8w_,
                      "PRJNA383777_Microglia_whole_brain" = PRJNA383777_Microglia_whole_brain_,
                      "PRJNA408225_Microglia_frontal_lobe" = PRJNA408225_Microglia_frontal_lobe_,
                      "PRJNA408225_Microglia_Hippocampus" = PRJNA408225_Microglia_Hippocampus_,
                      "GSE109099_Alveolar" = GSE109099_Alveolar_,
                      "GSE174207_Alveolar_CC" = GSE174207_Alveolar_CC_,
                      "GSE156799_Alveolar" = GSE156799_Alveolar_,
                      "GSE153299_OCP_MCSF_RANKL" = GSE153299_OCP_MCSF_RANKL_,
                      "GSE153299_OCP_MCSF" = GSE153299_OCP_MCSF_,
                      "BMDM_Cohort1" = BMDM_Cohort1_,
                      "BMDM_Cohort2" = BMDM_Cohort2_,
                      "GSE109099_Exudate" = GSE109099_Exudate_)

setwd("C:/Users/livis/Documents/Benayoun_Lab/11-5_rerun")

################################################################################################
######################## 2. Create Jitter plot
# code adapted from https://github.com/brunetlab/Leeman_et_al_2017/blob/master/kallisto_deseq2/Fig4A_stripplot_cell_type_colors.R

color_palette <- list(
  
  #peritoneal:
  "#495ac6",
  "#68a5ff",
  "#205883",
  "#76a6d6",
  "#016fad",
  
  #microglia:
  "#db3e68",
  "#892b36",
  "#da7579",
  "#d83c34",
  "#ff2400", 
  
  #alveolar:
  "#81d54b",
  "#8bd68e",
  "#4b7f39",
  
  #OCP:
  "#ff5722",
  "#fa9600",
  
  #BMDM:
  "#9c27b0",
  "#4c1273",
  
  #exudate:
  "#FFDA00")

# Get correspondence cell type/color
ct.unique <- data.frame("label" = names(my_deseq_list), "cols" = unlist(color_palette))
ct.unique
# Order by pvalue:
sex.results <- lapply(my_deseq_list,function(x) {x[order(x$padj),]})
n        <- sapply(sex.results, nrow)
names(n) <- names(sex.results)

# Assign Colors
cols <- list()
xlab <- character(length = length(sex.results))
for(i in seq(along = sex.results)){
  cols[[i]] <- rep(rgb(153, 153, 153, maxColorValue = 255, alpha = 70), n[i]) # grey60
  ind.sig.i <- sex.results[[i]]$padj < 0.05
  cols[[i]][ind.sig.i] <- ct.unique[i, "cols"]
  xlab[i] <- paste(names(sex.results)[i], "\n(", sum(ind.sig.i), " sig.)", sep = "")
}
names(cols) <- names(sex.results)

# Save to PDF
pdf(paste0(Sys.Date(),"_stripplot_DESeq2_with_cell_type_colors.pdf"), width = 9, height = 6)

# Create plot
par(mar = c(3.1, 4.1, 1, 1))
par(oma = c(6, 2, 1, 1))
plot(x = 1,
     y = 1,
     type = "n",
     xlim = c(0.5, 18.5),
     ylim = c(-15, 20),
     axes = FALSE,
     xlab = "",
     ylab = "Log2 fold change (F / M)"
)
x_range <- c(-0.5, 19.5)  # Get the x range of the plot
y_range <- c(-16.5, 21.5)  # Get the y range of the plot

# Create a polygon for the area above zero
polygon(
  c(x_range[1], x_range[2], x_range[2], x_range[1]), 
  c(0, 0, y_range[2], y_range[2]), 
  col = "#FDF3F8", border = NA
)

# Create a polygon for the area below zero
polygon(
  c(x_range[1], x_range[2], x_range[2], x_range[1]), 
  c(0, 0, y_range[1], y_range[1]), 
  col = "#F5FCFE", border = NA
)


abline(h = 0)
abline(h = seq(-20, 20, by = 5)[-5],
       lty = "dotted",
       col = "grey")
for(i in 1:length(sex.results)){
  set.seed(1234)
  points(x = jitter(rep(i, nrow(sex.results[[i]])), amount = 0.2),
         y = rev(sex.results[[i]]$log2FoldChange),
         pch = 16,
         col = rev(cols[[i]]),
         bg = rev(cols[[i]]),
         cex = 0.4)
}
axis(1,
     at = 1:18,
     tick = FALSE,
     las = 2,
     lwd = 0,
     labels = xlab,
     cex.axis = 0.7)
axis(2,
     las = 1,
     at = seq(-15, 20, 5))
box()
dev.off()
###############################################################################################

sink(file = paste(Sys.Date(),"SMAC_stripplot_RNAseq_analysis_session_Info.txt", sep =""))
sessionInfo()
sink()

