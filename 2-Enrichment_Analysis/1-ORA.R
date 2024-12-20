# 2024-9-30
# Find shared differentially expressed genes and perform over representation analysis

################################################################################################
######################## 1. Load in necessary data and packages

library(UpSetR)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)
library(stringr)
require(DOSE)

#read in FDR files
setwd("C:/Users/livis/Documents/Benayoun_Lab/11-5_rerun/fdr")
BMDM_Cohort1 <- read.csv("2024-11-05_Benayoun_BMDM_Cohort1_SEX_DIM_FDR5_genes_statistics.txt"  , header = T, sep = "\t", row.names = NULL)
BMDM_Cohort2 <- read.csv("2024-11-05_Benayoun_BMDM_Cohort2_SEX_DIM_FDR5_genes_statistics.txt"  , header = T, sep = "\t", row.names = NULL)
Peritoneal_Cohort1 <- read.csv("2024-11-05_Benayoun_Peritoneal_Cohort1_SEX_DIM_FDR5_genes_statistics.txt", header = T, sep = "\t", row.names = NULL)
Peritoneal_Cohort2 <- read.csv("2024-11-05_Benayoun_Peritoneal_Cohort2_SEX_DIM_FDR5_genes_statistics.txt"  , header = T, sep = "\t", row.names = NULL)
GSE41879_Peritoneal <- read.csv("2024-11-05_GSE41879_Peritoneal_SEX_DIM_FDR5_genes_statistics.txt"  , header = T, sep = "\t", row.names = NULL)
GSE99622_Microglia_whole_brain_P14 <- read.csv("2024-11-05_GSE99622_Microglia_whole_brain_P14_SEX_DIM_FDR5_genes_statistics.txt"  , header = T, sep = "\t", row.names = NULL)
GSE109099_Alveolar <- read.csv("2024-11-05_GSE109099_Alveolar_SEX_DIM_FDR5_genes_statistics.txt", header = T, sep = "\t", row.names = NULL)
GSE124829_Microglia_whole_brain_ImmGen_8w <- read.csv("2024-11-05_GSE124829_Microglia_whole_brain_ImmGen_8w_SEX_DIM_FDR5_genes_statistics.txt"  , header = T, sep = "\t", row.names = NULL)
GSE124829_Peritoneal_ImmGen_8w <- read.csv("2024-11-05_GSE124829_Peritoneal_ImmGen_6w_SEX_DIM_FDR5_genes_statistics.txt"  , header = T, sep = "\t", row.names = NULL)
GSE153299_OCP_MCSF_RANKL <- read.csv("2024-11-05_GSE153299_OCP_MCSF_RANKL_SEX_DIM_FDR5_genes_statistics.txt"  , header = T, sep = "\t", row.names = NULL)
GSE153299_OCP_MCSF <- read.csv("2024-11-05_GSE153299_OCP_MCSF_SEX_DIM_FDR5_genes_statistics.txt"  , header = T, sep = "\t", row.names = NULL)
GSE174207_Alveolar_CC <- read.csv("2024-11-05_GSE174207_Alveolar_CC_SEX_DIM_FDR5_genes_statistics.txt"  , header = T, sep = "\t", row.names = NULL)
GSE149014_Peritoneal <- read.csv("2024-11-05_GSE149014_Peritoneal_SEX_DIM_FDR5_genes_statistics.txt"  , header = T, sep = "\t", row.names = NULL)
PRJNA383777_Microglia_whole_brain <- read.csv("2024-11-05_PRJNA383777_Microglia_whole_brain_SEX_DIM_FDR5_genes_statistics.txt"  , header = T, sep = "\t", row.names = NULL)
PRJNA408225_Microglia_frontal_lobe <- read.csv("2024-11-05_PRJNA408225_Microglia_frontal_lobe_SEX_DIM_FDR5_genes_statistics.txt"  , header = T, sep = "\t", row.names = NULL)
GSE156799_Alveolar <- read.csv("2024-11-05_GSE156799_Alveolar_SEX_DIM_FDR5_genes_statistics.txt"  , header = T, sep = "\t", row.names = NULL)
GSE109099_Exudate <- read.csv("2024-11-05_GSE109099_Exudate_SEX_DIM_FDR5_genes_statistics.txt"  , header = T, sep = "\t", row.names = NULL)
PRJNA408225_Microglia_Hippocampus <- read.csv("2024-11-05_PRJNA408225_microglia_Hippocampus_SEX_DIM_FDR5_genes_statistics.txt"  , header = T, sep = "\t", row.names = NULL)

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

setwd("C:/Users/livis/Documents/Benayoun_Lab/11-5_rerun/ora")

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


my_DEGS <- list("GSE99622_Microglia_whole_brain_P14" = GSE99622_Microglia_whole_brain_P14$row.names,
                "GSE124829_Microglia_whole_brain_ImmGen_8w" = GSE124829_Microglia_whole_brain_ImmGen_8w$row.names, 
                "PRJNA383777_Microglia_whole_brain" = PRJNA383777_Microglia_whole_brain$row.names ,
                "PRJNA408225_Microglia_frontal_lobe" = PRJNA408225_Microglia_frontal_lobe$row.names,
                "PRJNA408225_Microglia_Hippocampus" = PRJNA408225_Microglia_Hippocampus$row.names, 
                "Peritoneal_Cohort1" = Peritoneal_Cohort1$row.names,
                "Peritoneal_Cohort2" = Peritoneal_Cohort2$row.names,
                "GSE41879_Peritoneal" = GSE41879_Peritoneal$row.names,
                "GSE149014_Peritoneal" = GSE149014_Peritoneal$row.names,
                "GSE124829_Peritoneal" = GSE124829_Peritoneal_ImmGen_8w$row.names,
                "BMDM_Cohort_1" = BMDM_Cohort1$row.names, 
                "BMDM_Cohort_2" = BMDM_Cohort2$row.names,
                "GSE153299_OCP_MCSF" = GSE153299_OCP_MCSF$row.names,
                "GSE153299_OCP_MCSF_RANKL" = GSE153299_OCP_MCSF_RANKL$row.names,
                "GSE109099_Exudate" = GSE109099_Exudate$row.names, 
                "GSE109099_Alveolar" = GSE109099_Alveolar$row.names, 
                "GSE156799_Alveolar" = GSE156799_Alveolar$row.names,
                "GSE174207_Alveolar_CC" = GSE174207_Alveolar_CC$row.names
)

all_DEGs <- unique(unlist(my_DEGS))
length(all_DEGs)
#[1] 11625

all_DEGs_with_dups <- unlist(my_DEGS)
length(all_DEGs_with_dups)
#[1] 23000

# Create matrix to find genes most commonly expressed across datasets
my.all.mat                 <- data.frame(matrix(0,length(all_DEGs),length(my_DEGS)))
colnames(my.all.mat)       <- names(my_DEGS)
rownames(my.all.mat)       <- all_DEGs

for (h in 1:length(my_DEGS)) {
  my.all.mat[all_DEGs %in% my_DEGS[[h]],h] <- 1
}

dim(my.all.mat)
#[1] 11625    18


# Find genes expressed across all datasets
sum(rowSums(my.all.mat) == 18)
#3
my.all.expressed <- all_DEGs[rowSums(my.all.mat) == 18]
print(my.all.expressed, quote = F)

#Xist    Eif2s3y Ddx3y  

################################################################################################
######################## 2. Now, split into body parts and repeat to find genes shared across body type macrophages
# then perform functional analysis to find affected pathways by body type

##########################################################################################
# Peritoneal

peritoneal_DEGs <- list( "Peritoneal_Cohort1" = Peritoneal_Cohort1$row.names,
                         "Peritoneal_Cohort2" = Peritoneal_Cohort2$row.names,
                         "GSE41879_Peritoneal" = GSE41879_Peritoneal$row.names,
                         "GSE149014_Peritoneal" = GSE149014_Peritoneal$row.names,
                         "GSE124829_Peritoneal" = GSE124829_Peritoneal_ImmGen_8w$row.names
)

all_peritoneal_genes <- unique(unlist(peritoneal_DEGs))
length(all_peritoneal_genes)
#[1] 3487

# Create matrix
my.peritoneal.mat                 <- data.frame(matrix(0,length(all_peritoneal_genes),length(peritoneal_DEGs)))
colnames(my.peritoneal.mat)       <- names(peritoneal_DEGs)
rownames(my.peritoneal.mat)       <- all_peritoneal_genes
# Populate matrix
for (h in 1:length(peritoneal_DEGs)) {
  my.peritoneal.mat[all_peritoneal_genes %in% peritoneal_DEGs[[h]],h] <- 1
}

dim(my.peritoneal.mat)
#[1] 3487    5


# Find genes shared across all datasets
sum(rowSums(my.peritoneal.mat) == 5)
#[1] 20
my.peritoneal.expressed <- all_peritoneal_genes[rowSums(my.peritoneal.mat) == 5]
print(my.peritoneal.expressed, quote = F)

# Export all shared DEGs to txt
write.table(my.peritoneal.expressed, file = "peritoneal_shared_DEGS_50.txt", sep = "\t" , row.names = F, quote=F)

##############################
# Divide up and down regulated for over representation analysis

peritoneal_dfs <- list( "Peritoneal_Cohort1" = Peritoneal_Cohort1,
                        "Peritoneal_Cohort2" = Peritoneal_Cohort2,
                        "GSE41879_Peritoneal" = GSE41879_Peritoneal,
                        "GSE149014_Peritoneal" = GSE149014_Peritoneal,
                        "GSE124829_Peritoneal" = GSE124829_Peritoneal_ImmGen_8w
)

peritoneal_all_expressed <- list("Peritoneal_Cohort1" = Peritoneal_Cohort1_$row.names,
                                 "Peritoneal_Cohort2" = Peritoneal_Cohort2_$row.names,
                                 "GSE41879_Peritoneal" = GSE41879_Peritoneal_$row.names,
                                 "GSE149014_Peritoneal" = GSE149014_Peritoneal_$row.names,
                                 "GSE124829_Peritoneal" = GSE124829_Microglia_whole_brain_ImmGen_8w_$row.names)

all_peritoneal_expressed <- unique(unlist(peritoneal_all_expressed))

# Create matrix
my.fc.peritoneal.mat                 <- data.frame(matrix(0,length(all_peritoneal_genes),length(peritoneal_DEGs)))
colnames(my.fc.peritoneal.mat)       <- names(peritoneal_DEGs)
rownames(my.fc.peritoneal.mat)       <- all_peritoneal_genes

# Populate matrix
for (h in 1:length(peritoneal_DEGs)) {
  my.fc.peritoneal.mat[all_peritoneal_genes %in% peritoneal_DEGs[[h]],h] <- 1
}

for(h in  1:length(peritoneal_DEGs)) {
  for(g in 1:length(all_peritoneal_genes)){
    if(my.fc.peritoneal.mat[g,h] == 1){
      gene <- row.names(my.fc.peritoneal.mat[g,])
      df <- peritoneal_dfs[[h]]
      subdf <- subset(df, df$row.names == gene)
      my.fc.peritoneal.mat[g,h] <- subdf$log2FoldChange
    }
  }
}

# Replace 0s with NAs
my.fc.peritoneal.mat[my.fc.peritoneal.mat == 0] <- NA

# Add new column with count of NAs
my.fc.peritoneal.mat$Total_NAs <- rowSums(is.na(my.fc.peritoneal.mat))

# Add new columns with counts of positive and negative log2fc values
my.fc.peritoneal.mat$up <- rowSums(my.fc.peritoneal.mat[,1:5] > 0, na.rm = TRUE)
my.fc.peritoneal.mat$down <- rowSums(my.fc.peritoneal.mat[,1:5] < 0, na.rm = TRUE)

# If the gene is upregulated in at least 3/5 of them and present in at least 4/5
my.up.peritoneal <- all_peritoneal_genes[my.fc.peritoneal.mat$up >= 3 & my.fc.peritoneal.mat$Total_NAs < 2]
my.up.peritoneal
length(my.up.peritoneal)
#[1] 78

my.down.peritoneal <- all_peritoneal_genes[my.fc.peritoneal.mat$down >= 3 & my.fc.peritoneal.mat$Total_NAs < 2]
my.down.peritoneal
length(my.down.peritoneal)
#[1] 42

# Perform ORA on upregulated genes
entrezID.genes.shared <- bitr(my.up.peritoneal , fromType="SYMBOL", toType="ENTREZID", 
                              OrgDb="org.Mm.eg.db")

# In bitr(my.up.peritoneal, fromType = "SYMBOL", toType = "ENTREZID",  :
# 1.28% of input gene IDs are fail to map...
# OK to proceed when low amount of input genes fail to map.

entrezID.genes.universe <- bitr(all_peritoneal_expressed ,fromType="SYMBOL", toType="ENTREZID", 
                                OrgDb="org.Mm.eg.db")
#In bitr(all_peritoneal_expressed, fromType = "SYMBOL", toType = "ENTREZID",  :
#4.99% of input gene IDs are fail to map...

# Run over-representation analysis
ego.shared <- enrichGO(gene   = entrezID.genes.shared$ENTREZID,
                       universe  = entrezID.genes.universe$ENTREZID,
                       OrgDb         = org.Mm.eg.db,
                       keyType       = 'ENTREZID',
                       ont           = "ALL",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE,
                       minGSSize     = 10  ,
                       maxGSSize     = 5000)

peri_up_ora <- as.data.frame(ego.shared)
#91 rows

# Export results to file
write.table(peri_up_ora, file = paste(Sys.Date(),"Peritoneal_upregulated_ORA_loosereqs_50.txt", sep ="_"), sep = "\t" , row.names = F, quote=F)

#sort and reduce for plot
df <- peri_up_ora
df <- df[order(df$qvalue),]
df <- df[1:10,]

# Make plot
#Make gene ratio into a decimal
df$GeneRatioNumeric <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), `[`, 1)) / 
  as.numeric(sapply(strsplit(df$GeneRatio, "/"), `[`, 2))

# Create a base plot
plot <- ggplot(df, aes(x = 1, y = Description, size = GeneRatioNumeric)) + 
  geom_point(alpha = 0.6) +
  scale_size_continuous(name = "Gene Ratio")

#Set color scale
plot <- plot + aes(color = -log10(pvalue)) + 
  scale_color_gradient2(low = "white", high = "deeppink", 
                        name = "-log10(p-value)" , limits=c(0,15))


# Final touches
plot <- plot + theme(axis.title.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.x = element_blank(),
                     panel.grid.major.y = element_line(color="#D3D3D3"),
                     panel.background = element_blank(),
                     panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(title = "Peritoneal UpRegulated ORA Bubble Plot")

pdf(paste(Sys.Date(),"peritoneal_upregulated_DEGs_bubbleplot.pdf", 
          sep = "_"))
print(plot)
dev.off()

# Perform ORA on downregulated genes

entrezID.genes.shared <- bitr(my.down.peritoneal , fromType="SYMBOL", toType="ENTREZID", 
                              OrgDb="org.Mm.eg.db")
#'select()' returned 1:1 mapping between keys and columns
#'
entrezID.genes.universe <- bitr(all_peritoneal_expressed ,fromType="SYMBOL", toType="ENTREZID", 
                                OrgDb="org.Mm.eg.db")
#'select()' returned 1:1 mapping between keys and columns
#Warning message:
# In bitr(all_peritoneal_expressed, fromType = "SYMBOL", toType = "ENTREZID",  :
#            4.99% of input gene IDs are fail to map...

# Run over-representation analysis
ego.shared <- enrichGO(gene   = entrezID.genes.shared$ENTREZID,
                       universe  = entrezID.genes.universe$ENTREZID,
                       OrgDb         = org.Mm.eg.db,
                       keyType       = 'ENTREZID',
                       ont           = "ALL",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE,
                       minGSSize     = 10  ,
                       maxGSSize     = 5000)


peri_down_ora <- as.data.frame(ego.shared)
#empty


#################################################################################
# Microglia

microglia_DEGs <- list("GSE99622_Microglia_whole_brain_P14" = GSE99622_Microglia_whole_brain_P14$row.names,
                       "GSE124829_Microglia_whole_brain_ImmGen_8w" = GSE124829_Microglia_whole_brain_ImmGen_8w$row.names, 
                       "PRJNA383777_Microglia_whole_brain" = PRJNA383777_Microglia_whole_brain$row.names ,
                       "PRJNA408225_Microglia_frontal_lobe" = PRJNA408225_Microglia_frontal_lobe$row.names,
                       "PRJNA408225_Microglia_Hippocampus" = PRJNA408225_Microglia_Hippocampus$row.names
)


all_microglia_genes <- unique(unlist(microglia_DEGs))
length(all_microglia_genes)
#[1] 5145

# Create matrix
my.microglia.mat                 <- data.frame(matrix(0,length(all_microglia_genes),length(microglia_DEGs)))
colnames(my.microglia.mat)       <- names(microglia_DEGs)
rownames(my.microglia.mat)       <- all_microglia_genes

# Populate matrix
for (h in 1:length(microglia_DEGs)) {
  my.microglia.mat[all_microglia_genes %in% microglia_DEGs[[h]],h] <- 1
}

dim(my.microglia.mat)
#[1] 2557    4

# Find genes shared across all datasets
sum(rowSums(my.microglia.mat) == 5)
#[1] 5
my.microglia.expressed <- all_microglia_genes[rowSums(my.microglia.mat) == 5]
my.microglia.expressed
print(my.microglia.expressed, quote = F)

# Export all shared DEGs to txt
write.table(my.microglia.expressed, file = "microglia_shared_DEGS_50.txt", sep = "\t" , row.names = F, quote=F)

#####################
# Divide up and downregulation for over representation analysis

microglia_dfs <- list("GSE99622_Microglia_whole_brain_P14" = GSE99622_Microglia_whole_brain_P14,
                      "GSE124829_Microglia_whole_brain_ImmGen_8w" = GSE124829_Microglia_whole_brain_ImmGen_8w, 
                      "PRJNA383777_Microglia_whole_brain" = PRJNA383777_Microglia_whole_brain,
                      "PRJNA408225_Microglia_frontal_lobe" = PRJNA408225_Microglia_frontal_lobe,
                      "PRJNA408225_Microglia_Hippocampus" = PRJNA408225_Microglia_Hippocampus_
)
#load in all genes for ORA background
microglia_all_expressed <- list("GSE99622_Microglia_whole_brain_P14" = GSE99622_Microglia_whole_brain_P14_$row.names,
                      "GSE124829_Microglia_whole_brain_ImmGen_8w" = GSE124829_Microglia_whole_brain_ImmGen_8w_$row.names, 
                      "PRJNA383777_Microglia_whole_brain" = PRJNA383777_Microglia_whole_brain_$row.names,
                      "PRJNA408225_Microglia_frontal_lobe" = PRJNA408225_Microglia_frontal_lobe_$row.names,
                      "PRJNA408225_Microglia_Hippocampus" = PRJNA408225_Microglia_Hippocampus_$row.names)

all_microglia_expressed <- unique(unlist(microglia_all_expressed))


# Create matrix
my.fc.microglia.mat                 <- data.frame(matrix(0,length(all_microglia_genes),length(microglia_DEGs)))
colnames(my.fc.microglia.mat)       <- names(microglia_DEGs)
rownames(my.fc.microglia.mat)       <- all_microglia_genes

# Populate matrix
for (h in 1:length(microglia_DEGs)) {
  my.fc.microglia.mat[all_microglia_genes %in% microglia_DEGs[[h]],h] <- 1
}
for(h in  1:length(microglia_DEGs)) {
  for(g in 1:length(all_microglia_genes)){
    if(my.fc.microglia.mat[g,h] == 1){
      gene <- row.names(my.fc.microglia.mat[g,])
      df <- microglia_dfs[[h]]
      subdf <- subset(df, df$row.names == gene)
      my.fc.microglia.mat[g,h] <- subdf$log2FoldChange
    }
  }
}

# Replace 0s with NAs
my.fc.microglia.mat[my.fc.microglia.mat == 0] <- NA

# Add new column with count of NAs
my.fc.microglia.mat$Total_NAs <- rowSums(is.na(my.fc.microglia.mat))

# Add new columns with counts of positive and negative log2fc values
my.fc.microglia.mat$up <- rowSums(my.fc.microglia.mat[,1:5] > 0, na.rm = TRUE)
my.fc.microglia.mat$down <- rowSums(my.fc.microglia.mat[,1:5] < 0, na.rm = TRUE)

my.up.microglia <- all_microglia_genes[my.fc.microglia.mat$up >= 3 & my.fc.microglia.mat$Total_NAs < 2]
my.up.microglia
length(my.up.microglia)
#[1] 9

my.down.microglia <- all_microglia_genes[my.fc.microglia.mat$down >= 3 & my.fc.microglia.mat$Total_NAs < 2]
my.down.microglia
length(my.down.microglia)
#[1] 4


# Perform ORA on upregulated genes

entrezID.genes.shared <- bitr(my.up.microglia , fromType="SYMBOL", toType="ENTREZID", 
                              OrgDb="org.Mm.eg.db")
#'select()' returned 1:1 mapping between keys and columns
#Warning message:
#  In bitr(my.up.microglia, fromType = "SYMBOL", toType = "ENTREZID",  :
#            11.11% of input gene IDs are fail to map...

entrezID.genes.universe <- bitr(all_microglia_expressed ,fromType="SYMBOL", toType="ENTREZID", 
                                OrgDb="org.Mm.eg.db")
#'select()' returned 1:1 mapping between keys and columns
#Warning message:
#  In bitr(all_microglia_expressed, fromType = "SYMBOL", toType = "ENTREZID",  :
 #           4.91% of input gene IDs are fail to map...

# Run over-representation analysis
ego.shared <- enrichGO(gene   = entrezID.genes.shared$ENTREZID,
                       universe  = entrezID.genes.universe$ENTREZID,
                       OrgDb         = org.Mm.eg.db,
                       keyType       = 'ENTREZID',
                       ont           = "ALL",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE,
                       minGSSize     = 10  ,
                       maxGSSize     = 5000)

micro_up_ora <- as.data.frame(ego.shared)
#empty


# Perform ORA on downregulated genes

entrezID.genes.shared <- bitr(my.down.microglia , fromType="SYMBOL", toType="ENTREZID", 
                              OrgDb="org.Mm.eg.db")
#'select()' returned 1:1 mapping between keys and columns
entrezID.genes.universe <- bitr(all_microglia_expressed ,fromType="SYMBOL", toType="ENTREZID", 
                                OrgDb="org.Mm.eg.db")
#'select()' returned 1:1 mapping between keys and columns
#Warning message:
#  In bitr(all_microglia_expressed, fromType = "SYMBOL", toType = "ENTREZID",  :
 #           4.91% of input gene IDs are fail to map...

# Run over-representation analysis
ego.shared <- enrichGO(gene   = entrezID.genes.shared$ENTREZID,
                       universe  = entrezID.genes.universe$ENTREZID,
                       OrgDb         = org.Mm.eg.db,
                       keyType       = 'ENTREZID',
                       ont           = "ALL",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE,
                       minGSSize     = 10  ,
                       maxGSSize     = 5000)

micro_down_ora <- as.data.frame(ego.shared)
#23 rows

# Export to file
write.table(micro_down_ora, file = paste(Sys.Date(),"microglia_downregulated_ORA_loosereqs_50_universe.txt", sep ="_"), sep = "\t" , row.names = F, quote=F)

#sort and reduce for plot
df <- micro_down_ora
df <- df[order(df$qvalue),]
df <- df[1:10,]

# Make plot
df$GeneRatioNumeric <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), `[`, 1)) / 
  as.numeric(sapply(strsplit(df$GeneRatio, "/"), `[`, 2))

# Create a base plot
plot <- ggplot(df, aes(x = 1, y = Description, size = GeneRatioNumeric)) + 
  geom_point(alpha = 0.6) +
  scale_size_continuous(name = "Gene Ratio")

#Set color scale
plot <- plot + aes(color = -log10(pvalue)) + 
  scale_color_gradient2(low = "white", high = "blue", 
                        name = "-log10(p-value)" , limits=c(0,5))


# Final touches
plot <- plot + theme(axis.title.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.x = element_blank(),
                     panel.grid.major.y = element_line(color="#D3D3D3"),
                     panel.background = element_blank(),
                     panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(title = "Microglia DownRegulated ORA Bubble Plot")+
  scale_y_discrete(labels = function(y) str_wrap(y, width = 10)) #for long descriptions, wrap the text

pdf(paste(Sys.Date(),"microglia_downregulated_DEGs_bubbleplot.pdf", 
          sep = "_"))
print(plot)
dev.off()

###########################################################################################33
#Alveolar

alveolar_DEGs <- list("GSE109099_Alveolar" = GSE109099_Alveolar$row.names, 
                      "GSE156799_Alveolar" = GSE156799_Alveolar$row.names,
                      "GSE174207_Alveolar_CC" = GSE174207_Alveolar_CC$row.names
)

all_alveolar_genes <- unique(unlist(alveolar_DEGs))
length(all_alveolar_genes)
#[1] 7198

# Create matrix
my.alveolar.mat                 <- data.frame(matrix(0,length(all_alveolar_genes),length(alveolar_DEGs)))
colnames(my.alveolar.mat)       <- names(alveolar_DEGs)
rownames(my.alveolar.mat)       <- all_alveolar_genes

# Populate matrix
for (h in 1:length(alveolar_DEGs)) {
  my.alveolar.mat[all_alveolar_genes %in% alveolar_DEGs[[h]],h] <- 1
}

dim(my.alveolar.mat)
#[1] 7198    3


# Find genes shared across all datasets
sum(rowSums(my.alveolar.mat) == 3)
#[1] 4
my.alveolar.expressed <- all_alveolar_genes[rowSums(my.alveolar.mat) == 3]
my.alveolar.expressed
print(my.alveolar.expressed, quote = F)

# Export all shared DEGs to txt
write.table(my.alveolar.expressed, file = "alveolar_shared_DEGS_50.txt", sep = "\t" , row.names = F, quote=F)


##################
# Divide up and downregulation for over representation analysis

alveolar_dfs <- list("GSE109099_Alveolar" = GSE109099_Alveolar, 
                     "GSE156799_Alveolar" = GSE156799_Alveolar,
                     "GSE174207_Alveolar_CC" = GSE174207_Alveolar_CC
)
alveolar_all_expressed <- list("GSE109099_Alveolar" = GSE109099_Alveolar_$row.names, 
                      "GSE156799_Alveolar" = GSE156799_Alveolar_$row.names,
                      "GSE174207_Alveolar_CC" = GSE174207_Alveolar_CC_$row.names
)

all_alveolar_expressed <- unique(unlist(alveolar_all_expressed))

# Create matrix
my.fc.alveolar.mat                 <- data.frame(matrix(0,length(all_alveolar_genes),length(alveolar_DEGs)))
colnames(my.fc.alveolar.mat)       <- names(alveolar_DEGs)
rownames(my.fc.alveolar.mat)       <- all_alveolar_genes

# Populate matrix
for (h in 1:length(alveolar_DEGs)) {
  my.fc.alveolar.mat[all_alveolar_genes %in% alveolar_DEGs[[h]],h] <- 1
}
for(h in  1:length(alveolar_DEGs)) {
  for(g in 1:length(all_alveolar_genes)){
    if(my.fc.alveolar.mat[g,h] == 1){
      gene <- row.names(my.fc.alveolar.mat[g,])
      df <- alveolar_dfs[[h]]
      subdf <- subset(df, df$row.names == gene)
      my.fc.alveolar.mat[g,h] <- subdf$log2FoldChange
    }
  }
}

# Replace 0s with NAs
my.fc.alveolar.mat[my.fc.alveolar.mat == 0] <- NA

# Add new column with count of NAs
my.fc.alveolar.mat$Total_NAs <- rowSums(is.na(my.fc.alveolar.mat))

# Add new columns with counts of positive and negative log2fc values
my.fc.alveolar.mat$up <- rowSums(my.fc.alveolar.mat[,1:3] > 0, na.rm = TRUE)
my.fc.alveolar.mat$down <- rowSums(my.fc.alveolar.mat[,1:3] < 0, na.rm = TRUE)

# If upregulated in at least 2/3
my.up.alveolar <- all_alveolar_genes[my.fc.alveolar.mat$up >= 2]
my.up.alveolar
length(my.up.alveolar)
#[1] 14

# If downregulated in at least 2/3
my.down.alveolar <- all_alveolar_genes[my.fc.alveolar.mat$down >= 2]
my.down.alveolar
length(my.down.alveolar)
#[1] 21

# Perform ORA on upregulated genes

entrezID.genes.shared <- bitr(my.up.alveolar , fromType="SYMBOL", toType="ENTREZID", 
                              OrgDb="org.Mm.eg.db")
#'select()' returned 1:1 mapping between keys and columns
entrezID.genes.universe <- bitr(all_alveolar_expressed ,fromType="SYMBOL", toType="ENTREZID", 
                                OrgDb="org.Mm.eg.db")
#'select()' returned 1:1 mapping between keys and columns
#Warning message:
#  In bitr(all_alveolar_expressed, fromType = "SYMBOL", toType = "ENTREZID",  :
#           5.07% of input gene IDs are fail to map...

# Run over-representation analysis
ego.shared <- enrichGO(gene   = entrezID.genes.shared$ENTREZID,
                       universe  = entrezID.genes.universe$ENTREZID,
                       OrgDb         = org.Mm.eg.db,
                       keyType       = 'ENTREZID',
                       ont           = "ALL",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE,
                       minGSSize     = 10  ,
                       maxGSSize     = 5000)


alveolar_up_ora <- as.data.frame(ego.shared)
#empty

# Perform ORA on downregulated genes

entrezID.genes.shared <- bitr(my.down.alveolar , fromType="SYMBOL", toType="ENTREZID", 
                              OrgDb="org.Mm.eg.db")
#'select()' returned 1:1 mapping between keys and columns
entrezID.genes.universe <- bitr(all_alveolar_expressed ,fromType="SYMBOL", toType="ENTREZID", 
                                OrgDb="org.Mm.eg.db")
#'select()' returned 1:1 mapping between keys and columns
#Warning message:
#  In bitr(all_alveolar_expressed, fromType = "SYMBOL", toType = "ENTREZID",  :
#            5.07% of input gene IDs are fail to map...
# run over-representation analysis
ego.shared <- enrichGO(gene   = entrezID.genes.shared$ENTREZID,
                       universe  = entrezID.genes.universe$ENTREZID,
                       OrgDb         = org.Mm.eg.db,
                       keyType       = 'ENTREZID',
                       ont           = "ALL",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE,
                       minGSSize     = 10  ,
                       maxGSSize     = 5000)


alveolar_down_ora <- as.data.frame(ego.shared)
#2 rows

# Export to file
write.table(alveolar_down_ora, file = paste(Sys.Date(),"alveolar_ORA_downregulated_50_background.txt", sep ="_"), sep = "\t" , row.names = F, quote=F)

#sort for plot
df <- alveolar_down_ora
df <- df[order(df$qvalue),]

# Make plot
df$GeneRatioNumeric <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), `[`, 1)) / 
  as.numeric(sapply(strsplit(df$GeneRatio, "/"), `[`, 2))

# Create a base plot
plot <- ggplot(df, aes(x = 1, y = Description, size = GeneRatioNumeric)) + 
  geom_point(alpha = 0.6) +
  scale_size_continuous(name = "Gene Ratio")

#Set color scale
plot <- plot + aes(color = -log10(pvalue)) + 
  scale_color_gradient2(low = "white", high = "blue", 
                        name = "-log10(p-value)" , limits=c(0,5))


# Final touches
plot <- plot + theme(axis.title.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.x = element_blank(),
                     panel.grid.major.y = element_line(color="#D3D3D3"),
                     panel.background = element_blank(),
                     panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(title = "Alveolar DownRegulated ORA Bubble Plot")+
  scale_y_discrete(labels = function(y) str_wrap(y, width = 10)) #for long descriptions, wrap the text

pdf(paste(Sys.Date(),"alveolar_downregulated_DEGs_bubbleplot.pdf", 
          sep = "_"))
print(plot)
dev.off()

#########################################################################3
#BMDM

bmdm_DEGs <- list( "BMDM_Cohort_1" = BMDM_Cohort1$row.names, 
                   "BMDM_Cohort_2" = BMDM_Cohort2$row.names)

all_bmdm_genes <- unique(unlist(bmdm_DEGs))
length(all_bmdm_genes)
#[1] 2541

# Create matrix
my.bmdm.mat                 <- data.frame(matrix(0,length(all_bmdm_genes),length(bmdm_DEGs)))
colnames(my.bmdm.mat)       <- names(bmdm_DEGs)
rownames(my.bmdm.mat)       <- all_bmdm_genes

# Populate matrix
for (h in 1:length(bmdm_DEGs)) {
  my.bmdm.mat[all_bmdm_genes %in% bmdm_DEGs[[h]],h] <- 1
}

dim(my.bmdm.mat)
#[1] 2541    2


# Find genes shared across all datasets
sum(rowSums(my.bmdm.mat) == 2)
#[1] 28
my.bmdm.expressed <- all_bmdm_genes[rowSums(my.bmdm.mat) == 2]
print(my.bmdm.expressed, quote = F)

# Export all shared DEGs to txt
write.table(my.bmdm.expressed, file = "bmdm_shared_DEGS_50.txt", sep = "\t" , row.names = F, quote=F)

#########################
# Divide up and downregulation for over representation analysis

bmdm_dfs <- list( "BMDM_Cohort_1" = BMDM_Cohort1, 
                  "BMDM_Cohort_2" = BMDM_Cohort2)

bmdm_expressed <- list( "BMDM_Cohort_1" = BMDM_Cohort1_$row.names, 
                   "BMDM_Cohort_2" = BMDM_Cohort2_$row.names)

all_bmdm_expressed <- unique(unlist(bmdm_expressed))

# Create matrix
my.fc.bmdm.mat                 <- data.frame(matrix(0,length(all_bmdm_genes),length(bmdm_DEGs)))
colnames(my.fc.bmdm.mat)       <- names(bmdm_DEGs)
rownames(my.fc.bmdm.mat)       <- all_bmdm_genes

# Populate matrix
for (h in 1:length(bmdm_DEGs)) {
  my.fc.bmdm.mat[all_bmdm_genes %in% bmdm_DEGs[[h]],h] <- 1
}
for(h in  1:length(bmdm_DEGs)) {
  for(g in 1:length(all_bmdm_genes)){
    if(my.fc.bmdm.mat[g,h] == 1){
      gene <- row.names(my.fc.bmdm.mat[g,])
      df <- bmdm_dfs[[h]]
      subdf <- subset(df, df$row.names == gene)
      my.fc.bmdm.mat[g,h] <- subdf$log2FoldChange
    }
  }
}

# Replace 0s with NAs
my.fc.bmdm.mat[my.fc.bmdm.mat == 0] <- NA

# Add new column with count of NAs
my.fc.bmdm.mat$Total_NAs <- rowSums(is.na(my.fc.bmdm.mat))

# Add new columns with counts of positive and negative log2fc values
my.fc.bmdm.mat$up <- rowSums(my.fc.bmdm.mat[,1:2] > 0, na.rm = TRUE)
my.fc.bmdm.mat$down <- rowSums(my.fc.bmdm.mat[,1:2] < 0, na.rm = TRUE)

my.up.bmdm <- all_bmdm_genes[my.fc.bmdm.mat$up >= 2]
my.up.bmdm
length(my.up.bmdm)
#[1] 12

my.down.bmdm <- all_bmdm_genes[my.fc.bmdm.mat$down >= 2]
my.down.bmdm
length(my.down.bmdm)
#[1] 4

# Perform ORA on upregulated genes

entrezID.genes.shared <- bitr(my.up.bmdm , fromType="SYMBOL", toType="ENTREZID", 
                              OrgDb="org.Mm.eg.db")
#'select()' returned 1:1 mapping between keys and columns
entrezID.genes.universe <- bitr(all_bmdm_expressed ,fromType="SYMBOL", toType="ENTREZID", 
                                OrgDb="org.Mm.eg.db")
#'select()' returned 1:1 mapping between keys and columns
#Warning message:
#  In bitr(all_bmdm_expressed, fromType = "SYMBOL", toType = "ENTREZID",  :
#            5.25% of input gene IDs are fail to map...

# Run over-representation analysis
ego.shared <- enrichGO(gene   = entrezID.genes.shared$ENTREZID,
                       universe  = entrezID.genes.universe$ENTREZID,
                       OrgDb         = org.Mm.eg.db,
                       keyType       = 'ENTREZID',
                       ont           = "ALL",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE,
                       minGSSize     = 10  ,
                       maxGSSize     = 5000)

bmdm_up_ora <- as.data.frame(ego.shared)
#28 rows

# Export to file
write.table(bmdm_up_ora, file = paste(Sys.Date(),"bmdm_ORA_upregulated_50.txt", sep ="_"), sep = "\t" , row.names = F, quote=F)

#sort and reduce for plot
df <- bmdm_up_ora
df <- df[order(df$qvalue),]
df <- df[1:10,]

# Make plot

df$GeneRatioNumeric <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), `[`, 1)) / 
  as.numeric(sapply(strsplit(df$GeneRatio, "/"), `[`, 2))

# Create a base plot
plot <- ggplot(df, aes(x = 1, y = Description, size = GeneRatioNumeric)) + 
  geom_point(alpha = 0.6) +
  scale_size_continuous(name = "Gene Ratio")

#Set color scale
plot <- plot + aes(color = -log10(pvalue)) + 
  scale_color_gradient2(low = "white", high = "red", 
                        name = "-log10(p-value)" , limits=c(0,5))


# Final touches
plot <- plot + theme(axis.title.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.x = element_blank(),
                     panel.grid.major.y = element_line(color="#D3D3D3"),
                     panel.background = element_blank(),
                     panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(title = "BMDM UpRegulated ORA Bubble Plot")

pdf(paste(Sys.Date(),"bmdm_upregulated_DEGs_bubbleplot.pdf", 
          sep = "_"))
print(plot)
dev.off()


# Perform ORA on downregulated genes

entrezID.genes.shared <- bitr(my.down.bmdm , fromType="SYMBOL", toType="ENTREZID", 
                              OrgDb="org.Mm.eg.db")
#'select()' returned 1:1 mapping between keys and columns
entrezID.genes.universe <- bitr(all_bmdm_expressed ,fromType="SYMBOL", toType="ENTREZID", 
                                OrgDb="org.Mm.eg.db")
#'select()' returned 1:1 mapping between keys and columns
#Warning message:
#  In bitr(all_bmdm_expressed, fromType = "SYMBOL", toType = "ENTREZID",  :
#            5.25% of input gene IDs are fail to map...
# Run over-representation analysis
ego.shared <- enrichGO(gene   = entrezID.genes.shared$ENTREZID,
                       universe  = entrezID.genes.universe$ENTREZID,
                       OrgDb         = org.Mm.eg.db,
                       keyType       = 'ENTREZID',
                       ont           = "ALL",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE,
                       minGSSize     = 10  ,
                       maxGSSize     = 5000)

bmdm_down_ora <- as.data.frame(ego.shared)
#16 rows

# Export to file
write.table(bmdm_down_ora, file = paste(Sys.Date(),"bmdm_ORA_downregulated_50.txt", sep ="_"), sep = "\t" , row.names = F, quote=F)

#sort and reduce for plot
df <- bmdm_down_ora
df <- df[order(df$qvalue),]
df <- df[1:10,]

df$GeneRatioNumeric <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), `[`, 1)) / 
  as.numeric(sapply(strsplit(df$GeneRatio, "/"), `[`, 2))

# Create a base plot
plot <- ggplot(df, aes(x = 1, y = Description, size = GeneRatioNumeric)) + 
  geom_point(alpha = 0.6) +
  scale_size_continuous(name = "Gene Ratio")

#Set color scale
plot <- plot + aes(color = -log10(pvalue)) + 
  scale_color_gradient2(low = "white", high = "blue", 
                        name = "-log10(p-value)" , limits=c(0,5))


# Final touches
plot <- plot + theme(axis.title.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.x = element_blank(),
                     panel.grid.major.y = element_line(color="#D3D3D3"),
                     panel.background = element_blank(),
                     panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(title = "BMDM DownRegulated ORA Bubble Plot") +
  scale_y_discrete(labels = function(y) str_wrap(y, width = 10)) #for long descriptions, wrap the text

pdf(paste(Sys.Date(),"bmdm_downregulated_DEGs_bubbleplot.pdf", 
          sep = "_"))
print(plot)
dev.off()

###################################3
#OCP
ocp_DEGs <- list("GSE153299_OCP_MCSF" = GSE153299_OCP_MCSF$row.names,
                 "GSE153299_OCP_MCSF_RANKL" = GSE153299_OCP_MCSF_RANKL$row.names)

all_ocp_genes <- unique(unlist(ocp_DEGs))
length(all_ocp_genes)
#[1] 148

# Create matrix
my.ocp.mat                 <- data.frame(matrix(0,length(all_ocp_genes),length(ocp_DEGs)))
colnames(my.ocp.mat)       <- names(ocp_DEGs)
rownames(my.ocp.mat)       <- all_ocp_genes

# Populate matrix
for (h in 1:length(ocp_DEGs)) {
  my.ocp.mat[all_ocp_genes %in% ocp_DEGs[[h]],h] <- 1
}

dim(my.ocp.mat)
#[1] 148    2


# Find genes shared across all datasets
sum(rowSums(my.ocp.mat) == 2)
#[1] 17
my.ocp.expressed <- all_ocp_genes[rowSums(my.ocp.mat) == 2]
print(my.ocp.expressed, quote = F)

# Export all shared DEGs to txt
write.table(my.ocp.expressed, file = "ocp_shared_DEGS_50.txt", sep = "\t" , row.names = F, quote=F)

#########################
# Divide up and downregulation for over representation analysis

ocp_dfs <- list("GSE153299_OCP_MCSF" = GSE153299_OCP_MCSF,
                "GSE153299_OCP_MCSF_RANKL" = GSE153299_OCP_MCSF_RANKL)

ocp_expressed <- list("GSE153299_OCP_MCSF" = GSE153299_OCP_MCSF_$row.names,
                 "GSE153299_OCP_MCSF_RANKL" = GSE153299_OCP_MCSF_RANKL_$row.names)

all_ocp_expressed <- unique(unlist(ocp_expressed))

# Create matrix
my.fc.ocp.mat                 <- data.frame(matrix(0,length(all_ocp_genes),length(ocp_DEGs)))
colnames(my.fc.ocp.mat)       <- names(ocp_DEGs)
rownames(my.fc.ocp.mat)       <- all_ocp_genes

# Populate matrix
for (h in 1:length(ocp_DEGs)) {
  my.fc.ocp.mat[all_ocp_genes %in% ocp_DEGs[[h]],h] <- 1
}
for(h in  1:length(ocp_DEGs)) {
  for(g in 1:length(all_ocp_genes)){
    if(my.fc.ocp.mat[g,h] == 1){
      gene <- row.names(my.fc.ocp.mat[g,])
      df <- ocp_dfs[[h]]
      subdf <- subset(df, df$row.names == gene)
      my.fc.ocp.mat[g,h] <- subdf$log2FoldChange
    }
  }
}

# Replace 0s with NAs
my.fc.ocp.mat[my.fc.ocp.mat == 0] <- NA

# Add new column with count of NAs
my.fc.ocp.mat$Total_NAs <- rowSums(is.na(my.fc.ocp.mat))

# Add new columns with counts of positive and negative log2fc values
my.fc.ocp.mat$up <- rowSums(my.fc.ocp.mat[,1:2] > 0, na.rm = TRUE)
my.fc.ocp.mat$down <- rowSums(my.fc.ocp.mat[,1:2] < 0, na.rm = TRUE)

my.up.ocp <- all_ocp_genes[my.fc.ocp.mat$up > 1 ]
my.up.ocp
length(my.up.ocp)
#[1] 13

my.down.ocp <- all_ocp_genes[my.fc.ocp.mat$down > 1 ]
my.down.ocp
length(my.down.ocp)
#[1] 4

# Perform ORA on upregulated genes

entrezID.genes.shared <- bitr(my.up.ocp , fromType="SYMBOL", toType="ENTREZID", 
                              OrgDb="org.Mm.eg.db")
#'select()' returned 1:1 mapping between keys and columns
entrezID.genes.universe <- bitr(all_ocp_expressed ,fromType="SYMBOL", toType="ENTREZID", 
                                OrgDb="org.Mm.eg.db")
#'select()' returned 1:1 mapping between keys and columns
#Warning message:
#  In bitr(all_ocp_expressed, fromType = "SYMBOL", toType = "ENTREZID",  :
#            4.71% of input gene IDs are fail to map...

# Run over-representation analysis
ego.shared <- enrichGO(gene   = entrezID.genes.shared$ENTREZID,
                       universe  = entrezID.genes.universe$ENTREZID,
                       OrgDb         = org.Mm.eg.db,
                       keyType       = 'ENTREZID',
                       ont           = "ALL",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE,
                       minGSSize     = 10  ,
                       maxGSSize     = 5000)

ocp_up_ora <- as.data.frame(ego.shared)
#41 rows

# Export to file
write.table(ocp_up_ora, file = paste(Sys.Date(),"ocp_ORA_upregulated_50.txt", sep ="_"), sep = "\t" , row.names = F, quote=F)

#sort and reduce for plot
df <- ocp_up_ora
df <- df[order(df$qvalue),]
df <- df[1:10,]

# Make plot
df$GeneRatioNumeric <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), `[`, 1)) / 
  as.numeric(sapply(strsplit(df$GeneRatio, "/"), `[`, 2))

# Create a base plot
plot <- ggplot(df, aes(x = 1, y = Description, size = GeneRatioNumeric)) + 
  geom_point(alpha = 0.6) +
  scale_size_continuous(name = "Gene Ratio")

#Set color scale
plot <- plot + aes(color = -log10(pvalue)) + 
  scale_color_gradient2(low = "white", high = "red", 
                        name = "-log10(p-value)" , limits=c(0,10))


# Final touches
plot <- plot + theme(axis.title.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.x = element_blank(),
                     panel.grid.major.y = element_line(color="#D3D3D3"),
                     panel.background = element_blank(),
                     panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(title = "OCP UpRegulated ORA Bubble Plot")

pdf(paste(Sys.Date(),"ocp_upregulated_DEGs_bubbleplot.pdf", 
          sep = "_"), width = 5, height = 5)
print(plot)
dev.off()


# Perform ORA on downregulated genes

entrezID.genes.shared <- bitr(my.down.ocp , fromType="SYMBOL", toType="ENTREZID", 
                              OrgDb="org.Mm.eg.db")
#'select()' returned 1:1 mapping between keys and columns
entrezID.genes.universe <- bitr(all_ocp_expressed ,fromType="SYMBOL", toType="ENTREZID", 
                                OrgDb="org.Mm.eg.db")
#'select()' returned 1:1 mapping between keys and columns
#Warning message:
#  In bitr(all_ocp_expressed, fromType = "SYMBOL", toType = "ENTREZID",  :
#            4.71% of input gene IDs are fail to map...

# Run over-representation analysis
ego.shared <- enrichGO(gene   = entrezID.genes.shared$ENTREZID,
                       universe  = entrezID.genes.universe$ENTREZID,
                       OrgDb         = org.Mm.eg.db,
                       keyType       = 'ENTREZID',
                       ont           = "ALL",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE,
                       minGSSize     = 10  ,
                       maxGSSize     = 5000)

ocp_down_ora <- as.data.frame(ego.shared)
#16 rows

# Export to file
write.table(ocp_down_ora, file = paste(Sys.Date(),"ocp_ORA_downregulated_50.txt", sep ="_"), sep = "\t" , row.names = F, quote=F)

#sort and reduce for plot
df <- ocp_down_ora
df <- df[order(df$qvalue),]
df <- df[1:10,]

# Make plot
df$GeneRatioNumeric <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), `[`, 1)) / 
  as.numeric(sapply(strsplit(df$GeneRatio, "/"), `[`, 2))

# Create a base plot
plot <- ggplot(df, aes(x = 1, y = Description, size = GeneRatioNumeric)) + 
  geom_point(alpha = 0.6) +
  scale_size_continuous(name = "Gene Ratio")

#Set color scale
plot <- plot + aes(color = -log10(pvalue)) + 
  scale_color_gradient2(low = "white", high = "blue", 
                        name = "-log10(p-value)" , limits=c(0,5))


# Final touches
plot <- plot + theme(axis.title.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.x = element_blank(),
                     panel.grid.major.y = element_line(color="#D3D3D3"),
                     panel.background = element_blank(),
                     panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(title = "OCP DownRegulated ORA Bubble Plot")

pdf(paste(Sys.Date(),"ocp_downregulated_DEGs_bubbleplot.pdf", 
          sep = "_"))
print(plot)
dev.off()

##################################
#Exudate
#only one dataset here

ex_DEGs <- list(GSE109099_Exudate$row.names)

all_ex_genes <- unlist(ex_DEGs)
length(all_ex_genes)
#[1] 1323

#########################
# Divide up and downregulation for over representation analysis

ex_dfs <-  list("GSE109099_Exudate" = GSE109099_Exudate)

all_ex_expressed <- unlist(GSE109099_Exudate_$row.names)
length(all_ex_expressed)
#10981


# Create matrix
my.fc.ex.mat                 <- data.frame(matrix(GSE109099_Exudate$log2FoldChange))
row.names(my.fc.ex.mat) <- GSE109099_Exudate$row.names

# Add new columns with counts of positive and negative log2fc values
my.fc.ex.mat$up <- as.numeric(my.fc.ex.mat[,1] > 0)
my.fc.ex.mat$down <- as.numeric(my.fc.ex.mat[,1] < 0)

my.up.ex <- all_ex_genes[my.fc.ex.mat$up == 1 ]
my.up.ex
length(my.up.ex)
#[1] 741

my.down.ex <- all_ex_genes[my.fc.ex.mat$down == 1 ]
my.down.ex
length(my.down.ex)
#[1] 582

# Perform ORA on upregulated genes

entrezID.genes.shared <- bitr(my.up.ex , fromType="SYMBOL", toType="ENTREZID", 
                              OrgDb="org.Mm.eg.db")
#'select()' returned 1:1 mapping between keys and columns
#Warning message:
#  In bitr(my.up.ex, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db") :
 # 4.18% of input gene IDs are fail to map...
entrezID.genes.universe <- bitr(all_ex_expressed ,fromType="SYMBOL", toType="ENTREZID", 
                                OrgDb="org.Mm.eg.db")
#'select()' returned 1:1 mapping between keys and columns
#Warning message:
#  In bitr(all_ex_expressed, fromType = "SYMBOL", toType = "ENTREZID",  :
#            4.81% of input gene IDs are fail to map...

# Run over-representation analysis
ego.shared <- enrichGO(gene   = entrezID.genes.shared$ENTREZID,
                       universe  = entrezID.genes.universe$ENTREZID,
                       OrgDb         = org.Mm.eg.db,
                       keyType       = 'ENTREZID',
                       ont           = "ALL",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE,
                       minGSSize     = 10  ,
                       maxGSSize     = 5000)

ex_up_ora <- as.data.frame(ego.shared)
#empty

# Perform ORA on downregulated genes

entrezID.genes.shared <- bitr(my.down.ex , fromType="SYMBOL", toType="ENTREZID", 
                              OrgDb="org.Mm.eg.db")
#'select()' returned 1:1 mapping between keys and columns
#Warning message:
#  In bitr(my.down.ex, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db") :
#  4.64% of input gene IDs are fail to map...
entrezID.genes.universe <- bitr(all_ex_expressed ,fromType="SYMBOL", toType="ENTREZID", 
                                OrgDb="org.Mm.eg.db")
#'select()' returned 1:1 mapping between keys and columns
#Warning message:
#  In bitr(all_ex_expressed, fromType = "SYMBOL", toType = "ENTREZID",  :
#            4.81% of input gene IDs are fail to map...


# Run over-representation analysis
ego.shared <- enrichGO(gene   = entrezID.genes.shared$ENTREZID,
                       universe  = entrezID.genes.universe$ENTREZID,
                       OrgDb         = org.Mm.eg.db,
                       keyType       = 'ENTREZID',
                       ont           = "ALL",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE,
                       minGSSize     = 10  ,
                       maxGSSize     = 5000)

ex_down_ora <- as.data.frame(ego.shared)
#10 rows

# Export to file
write.table(ex_down_ora, file = paste(Sys.Date(),"exudate_ORA_downregulated_50.txt", sep ="_"), sep = "\t" , row.names = F, quote=F)

#sort and reduce for plot
df <- ex_down_ora
df <- df[order(df$qvalue),]

# Make plot
df$GeneRatioNumeric <- as.numeric(sapply(strsplit(df$GeneRatio, "/"), `[`, 1)) / 
  as.numeric(sapply(strsplit(df$GeneRatio, "/"), `[`, 2))

# Create a base plot
plot <- ggplot(df, aes(x = 1, y = Description, size = GeneRatioNumeric)) + 
  geom_point(alpha = 0.6) +
  scale_size_continuous(name = "Gene Ratio")

#Set color scale
plot <- plot + aes(color = -log10(pvalue)) + 
  scale_color_gradient2(low = "white", high = "blue", 
                        name = "-log10(p-value)" , limits=c(0,5))


# Final touches
plot <- plot + theme(axis.title.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.x = element_blank(),
                     panel.grid.major.y = element_line(color="#D3D3D3"),
                     panel.grid.major.x = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
                     panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(title = "Exudate DownRegulated ORA Bubble Plot")

pdf(paste(Sys.Date(),"exudate_downregulated_DEGs_bubbleplot.pdf", 
          sep = "_"))
print(plot)
dev.off()


################################################################################################
######################## 3. Export list of shared DEGs from each body type for venn diagram

#setwd("C:/Users/livis/Documents/Benayoun_Lab/deg_comparisons_50/DEG_lists")
write.table(GSE109099_Exudate$row.names, file="exudate_DEGs_list.txt", row.names = F, quote = F)
write.table(my.peritoneal.expressed, file="peritoneal_DEGs_list.txt", row.names = F, quote = F)
write.table(my.microglia.expressed, file="microglia_DEGs_list.txt", row.names = F, quote = F)
write.table(my.alveolar.expressed, file="alveolar_DEGs_list.txt", row.names = F, quote = F)
write.table(my.ocp.expressed, file="ocp_DEGs_list.txt", row.names = F, quote = F)
write.table(my.bmdm.expressed, file="bmdm_DEGs_list.txt", row.names = F, quote = F)


#######################
sink(file = paste(Sys.Date(),"SMAC_ORA_RNAseq_analysis_session_Info.txt", sep =""))
sessionInfo()
sink()

