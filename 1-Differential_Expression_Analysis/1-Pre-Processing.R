# 2023-4-26
# Take raw STAR counts, run SVA, normalize with DESeq2, perform MDS and PCA analysis, and gather 
# differentially expressed genes with FDR < 0.05

################################################################################################
######################## 1. Load in necessary data and packages

options(stringsAsFactors = F)
library(DESeq2)
library(pheatmap)
library('pvclust')
library('bitops')
library('sva')
library('limma')
library(RColorBrewer)
library(bitops)
library(UpSetR)
library("ComplexHeatmap")
library(gplots)
library(circlize)
library(grid)
library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Hs.eg.db)


#read in raw datasets

setwd("C:/Users/livis/Dropbox/Olivia_SexDim_Macrophage/STAR_count_Matrices")
BMDM_Cohort1 <- read.csv("2021-06-02_BenayounLab_BMDM_Cohort1_SEX_counts.txt", sep = "\t", header = T, skip = 1)
BMDM_Cohort2 <- read.csv("2021-06-02_BenayounLab_BMDM_Cohort2_SEX_counts.txt" ,sep = "\t", header = T, skip = 1)
Peritoneal_Cohort2 <- read.csv("2021-06-02_BenayounLab_Peritoneal_Macrophages_Cohort2_SEX_counts.txt" , sep = "\t", header = T, skip = 1)
GSE41879_Peritoneal <- read.csv("2021-06-02_GSE41879_Peritoneal_macrophages_SEX_counts.txt", sep = "\t", header = T, skip = 1)
GSE99622_Microglia_whole_brain_P14 <- read.csv("2021-06-02_GSE99622_microglia_whole_brain_P14_SEX_counts.txt", sep = "\t", header = T, skip = 1)
GSE124829_Microglia_whole_brain_ImmGen_8w <- read.csv("2021-06-02_GSE124829_microglia_whole_brain_ImmGen_8w_SEX_counts.txt", sep = "\t", header = T, skip = 1)
GSE124829_Peritoneal_ImmGen_6w <- read.csv("2021-06-02_GSE124829_Peritoneal_Macrophages_ImmGen_6w_SEX_counts.txt", sep = "\t", header = T, skip = 1)
GSE124829_Spleen_ImmGen_6w <- read.csv("2021-06-02_GSE124829_Spleen_Macrophages_ImmGen_6w_SEX_counts.txt", sep = "\t", header = T, skip = 1)
GSE149014_Peritoneal <- read.csv("2021-06-02_GSE149014_Peritoneal_Macrophages_SEX_counts.txt", sep = "\t", header = T, skip = 1)
GSE149014_Pleural <- read.csv("2021-06-02_GSE149014_Pleural_Macrophages_SEX_counts.txt", sep = "\t", header = T, skip = 1)
GSE153299_OCP_MCSF_RANKL <- read.csv("2021-06-02_GSE153299_OCP_MCSF_RANKL_SEX_counts.txt", sep = "\t", header = T, skip = 1)
GSE153299_OCP_MCSF <- read.csv("2021-06-02_GSE153299_OCP_MCSF_SEX_counts.txt", sep = "\t", header = T, skip = 1)
PRJNA383777_Microglia_whole_brain <- read.csv("2021-06-02_PRJNA383777_microglia_whole_brain_SEX_counts.txt", sep = "\t", header = T, skip = 1)
PRJNA408225_Microglia_frontal_lobe <- read.csv("2021-06-02_PRJNA408225_microglia_Frontal_Lobe_SEX_counts.txt", sep = "\t", header = T, skip = 1)
PRJNA408225_Microglia_hippocampus <- read.csv("2021-06-02_PRJNA408225_microglia_Hippocampus_SEX_counts.txt", sep = "\t", header = T, skip = 1)

setwd("C:/Users/livis/Dropbox/Olivia_SexDim_Macrophage/New_count_matrices")
GSE174207_Alveolar_CC <- read.csv("2022-03-22_GSE174207_Alveolar_Macrophages_CC_SEX_counts.txt", sep = "\t", header = T, skip = 1)
GSE156799_Alveolar <- read.csv("2022-03-22_GSE156799_Alveolar_Macrophages_SEX_counts.txt", sep = "\t", header = T, skip = 1)
GSE109099_Exudate <- read.csv("2022-03-22_GSE109099_Exudate_Macrophages_SEX_counts.txt", sep = "\t", header = T, skip = 1)
GSE109099_Alveolar <- read.csv("2022-03-22_GSE109099_Alveolar_Macrophages_SEX_counts.txt", sep = "\t", header = T, skip = 1)

setwd("C:/Users/livis/Dropbox/2023_Sex_Diff_macrophages_Cassie_Olivia/STAR_Files/DS")
Peritoneal_Cohort1 <- read.csv("2023-04-20_BenayounLab_Peritoneal_Macrophages_Cohort1_SEX_DS_CLEAN_counts.txt", sep = "\t", header = T, skip = 1)
GSE99622_Microglia_whole_brain_P60 <- read.csv("2023-04-20_GSE99622_microglia_whole_brain_P60_SEX_DS_CLEAN_counts.txt", sep = "\t", header = T, skip = 1)

#21 total datasets

#put into list
my.raw.data.list <- list("Benayoun_BMDM_Cohort1" = BMDM_Cohort1, 
                         "Benayoun_BMDM_Cohort2" = BMDM_Cohort2, 
                         "Benayoun_Peritoneal_Cohort1" = Peritoneal_Cohort1, 
                         "Benayoun_Peritoneal_Cohort2" = Peritoneal_Cohort2, 
                         "GSE41879_Peritoneal" = GSE41879_Peritoneal, 
                         "GSE99622_Microglia_whole_brain_P14" = GSE99622_Microglia_whole_brain_P14,
                         "GSE99622_Microglia_whole_brain_P60" = GSE99622_Microglia_whole_brain_P60, 
                         "GSE124829_Microglia_whole_brain_ImmGen_8w" = GSE124829_Microglia_whole_brain_ImmGen_8w, 
                         "GSE124829_Peritoneal_ImmGen_6w" = GSE124829_Peritoneal_ImmGen_6w, 
                         "GSE124829_Spleen_ImmGen_6w" = GSE124829_Spleen_ImmGen_6w, 
                         "GSE149014_Peritoneal" = GSE149014_Peritoneal, 
                         "GSE149014_Pleural" = GSE149014_Pleural, 
                         "GSE153299_OCP_MCSF_RANKL" = GSE153299_OCP_MCSF_RANKL, 
                         "GSE153299_OCP_MCSF" = GSE153299_OCP_MCSF, 
                         "PRJNA383777_Microglia_whole_brain" = PRJNA383777_Microglia_whole_brain, 
                         "PRJNA408225_Microglia_frontal_lobe" = PRJNA408225_Microglia_frontal_lobe, 
                         "PRJNA408225_Microglia_hippocampus" = PRJNA408225_Microglia_hippocampus,
                         "GSE174207_Alveolar_CC" = GSE174207_Alveolar_CC, 
                         "GSE156799_Alveolar" = GSE156799_Alveolar, 
                         "GSE109099_Alveolar" = GSE109099_Alveolar, 
                         "GSE109099_Exudate" = GSE109099_Exudate)



################################################################################################
######################## 2. Process each dataset for further analysis

#set up new wd
setwd("C:/Users/livis/Documents/Benayoun_Lab/11-5_rerun/preprocess2")

#set up lists to be added to (just as an added way to check values)
my.percent.difs <- list()
my.degs <- list()

#adjust the code to work for your working directories

#run preprocessing for each dataset
for(i in 1:21){
  my.outprefix <-paste(Sys.Date(),names(my.raw.data.list[i]), sep="_")
  print(my.outprefix)
  curr_name <- names(my.raw.data.list[i])
  my.initial <- my.raw.data.list[[i]]
  
  #remove extra data
  my.data <- my.initial[ , -c(2:6)]
  
  
  #look for females
  femalecols <- grep("_female" , colnames(my.data), value = FALSE)
  if (length(femalecols) == 0)
    femalecols <- grep("_Female", colnames(my.data), value = FALSE)
  if (length(femalecols) == 0)
    femalecols <- grep("_F", colnames(my.data), value = FALSE)
  
  #rename female columns
  q = 1
  for (j in femalecols){
    colnames(my.data)[j] <- paste0("Female" , q)
    q <- q + 1
  }
  
  #look for males (cannot search for _M as most colnames contain _Macrophages)
  countmales = length(my.data) - length(femalecols) - 1
  malecols = c(2:length(my.data))
  malecols = subset(malecols, !(malecols %in% femalecols))
  
  #rename male columns
  p = 1
  for (k in malecols) {
    colnames(my.data)[k] <- paste0("Male", p)
    p <- p + 1
  }
  
  #get the genes with no reads out
  my.mp3 <- my.data[rowSums(my.data[,-1])>0,]
  my.mp3 <- my.mp3[,-1]
  my.goodlist <- which(rowSums(my.data[,-1])>0)
  rownames(my.mp3) <- my.data$Geneid[my.goodlist]
  
  #fix male and female cols
  malecols = malecols - 1
  femalecols = femalecols - 1
  
  #############################################################################
  #first round of QC here - filtering out datasets with read depth discrepancies
  #need to sum read depth for each sex then find percent difference
  
  male.rd <- colSums(my.mp3[malecols])
  female.rd <- colSums(my.mp3[femalecols])
  avg.male.rd <- sum(male.rd) / length(malecols)
  avg.female.rd <- sum(female.rd) / length(femalecols)
  percent.dif <- (avg.female.rd - avg.male.rd) / avg.female.rd
  
  #remove percent difference discrepancies
  if(percent.dif < -0.5 || percent.dif > 0.5){
    print(paste(names(my.raw.data.list[i]), " excluded for discrepancy in read depth between males and females"))
    next
  }
  
  #if not removed, add to list of read depths
  my.percent.difs[curr_name] <- percent.dif
  
  #remove all individuals with read depth < 1 mil
  my.best <- colSums(my.mp3) > 1000000 
  my.filtered.matrix <- my.mp3[,my.best] 
  print(colnames(my.filtered.matrix))
  
  #make sex array
  filtered.arr <- grep("Female" , colnames(my.filtered.matrix), value = TRUE)
  num.males <- length(my.filtered.matrix) - length(filtered.arr)
  my.Sex <- colnames(my.filtered.matrix)
  for(m in 1:length(my.filtered.matrix)){
    if(grepl('Female',my.Sex[m])){
      my.Sex[m] <- "F"
    }
    else{
      my.Sex[m] <- "M"
    }
  }
  print(my.Sex)
  
  
  
  ###########################################################################
  # Run SVA to remove unwanted variation
  
  
  # build design matrix
  dataDesign = data.frame( row.names = colnames( my.filtered.matrix ), 
                           sex = my.Sex)
  
  # Set null and alternative models (ignore batch)
  mod1 = model.matrix(~ sex, data = dataDesign)
  n.sv.be = num.sv(my.filtered.matrix,mod1,method="be") 
  
  #check n.sv.be
  n.sv.be
  
  if (n.sv.be > 0){
    #apply SVAseq algorithm
    my.svseq = svaseq(as.matrix(my.filtered.matrix), mod1, n.sv=n.sv.be, constant = 0.1)
    
    # remove RIN and SV, preserve sex
    my.clean <- removeBatchEffect(log2(my.filtered.matrix + 0.1), 
                                  batch=NULL, 
                                  covariates=cbind(my.svseq$sv),
                                  design=mod1)
    
    # delog and round data for DEseq2 processing
    my.filtered.sva <- round(2^my.clean-0.1)
    
    write.table(my.filtered.sva, file = paste(my.outprefix,"postSVA_Counts.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)
    
  }
  if (n.sv.be == 0) {
    my.filtered.sva <- (my.filtered.matrix)
    
  }
  ############################################################################
  #DESeq2 on cleaned data
  
  # design matrix
  dataDesignDE = data.frame( row.names = colnames( my.filtered.sva ), 
                             sex = my.Sex)
  
  # get matrix
  dds <- DESeqDataSetFromMatrix(countData = my.filtered.sva,
                                colData = dataDesign,
                                design = ~ sex)
  
  # run DESeq normalizations and export results
  dds.deseq <- DESeq(dds)
  
  # plot dispersion
  my.disp.out <- paste(my.outprefix,"dispersion_plot.pdf",sep="_")
  
  pdf(my.disp.out)
  plotDispEsts(dds.deseq)
  dev.off()
  
  # normalized expression value
  tissue.cts <- getVarianceStabilizedData(dds.deseq)
  
  #get file
  write.table(tissue.cts, file = paste(my.outprefix,"normalized_expression.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)
  
  
  # color-code 
  my.colors <- c(1:2)
  my.colors[grep("Female",colnames(tissue.cts))] <- "deeppink"
  my.colors[grep("Male",colnames(tissue.cts))]    <- "deepskyblue"
  
  # do MDS analysis
  mds.result <- cmdscale(1-cor(tissue.cts,method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
  x <- mds.result[, 1]
  y <- mds.result[, 2]
  
  #setwd("C:/Users/livis/Documents/Benayoun_Lab/11-5_rerun/preprocess2/mds")
  
  my.mds.out <- paste(my.outprefix,"MDS_plot.pdf",sep="_")
  pdf(my.mds.out)
  plot(x, y,
       xlab = "MDS dimension 1", ylab = "MDS dimension 2",
       main="Multi-dimensional Scaling",
       cex=3, pch = 16, col = my.colors,
       xlim = c(-0.06,0.06),
       ylim = c(-0.06,0.06),
       cex.lab = 1.5,
       cex.axis = 1.5)
  dev.off()
  
  #setwd("C:/Users/livis/Documents/Benayoun_Lab/11-5_rerun/preprocess2")
  
  # expression range
  pdf(paste(my.outprefix,"_Normalized_counts_boxplot.pdf", sep="_"))
  boxplot(tissue.cts,col=my.colors,cex=0.5,ylab="Log2 DESeq2 Normalized counts", las = 2)  
  dev.off()
  
  #setwd("C:/Users/livis/Documents/Benayoun_Lab/11-5_rerun/preprocess2/xist-ddx3y")
  
  # plot Xist and Ddx3y expression 
  pdf(paste(my.outprefix,"_Normalized_counts_Xist_DdX3y_scatter.pdf" , sep="_"))
  plot(
    tissue.cts["Xist",], tissue.cts["Ddx3y",],
    xlab = "Normalized log2(counts) Xist expression",
    ylab = "Normalized log2(counts) Ddx3y expression",
    cex=3, pch = 16, col = my.colors,
    cex.lab = 1.5,
    cex.axis = 1.5
  )
  dev.off()
  #setwd("C:/Users/livis/Documents/Benayoun_Lab/11-5_rerun/preprocess2")
  
  ################################################################################
  #Round 2 of QC - must visually inspect above graphs to ensure proper labeling of males and females
  
  ###############################################################################
  ## c. sex with age as covariate
  res.sex <- results(dds.deseq, contrast = c("sex","F","M")) # FC in females over Males
  
  ### get the genes with sex dimorphic changes at FDR5; exclude NA
  res.sex <- res.sex[!is.na(res.sex$padj),]
  
  genes.sex <- rownames(res.sex)[res.sex$padj < 0.05]
  my.num.sex <- length(genes.sex)
  
  save(res.sex, file = paste(my.outprefix,"_SEX.RData", sep ="_"))
  
  
  # output result tables of combined analysis to text files
  my.out.ct.mat <- paste(my.outprefix,"_log2_counts_matrix_DEseq2_SVA.txt",sep = "_")
  write.table(tissue.cts, file = my.out.ct.mat , sep = "\t" , row.names = T, quote=F)
  
  my.out.stats.sex <- paste(my.outprefix,"SEX_DIM_all_genes_statistics.txt",sep = "_")
  write.table(res.sex, file = my.out.stats.sex , sep = "\t" , row.names = T, quote=F)
  
  my.out.fdr5.sex <- paste(my.outprefix,"SEX_DIM_FDR5_genes_statistics.txt",sep = "_")
  write.table(res.sex[genes.sex,], file = my.out.fdr5.sex, sep = "\t" , row.names = T, quote=F)
  ################################################################################################
  #Round 3 of QC - exclude all datasets with fewer than 50 differentially expressed genes
  degs <- length(res.sex[genes.sex,1])
  if(degs < 50){
    print(paste(names(my.raw.data.list[i]), " excluded for too few differentially expressed genes: ", degs))
    next
  }
  
  #if not removed add to list of number of DEGs
  my.degs[curr_name] <- degs

}
# 18 datasets left after QC
# 21 warnings are converting characters to factors

#export list of degs
my.degs.df <- data.frame(name = names(my.degs), degs = unlist(my.degs), row.names = NULL)
write.csv(my.degs.df, "my_degs.csv", row.names = FALSE)

#######################
sink(file = paste(Sys.Date(),"SMAC_RNAseq_preprocessing_analysis_session_Info.txt", sep =""))
sessionInfo()
sink()

