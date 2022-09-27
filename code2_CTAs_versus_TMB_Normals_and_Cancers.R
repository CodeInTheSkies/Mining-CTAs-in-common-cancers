# NOTE: This script needs to be run in R version 4.2.1 or later, so that the
# latest ComplexHeatmap features can be fully used. For example, adding text on
# top of the bottomAnnotation bar plots.
#
# Here, we load the saved Rds file that has the subset expression matrix of the
# GTEx and TCGA samples only for the chosen genes of interest that are expected
# to be analyzed. If any other genes need to be analyzed, then
# "code1_preProcess_and_subset" has to be run first to read the whole raw data
# matrix and then create the Rds subset. This script can then be run to see the
# results with the new genes included.
#
# NOTE: The downloaded expression values stored in Rds are from Xena, and are
# already log2 transformed. So, the values are log2(norm_count+1) from:
# https://xenabrowser.net/datapages/?dataset=TcgaTargetGtex_RSEM_Hugo_norm_count&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
#
# The analysis involves systematically matching the clinical data for each of
# the samples, bin the cancer samples based on predefined tumor mutation burden
# (TMB) ranges or other predefined groups, and then grouping the samples at a
# high level based on normal tissues and cancers. Z-scoring for the heatmaps are
# done for the whole final matrix (just before plotting the heatmap), instead of
# the usual row-wise z-scoring. As discussed and decided with the team,
# full-matrix z-scoring is advantageous in mining and bringing the most
# interesting CTAs (genes) to the top based on their absolute expression values
# with respect to the full matrix (including all normals and chosen cancers in a
# given heatmap).
#
# For certain cancers, we may need to do things a bit differently, as for some
# cancers such as breast, we may need to further divide the patients into
# subgroups (breast cancer subtypes). So, we use clinical data accordingly to
# divide into subgroups as needed, depending on the type of cancer being
# analyzed.
#
# This code is designed such that we can do multiple cancers all in one go,
# using minor modifications.
# 

library(ggplot2)
library(ggpubr)
library(tools)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# Specifying the codes for the cancers to include. Refer to the table below for codes for all available cancers in TCGA.

# cnames <- toupper(c("PAAD")) PAAD is problematic, as not all 3 intervals have
# at least one tmb value; so pending some thinking as to how to deal with these;
# probably we can ignore the TMB-based grouping for such cases, and instead,
# simply pool all samples into one main group for the particular cancer, or
# divide using subgroups in some other clinical column
#
# We will work on dealing with such cases and release an update in a future
# version of the code
#
# cnames <- toupper(c("blca", "gbm", "luad", "skcm", "OV", "COAD", "LIHC", "STAD", "LAML", "DLBC", "KIRC"))
# 
cnames <- toupper(c("blca","brca"))
subtypeFlag="Yes"
extrastrng=""

# ---- TABLE of Cancer Codes, for reference
# -------------------------------------------
# LAML	Acute Myeloid Leukemia
# ACC   Adrenocortical carcinoma
# BLCA	Bladder Urothelial Carcinoma
# LGG	  Brain Lower Grade Glioma
# BRCA	Breast invasive carcinoma
# CESC	Cervical squamous cell carcinoma and endocervical adenocarcinoma
# CHOL	Cholangiocarcinoma
# LCML	Chronic Myelogenous Leukemia
# COAD	Colon adenocarcinoma
# CNTL	Controls
# ESCA	Esophageal carcinoma
# FPPP	FFPE Pilot Phase II
# GBM	  Glioblastoma multiforme
# HNSC	Head and Neck squamous cell carcinoma
# KICH	Kidney Chromophobe
# KIRC	Kidney Clear Cell Carcinoma
# KIRP	Kidney renal papillary cell carcinoma
# LIHC	Liver hepatocellular carcinoma
# LUAD	Lung adenocarcinoma
# LUSC	Lung squamous cell carcinoma
# DLBC	Diffuse Large B-cell Lymphoma
# MESO	Mesothelioma
# MISC	Miscellaneous
# OV	  Ovarian serous cystadenocarcinoma
# PAAD	Pancreatic adenocarcinoma
# PCPG	Pheochromocytoma and Paraganglioma
# PRAD	Prostate adenocarcinoma
# READ	Rectum adenocarcinoma
# SARC	Sarcoma
# SKCM	Skin Cutaneous Melanoma
# STAD	Stomach adenocarcinoma
# TGCT	Testicular Germ Cell Tumors
# THYM	Thymoma
# THCA	Thyroid carcinoma
# UCS	  Uterine Carcinosarcoma
# UCEC	Uterine Corpus Endometrial Carcinoma
# UVM	  Uveal Melanoma

# The following line automatically gets the current script's directory, and sets
# it to the current working directory. However, this may work only within
# RStudio's editor. Once the working directory is successfully set using this
# line, the relative paths below will work
cwd <- dirname(rstudioapi::getActiveDocumentContext()$path)
cat("The current directory of the source script is: ", cwd)
setwd(cwd)

# Read the subset Rds matrix of normalized TOIL-reCompute expressions
submat <- readRDS(file="./bigdata/subset_GTEx_TCGA_expressions.Rds")

# Replace all dots by hyphens for consistency of sample IDs with clinical phenotype data
colnames(submat) <- gsub(".", "-", colnames(submat), fixed=TRUE)

# Reading the phenotype data
phenotype <- read.delim("./smalldata/TcgaTargetGTEX_phenotype.txt")

# Excluding all TARGET entries from the Phenotype data
phenotype <- dplyr::filter(phenotype, X_study == "TCGA" | X_study == "GTEX")
row.names(phenotype) <- phenotype$sample

# Creating a clinical data matrix for the subset expr data, with matching order
clindat <- phenotype[colnames(submat), ]

# Checking if rownames of clindat and colnames of submat perfectly match
cat("\nChecking if rownames of clindat and colnames of submat perfectly match each other:", identical(row.names(clindat), colnames(submat)), "\n")
cat("Dimensions:\n")
cat("Subset matrix, submat (genes x samples):", dim(submat)[1], "x", dim(submat)[2])
cat("Clinical data, clindat (samples x clinical):", dim(clindat)[1], "x", dim(clindat)[2])

# The TMB data is downloadable in terms of individual files for each cancer. A
# single file containing all the TMB values for all cancer samples from TCGA
# could not be found for easy one-time download. Hence, we use the individual
# files corresponding to single cancers and then read them one by one and then
# match using the samples ID to add a new column in "clindat" for the TMB
# values.
#
# Main link used for download:  https://www.cbioportal.org/datasets
#
# From the above link, search for "pancan". All the cancers available will be
# listed, and the links for each cancer will lead to its corresponding data.
#
# For example, for lung, the particular link would be:
# https://www.cbioportal.org/study/summary?id=luad_tcga_pan_can_atlas_2018
#
# From this main lung link, choose the "clinical data" sub-tab and then choose
# "TMB(nonsynomymous)" and "Tumor Type" columns only. Deselect all else. Then
# download using the "down arrow" button below the "custom selection" button
# near the left of the bottom search bar (it should look like a cloud). This
# will download a tsv table as a separate, small-size text file with a few
# relevant clinical columns.
#
# Once all files are successfully downloaded, if they are placed in a subfolder
# called "clindata", then the below line would read all TMB data.
fnames <- list.files(path="./bigdata/clindata/", pattern="_clinical_data.tsv")

# Populating a grandtmb dataframe that collects all the TMB values from individual files
grandtmb <- data.frame(sample=character(), tmb=double())
for (kii in 1:length(fnames)) {
  clinical_data <- read.delim(paste0("./bigdata/clindata/", fnames[kii]), header=TRUE, stringsAsFactors = FALSE)
  subcdat <- data.frame(sample=clinical_data$Sample.ID, tmb=clinical_data$TMB..nonsynonymous.)
  row.names(subcdat) <- subcdat$sample
  subcdat <- subcdat[row.names(subcdat) %in% row.names(clindat), ]
  grandtmb <- rbind(grandtmb, subcdat)
}

# Merging grandtmb, and cleaning up
clindat <- merge(clindat, grandtmb, all.x=TRUE, by='row.names')
row.names(clindat) <- clindat$Row.names
clindat <- clindat[, !(colnames(clindat) %in% c("sample.x", "Row.names", "sample.y"))]

# For ease of use, making two more simple SampleCategory and TissueType columns
clindat$SampleCategory <- NA
clindat$TissueType <- NA
clindat$SampleCategory <- ifelse(clindat$X_study == "GTEX", clindat$X_primary_site, clindat$primary.disease.or.tissue)
clindat$TissueType <- ifelse(clindat$X_study == "GTEX", "Normal", "Cancer")

# Changing to TitleCase
clindat$SampleCategory <- tools::toTitleCase(clindat$SampleCategory)

# making sure the rownames of clindat and colnames of submat perfectly match each other
# NOTE! Since we use merge above and replace the merged output in clindat, the row order is random now!
# So, there is need to again make sure the rownames of clindat and colnames of submat perfectly match each other
clindat <- clindat[colnames(submat), ]

# Checking if rownames of clindat and colnames of submat perfectly match
cat("\nONCE AGAIN Checking if rownames of clindat and colnames of submat perfectly match each other:", identical(row.names(clindat), colnames(submat)), "\n")
cat("Dimensions:\n")
cat("Subset matrix, submat (genes x samples):", dim(submat)[1], "x", dim(submat)[2])
cat("Clinical data, clindat (samples x clinical):", dim(clindat)[1], "x", dim(clindat)[2])

# Now, at this point, we have all needed data, i.e., the needed expression values, the clinical data to use for categorization, and also the TMB values for the cancer samples
cat("\nNow, at this point, we have all needed data, i.e., the needed expression values, the clinical data to use for categorization, and also the TMB values for the cancer samples!\n")

# Read the tcga cancers abbreviations table. This table is used to decode the
# short abbreviations used above to specify the particular cancers for analysis
tcga_abbr <- read.delim(paste0("./smalldata/tcga_abbr.txt"), header=FALSE, stringsAsFactors = FALSE)
# keep only those in the table for which we have RNAseq data (only 32 out of a larger set, not all have rnaseq data)
# And, only those for which we need to see the plots at this time (cnames)
tcga_abbr <- tcga_abbr[tcga_abbr$V1 %in% cnames, ]
tcga_abbr$V2 <- tools::toTitleCase(tcga_abbr$V2)

# Now, going through the cancer names one by one, and then making the plots
for (curnam in cnames){
  
  cat("\nProcessing for cancer : ",curnam,"\n")
  
  # This is just to preserve certain variables whose values should not change within the loop
  rm(list=setdiff(ls(), c("curnam", "clindat", "clinical_data", "cnames", "cwd", "fnames", "grandtmb", "subtypeFlag", "extrastrng",
                          "phenotype", "subcdat", "submat", "tcga_abbr", "TCGA_GTEX_category", "TcgaTargetGTEX_phenotype") ))
  
  strcanType <- tcga_abbr$V2[tcga_abbr$V1 == curnam]
  
  subclin <- dplyr::filter(clindat, SampleCategory == strcanType | TissueType == "Normal")
  
  # Excluding the few rows where SampleCategory is empty
  subclin <- subclin[!(subclin$SampleCategory == ""), ]
  
  # Sorting subclin first by TissueType, then by Tissue
  subclin <- dplyr::arrange(subclin, desc(TissueType), SampleCategory)
  
  subrsem <- submat[, row.names(subclin)]
  # at this point, identical(row.names(subclin), names(subrsem)) should be TRUE
  # the rownames of subclin should be identical to column names of subrsem
  cat("\nNow inside the loop checking if rownames of subclin and colnames of subrsem perfectly match each other:", identical(row.names(subclin), colnames(subrsem)), "\n")
  
  # Setting flag for htmapdatMedian matrix initialization and pre-allocation
  intlzd=FALSE
  loopcount=0
  pltw=20
  pltht=50
  for (curg in sort(row.names(subrsem))){
    loopcount=loopcount+1
    # cat("\nDoing iteration",loopcount,"| Gene:",curg,"\n")
    hmdat <- as.data.frame(t(subrsem[curg,]))
    # Already Log2 transformed! So, no need to do again.
    names(hmdat)[names(hmdat) == curg] <- "value"
    hmdat <- merge(hmdat, subclin, by="row.names")
    row.names(hmdat) <- hmdat$Row.names
    hmdat <- hmdat[, !colnames(hmdat)=='Row.names']
    
    # Dropping any Cancer samples for which "mutation count" is NA
    hmdat <- hmdat[!(hmdat$TissueType=="Cancer" & is.na(hmdat$tmb)), ]
    
    hmdat$Group <- NA
    hmdat$Group <- ifelse(hmdat$tmb<=5, "low TMB", hmdat$Group)
    hmdat$Group <- ifelse((hmdat$tmb>5 & hmdat$tmb<10), "inbet TMB", hmdat$Group)
    hmdat$Group <- ifelse(hmdat$tmb>=10, "high TMB", hmdat$Group)
    
    hmdat$Group <- ifelse(hmdat$TissueType=="Normal", hmdat$SampleCategory, hmdat$Group)
    
    # If Breast cancer, then further divide into subtypes
    # To achieve this with minimal code, we just work on changing the labels in the "Group" column
    # Ultimately, that is what matters for the final heatmap subgroups
    if (curnam == "BRCA" & subtypeFlag=="Yes"){
      extrastrng <- "_Subtyped_"
      # First, load the subtype clinical data column
      clinical_data <- read.delim(paste0("./bigdata/clindata/brca_tcga_pan_can_atlas_2018_clinical_data.tsv"), header=TRUE, stringsAsFactors = FALSE)
      subcdat <- data.frame(sample=clinical_data$Sample.ID, Subtype=clinical_data$Subtype)
      row.names(subcdat) <- subcdat$sample
      hmdat <- merge(hmdat, subcdat, all.x=TRUE, by="row.names")
      row.names(hmdat) <- hmdat$Row.names
      hmdat <- hmdat[, !colnames(hmdat)=='Row.names']
      # dropping all rows with NAs in the Subtype column
      hmdat <- dplyr::filter(hmdat, TissueType=="Normal" | (TissueType=="Cancer" & !is.na(Subtype)))
      # Just to store the Group column to backup, before changing it
      hmdat$Group_bckup <- hmdat$Group
      hmdat$Group <- ifelse(hmdat$TissueType=="Cancer", paste0(hmdat$Subtype, " | ", hmdat$Group_bckup), hmdat$Group_bckup)
      pltw=25
      pltht=50
    }
    
    cancerGrps <- sort(unique(hmdat$Group[hmdat$TissueType=="Cancer"]), decreasing = TRUE)
    nrmlGrps <- sort(unique(hmdat$Group[hmdat$TissueType=="Normal"]))
    GroupsList <- c(nrmlGrps, cancerGrps)
    # Taking the extra general cancer category out if any, as that is split into the TMB groups
    GroupsList <- setdiff(GroupsList, strcanType)
    # 
    splitvec <- c(rep("Normal", length(nrmlGrps)), rep("Cancer", length(cancerGrps)))
    splitvec <- factor(splitvec, levels=c("Normal", "Cancer"))
    
    # pre-allocate data frame
    # median matrix
    if (!intlzd){
      htmapdatMedian <- matrix(nrow=dim(subrsem)[1], ncol=length(GroupsList))
      row.names(htmapdatMedian) <- sort(row.names(subrsem))
      colnames(htmapdatMedian) <- GroupsList
      intlzd=TRUE
    }
    
    hmdat$Group <- factor(hmdat$Group, levels=GroupsList)
    
    hmdat <- hmdat[, c("value", "tmb", "Group")]
    
    # Populate heatmap matrix with Median values
    summat <- aggregate(hmdat$value, by=list(hmdat$Group), FUN=median)
    row.names(summat) <- summat$Group.1
    summat <- summat[, !(names(summat) %in% "Group.1"), drop=F]
    for (kingg in GroupsList){
      htmapdatMedian[curg, kingg] <- summat[kingg,]
    }
  }
  
  # Number of patients in each group
  numpat <- dplyr::count(hmdat, Group)
  # row.names(numpat) <- numpat$Group
  # numpat <- numpat[, !(names(numpat) %in% "Group"), drop=T]
  numpat <- with(numpat, setNames(n, Group))
  numpat <- numpat[GroupsList] # making sure the order is right
  
  # Median Heatmap
  # currt <- t(scale(t(htmapdatMedian)))  # this line is for z-scoring by rows
  #
  # The line below does z-scoring for the entire matrix!
  currt <- (htmapdatMedian - mean(as.vector(htmapdatMedian)))/sd(as.vector(htmapdatMedian))
  currt[is.nan(currt)] <- 0 # This line deals with NaN z-scores created in rows that have no variation
  
  # Next, we eliminate all rows that have ONLY zeros
  currt <- currt[!(MatrixGenerics::rowMaxs(currt)==0 & MatrixGenerics::rowMins(currt)==0), ]
  
  # Note: just doing rowSums and eliminating all rows with rowSums==0 took out
  # too many; since there is a chance that even two equal nonzero values can
  # cancel and generate zero sums. So, above, we are looking at cases that have
  # both rowMins and rowMaxs as zero; this can happen ONLY if all values are
  # zero.  This is a better way.
  
  # # sorted by "High" group
  # ---------------------------------
  # CODE currently not being used.  Uncomment if needed.
  # currt <- currt[sort(currt[,"high | >= 20 and < 50 Mut/Mb"], decreasing = T, index.return=T)$ix,]
  #
  # col_fun = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  # hh <- ComplexHeatmap::Heatmap(currt, name = paste0("CTA expr vs TMBnonsynonymous - TCGA PanCan ",strcanType), col = col_fun, rect_gp = gpar(col= "white"),
  #                               heatmap_legend_param = list(title="Z-scored Expr"), border=T, row_names_gp = gpar(fontsize = 14),
  #                               column_title = paste0("CTA expr vs TMBnonsynonymous - TCGA PanCan\n",strcanType),
  #                               column_title_gp = gpar(fontsize = 14),
  #                               row_title = "", cluster_columns = F, cluster_rows = F, show_row_dend = F, show_column_dend = F,
  #                               column_names_rot = 30, column_names_side=c("bottom", "top"), column_names_gp = gpar(fontsize = 12))
  # pdf(file = paste(spath,"OverViewHtMap_CTAvsTMBnonsynonymous_srtdByHigh_",gsub(" ","_",strcanType,fixed=TRUE),".pdf",sep=""),
  #     useDingbats = FALSE, width=6, height=30)
  # draw(hh, padding = unit(c(12, 28, 5, 7), "mm")) #bottom, left, top, right paddings
  # dev.off()
  
  # # The below line can be uncommented if the rows need to be sorted in decreasing order of the
  # # "high" TMB group.
  # # Currently, it is commented because we are sorting rows based on the difference of the cancer avg. to the normal avg. (sans testis)
  # currt <- currt[sort(currt[,"high TMB | >= 10 Mut/Mb"], decreasing = TRUE, index.return=TRUE)$ix,]
  
  col_fun = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  
  column_ha = HeatmapAnnotation(CohortSize = anno_barplot(numpat, add_numbers = TRUE, axis_param=list(gp=gpar(fontsize = 12)),
                                                          numbers_gp=gpar(fontsize = 12, fontface="bold"), height = unit(5, "cm")))
  
  # Preparing two row annotation for boxplots for RSEM expressions
  # (1) First, rowanno1, the normal expression boxplots
  # preparing the matrix for the normals - the idea here is to leave out the testis, as
  # the testis expression is obvious and way high, which will skew other organs
  # So, we leave out the testis samples, and then include all other organs from GTEX
  curclin <- dplyr::filter(subclin, TissueType == "Normal", X_primary_site != "Testis")
  ramat1 <- as.matrix(subrsem[row.names(currt), row.names(curclin)])
  # (2) Next, we do rowanno2, which makes boxplots for just the cancer expressions
  curclin <- dplyr::filter(subclin, TissueType == "Cancer")
  ramat2 <- as.matrix(subrsem[row.names(currt), row.names(curclin)])
  
  # Now, we need to sort rows based on the difference of cancer medians vs normals medians (sans testis).
  # This order is so that we are able to raise to the top those genes that are highly expressed in cancer but least
  # expressed in the normals on average (sans testis)
  CvsNratios <- matrixStats::rowMedians(ramat2) - matrixStats::rowMedians(ramat1)
  curordr <- sort(CvsNratios, decreasing = TRUE, index.return=TRUE)$ix
  currt <-  currt[curordr,  ]
  ramat1 <- ramat1[curordr, ]
  ramat2 <- ramat2[curordr, ]
  
  # Now, we prepare the combined right annotation module
  row_ha = rowAnnotation(NormalsLogExprs = anno_boxplot(ramat1, pch=20, size = unit(.5, "mm"), gp = gpar(fill = "#a9cc7a", lwd=0.5),
                                                        axis_param=list(gp=gpar(fontsize = 12)), width = unit(4, "cm")),
                         CancerLogExprs = anno_boxplot(ramat2, pch=20, size = unit(.5, "mm"), gp = gpar(fill = "#7aa2cc", lwd=0.5),
                                                       axis_param=list(gp=gpar(fontsize = 12)), width = unit(4, "cm")),
                         annotation_label = c("Normals LogExprs \n (sans Testis)", "Cancers LogExprs"),
                         annotation_name_rot=0,
                         gap = unit(3, "mm"))
  
  # spacing between annotations
  ht_opt$ROW_ANNO_PADDING = unit(3, "mm")
  
  hh <- ComplexHeatmap::Heatmap(currt, name = paste0("CTA expr vs TMB - TCGA GTEx (Xena) ",strcanType), col = col_fun, rect_gp = gpar(col= "white"),
                                heatmap_legend_param = list(title="Expr \n (Whole matrix \n Z-score)"), border=T, row_names_gp = gpar(fontsize = 14),
                                column_title = paste0("CTA expr vs TMB - TCGA GTEx (Xena)\n",strcanType),
                                column_title_gp = gpar(fontsize = 14), column_split = splitvec, column_gap = unit(3, "mm"),
                                row_title = "", cluster_columns = F, cluster_rows = F, show_row_dend = F, show_column_dend = F,
                                column_names_rot = 30, column_names_side=c("bottom", "top"), column_names_gp = gpar(fontsize = 12),
                                bottom_annotation = column_ha, right_annotation = row_ha,
                                top_annotation = HeatmapAnnotation(groups = anno_text(colnames(currt), rot=30, location=0, gp=gpar(fontsize = 12), just='left')))
  pdf(file = paste("./plots/HtMapTCGA_andNormals_CTAvsTMB_srtdByHigh_",gsub(" ","_",strcanType,fixed=TRUE),extrastrng,"_ZScrWholeMtrx.pdf",sep=""),
      useDingbats = FALSE, width=pltw, height=pltht)
  draw(hh, padding = unit(c(12, 28, 5, 7), "mm")) #bottom, left, top, right paddings
  dev.off()
  
}
