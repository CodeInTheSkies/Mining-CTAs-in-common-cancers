# This script reads the BIG Raw Data file (7GB) of TOIL reCompute of
# GTEx-TARGET-TCGA data, downloaded from:
# https://xenabrowser.net/datapages/?dataset=TcgaTargetGtex_RSEM_Hugo_norm_count&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
# "rsem_Hugo_norm_count" file
#
# After reading the raw file, the script then subsets the data to just the
# "genes of interest" list. It also provides provision to add any genes outside
# that list as necessary. We also exclude all the TARGET data (which are
# pediatric patient data that are outside the scope of this analysis). We then
# save this subset data as an Rds file. This way, we can simply load the smaller
# subset Rds easily and proceed with further processing in other scripts
#
# NOTE: If in the future, new genes need to be added to the full analysis, then
# first this script needs to be modified to include those genes and then run, so
# that the stored subset matrix in the Rds file gets updated first. The new and
# updated Rds will contain the new genes, which can then be used by other
# scripts for further processing. In other words, the addition of new genes is a
# 2-step process.
#
# NOTE also that the raw files in bigdata are NOT included (see .gitignore) as
# part of the folder or repo due to their large sizes. Before running this
# script, the raw files need to be downloaded from the direct links specified in
# the "Readme".
# 
# However, the CTA gene list (in smalldata) are included in the repo for convenience

library(rstudioapi)

# The following line automatically gets the current script's directory, and sets
# it to the current working directory. However, this may work only within
# RStudio's editor. Once the working directory is successfully set using this
# line, the relative paths below will work
cwd <- dirname(rstudioapi::getActiveDocumentContext()$path)
cat("The current directory of the source script is: ", cwd)
setwd(cwd)

# Reading the big raw data matrix; this line may take some time to run!
rawmat <- read.delim("./data/TcgaTargetGtex_RSEM_Hugo_norm_count")

# Reading the list of "genes of interest"
glist <- read.table("./smalldata/CompleteCTAList.txt", quote="", comment.char="")
# convert into a vector, from a single-column dataframe
glist <- glist$V1
# Here, if you like to add other genes, or remove any genes already in the list,
# simply use the lines below
# Removing any gene not needed
# Specify here, any gene to be removed
rmvGenesList <- c()  
glist <- glist[!(glist %in% rmvGenesList)]
# Add any genes that are needed in addition to what is on the "CompleteCTAList.txt" list
# Specify here, any gene to be added
addGenesList <- c("PLK4",
                  "CCNE1",
                  "CDK2",
                  "CCND1",
                  "CDK4",
                  "CDK6",
                  "CDK1",
                  "CCNA1",
                  "CCNA2",
                  "CCNB1",
                  "CDKN2A")
glist <- c(glist, addGenesList)
# Changing all in glist to upperCase for consistency
glist=unique(toupper(glist))
spath="./plots/"

# Checking for genes in list to be present in data
ntfnd <- glist[!(glist %in% rawmat$sample)]
cat("\nThese", length(ntfnd), "genes from the list of", length(glist), "'genes of interest' are NOT found in the data:\n", ntfnd, "\n")

# Subsetting to ONLY the "genes of interest" (provided the gene occurs in the data)
rawmat <- rawmat[(rawmat$sample %in% glist) , ]

# Setting gene names to be row names, and then removing the "sample" column
row.names(rawmat) <- rawmat$sample
rawmat <- rawmat[, !(colnames(rawmat) %in% "sample")]

# Excluding TARGET data
cnames <- sapply(strsplit(colnames(rawmat), split='.', fixed=TRUE), function(x) (x[1]))
rawmat <- rawmat[, (cnames == "GTEX" | cnames == "TCGA") ]

# At this point, the rawmat matrix is clean, has ONLY the genes or rows that are of interest, and then also has ONLY the columns of GTEx and TCGA.
# So, we simply store this matrix as Rds
saveRDS(rawmat, file="./data/subset_GTEx_TCGA_expressions.Rds")

