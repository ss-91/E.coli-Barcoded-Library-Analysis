## this script analyzise the reuslts of a barcode screen (single condition)

## Load libraries (suppresss warnings)
suppressPackageStartupMessages({
  library("dplyr") # Part of the tidyverse module that has some functionalities absent in the main package
  library("tidyverse") # Data wrangling package (selecting, filtering, subsetting dataframe)
  library("reshape2") # Similar to melt, used to reshape data to fit a certain structure
  library("pheatmap") # For plotting heatmaps, can generate correlation matrix from raw counts and develop heatmap
  library("DESeq2") # library for differential expression analysis
  library("ggpubr") # For plotting, just like ggplot. Not used anyway in the pipeline
  library("org.EcK12.eg.db") # Updated database of E.coli K-12 librarry from Bioconductor
  library("GO.db") # Gene Ontology (GO) database 
  library("gage") # Library for functional enrichment of gene sets
  library("AnnotationDbi") # ?
  library("BiocGenerics") # ?
  library("DEBRA") # testing if DESeq2 statistics can be applied to barcode data
  library("pathview") # a library for plotting KEGG pathways
})

# load the relevant functions (after switching to the correct working directory)
## this is the main script for analyzing screens with the ASKA barcoded library
# Note - before running, the the barcode counts should be prepared with a matlab script

## USAGE
# specify the names of two csv files with counts (one for exp replicate and one for 
# control replicates) and string that will be used for naming the output files
strDir = "/Users/serkansayin/Dropbox (UMass Medical School)/Bacteria_Host (gemcitabine)/Raw_results/ASKA barcode screen"; setwd(strDir);
source("barcodeAnalysis.R") # load the barcodeAnalysis function to the memory
source("plotKEGG.R") # load the plotKEGG function to the memory

# input for a single run (biological condition)
#strCountsExpFile = "./COUNTS/COUNTS_LB5FU.csv" # name of CSV file with all barcode counts for exp (single condition)
#strCountsContFile = "./COUNTS/COUNTS_LBND.csv" # name of CSV file with all barcode counts for control (single condition)
#strTitle = "LB-5FU" # title for all output file names
#barcodeAnalysis(strCountsExpFile, strCountsContFile, "LB-5FU");

# OR input for multiple runs (run each biological condition)
barcodeAnalysis("./COUNTS/Counts_Gemcitabine.csv", "./COUNTS/Counts_NoDrug.csv", "NoDrug-Gemcitabine");
#barcodeAnalysis("./COUNTS/COUNTS_LBFUDR.csv", "./COUNTS/COUNTS_LBND.csv", "LB-FUDR");
#barcodeAnalysis("./COUNTS/COUNTS_M9FUDR.csv", "./COUNTS/COUNTS_M9ND.csv", "M9-FUDR");
#barcodeAnalysis("./COUNTS/COUNTS_M95FU.csv","./COUNTS/COUNTS_M9ND.csv", "M9-5FU");
#barcodeAnalysis("./COUNTS/COUNTS_LBND.csv", "./COUNTS/COUNTS_M9ND.csv", "LB-M9");

print("Done")

