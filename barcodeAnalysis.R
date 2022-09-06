barcodeAnalysis <- function(strCountsExpFile, strCountsContFile, strTitle){
  
  #####################
  ### INITIALIZING  ###
  #####################
  
  
  ## User definitions
  #strDir = "/Users/mitchela2/Downloads/R analysis/"; setwd(strDir);
  #strCountsExpFile = "COUNTS_LB5FU.csv" # name of CSV file with all barcode counts for exp (single condition)
  #strCountsContFile = "COUNTS_LBND.csv" # name of CSV file with all barcode counts for control (single condition)
  #strTitle = "LB-5FU" # title for all output file names
  
  # cutoffs for single barcode statstics
  log2cutoff = log2(4); # cutoff for enrichment in log2 (e.g., 2 is 4-fold)
  fdrBarcodeCutoff = 0.05 # cutoff for FDR adjusted p-value for over/under representation of a barcode
  minCountPerGene = 10 # minimal number of count for each gene to be included in analysis (sum across all conditions)
  
  # cutoff for Gene set enrichement (GO, Kegg)
  fdrGeneSetCutoff = 0.1 #0.25 # cutoff for FDR adjusted p-value for over/under representation of a barcode
  
  ## Load dataset (counts) and prepare for further analysis
  counts_exp <- as.matrix(read.table(strCountsExpFile, header = T, row.names = 1, sep=",")) # load the counts
  counts_cont <- as.matrix(read.table(strCountsContFile, header = T, row.names = 1, sep=",")) # load the counts
  
  total <- cbind(counts_exp,counts_cont) # bind the two dataframes
  exp.idx <- seq(1,ncol(counts_exp)) # column indexes of "experiment" (e.g., with drug) 
  cont.idx <- seq(ncol(counts_exp)+1,ncol(counts_exp)+ncol(counts_cont)) # column indexes of controls (e.g., w/o drug) 
  
  # extract replicate names
  allColNames = colnames(counts_exp)
  replicatesNamesExp = allColNames[seq(1,ncol(counts_exp))]
  
  allColNames = colnames(counts_cont)
  replicatesNamesCont = allColNames[seq(1,ncol(counts_cont))]
  
  ########################
  ### Find screen hits ###
  ########################
  
  ## Run DEBRA (a DESeq2 wrapper)
  drb=DEBRA(total, replicatesNamesCont, replicatesNamesExp, method="DESeq2(Wald)")
  results=resultsDRB(drb); # make a result table from the DEBRA output
  results <- na.omit(results) # remove missing values
  # need to correct pval and padjust if they are too small
  results <- transform(results, pvalue = ifelse(pvalue == 0, 6.063463e-172, pvalue))
  results <- transform(results, padj = ifelse(padj == 0, 6.900220e-169, padj))
  
  # logical vector for significant hits (by pval and log-fold) 
  tf <- ((abs(results$log2FoldChange) > log2cutoff) & (results$padj < fdrBarcodeCutoff)) # over or under represented 
  
  ## export hit to csv file
  hits <- subset(results,tf)
  hits <- hits[order(hits$log2FoldChange),]
  write.csv(hits, file = paste(strTitle,".csv")) # save list of significant barcodes
  
  ## Scatter plots of hits
  pdf(paste(strTitle,".pdf")); 
  
  # prepare x and y vectors
  x = results$log2FoldChange; 
  y= -log10(results$padj);
  
  #cols <- densCols(results$log2FoldChange, -log10(results$pvalue)) # calculate colors for point the scatter plot by local density 
  
  #plot(x, y, col="gray", panel.first = grid(), xlim=c(-abs(max(x)),abs(max(x))),pch = 20, cex = 0.2)
  
  plot(x, y, col="gray", panel.first = grid(), xlim=c(-9, 9),ylim=c(0, 170),pch = 20, cex = 1)
  title(main = strTitle, xlab = "log2(fold-change)", ylab = "-log10(adjusted p-value)")
  points(x[tf], y[tf],pch = 20, cex = 1, col="blue")
  abline(v=0, col="gray")
  abline(v=c(-log2cutoff,log2cutoff), col = "blue")
  abline(h=-log10(fdrBarcodeCutoff), col = "blue")
  
  #text(x[tf],y[tf],lab = rownames(results)[tf], cex = 0.5, col = "red", adj=c(0,0))
  
  dev.off() # close the pdf file

# write the log2 and pval into a file
  strFile = paste(strTitle,"- log2 and pval.csv");
  write.table(results, strFile,sep = ",");
  
  
  ###########################
  ### Gene set enrichment ###
  ###########################
  
  cnts <- total # make a copy of the counts matrix
  sel.rn=rowSums(cnts) > minCountPerGene # remove rows with low counts from the analysis
  cnts=cnts[sel.rn,]
  
  cat(nrow(total)-nrow(cnts),"barcodesremoved by filter low counts, total below" ,minCountPerGene, "reads")
  
  # make a normalized version of counts (stabilizes variance?) follow gage documentation
  libsizes=colSums(cnts)
  size.factor=libsizes/exp(mean(log(libsizes)))
  cnts.norm=t(t(cnts)/size.factor)
  range(cnts.norm)
  
  cnts.norm=log2(cnts.norm+8)
  range(cnts.norm)
  
  # change row names to Entrez ids
  symbol2entrez <- mapIds(org.EcK12.eg.db, rownames(cnts), "ENTREZID", keytype="ALIAS")
  symbol2entrez.df <- as.data.frame(symbol2entrez) # make a data frame version of the map
  
  # multiple steps to replace row names by the matching entrez ids
  onlyEntrez <- data.frame(gene = symbol2entrez.df$symbol2entrez)
  cnts.norm.2 <- cbind(onlyEntrez, cnts.norm)
  cnts.norm.2 <- cnts.norm.2[!is.na(cnts.norm.2$gene),]
  cnts.norm.3 <- cnts.norm.2[!duplicated(cnts.norm.2[,c('gene')]),]
  rownames(cnts.norm.3) <- cnts.norm.3$gene # replace row names by gene name
  cnts.norm.4 <- cnts.norm.3[,-1] # remove the gene column 
  
  # Generate kegg database for E. coli (using the gage library)
  #kg.eco=kegg.gsets("eco") # bnames
  kg.eco.eg=kegg.gsets("eco", id.type="entrez") # entrez gene
  
  # run GAGE (compare cont/exp groups, use KS for statistics and DONT allow changes in diffent directions)
  kg.data <- gage(cnts.norm.4, gsets = kg.eco.eg$kg.sets, ref = cont.idx,samp = exp.idx, compare ="as.group", saaTest=gs.KSTest, same.dir=TRUE)
  keggGreater <- as.data.frame(kg.data$greater)
  keggGreater <- keggGreater[keggGreater$q.val < fdrGeneSetCutoff & !is.na(keggGreater$q.val),]
  #rownames(keggGreater)
  keggLess <- as.data.frame(kg.data$less)
  keggLess <- keggLess[keggLess$q.val < fdrGeneSetCutoff & !is.na(keggLess$q.val),]

  # Overlay log2 data on enriched/depleted KEGG pathways
   keggAll <- rbind(keggLess,keggGreater)
   keggNames <- rownames(keggAll)
   for(i in 1:length(keggNames)){
     curKegg = keggNames[i];
     if(keggAll$set.size[i]<200){
       tokens <- str_split(curKegg, " ", 2, simplify = TRUE)
       plotKEGG(results, tokens[1], "eco", strTitle)
     }
   }
  
  # GO enrichment
  go=go.gsets("E coli strain K12", id.type="entrez") # entrez gene
  go.bp=go$go.sets[go$go.subs$BP] # subset of GO (Bio. Process)
  go.mf=go$go.sets[go$go.subs$MF] # subset of GO (Mol. Function)
  go.cc=go$go.sets[go$go.subs$CC] # subset of GO (Cellular componenet)
  
  # run GAGE (compare cont/exp groups, use KS for statistics and DONT allow changes in diffent directions)
  go.data <- gage(cnts.norm.4, gsets = go.bp, ref = cont.idx,samp = exp.idx, compare ="as.group", saaTest=gs.KSTest, same.dir=TRUE)
  attributes(go.data)
  goGreater <- as.data.frame(go.data$greater)
  goGreater <- goGreater[goGreater$q.val < fdrGeneSetCutoff & !is.na(goGreater$q.val),]
  rownames(goGreater)
  goLess <- as.data.frame(go.data$less)
  goLess <- goLess[goLess$q.val < fdrGeneSetCutoff & !is.na(goLess$q.val),]
  
  # EcoCyc enrichment
  ecocyc=readList("EcoCyc Pathways Entrez GMT.csv")
  
  # run GAGE (compare cont/exp groups, use KS for statistics and DONT allow changes in diffent directions)
  eco.data <- gage(cnts.norm.4, gsets = ecocyc, ref = cont.idx,samp = exp.idx, compare ="as.group", saaTest=gs.KSTest, same.dir=TRUE)
  attributes(eco.data)
  ecoGreater <- as.data.frame(eco.data$greater)
  ecoGreater <- ecoGreater[ecoGreater$q.val < fdrGeneSetCutoff & !is.na(ecoGreater$q.val),]
  rownames(ecoGreater)
  ecoLess <- as.data.frame(eco.data$less)
  ecoLess <- ecoLess[ecoLess$q.val < fdrGeneSetCutoff & !is.na(ecoLess$q.val),]
  
  
  ## Save the enrichment and depletion results to a CSV file
  strFile = paste(strTitle,"- Enrichment.csv");
  strFileDep = paste(strTitle,"- Depletion.csv");
  
  if(nrow(keggGreater)) {
    write.table(keggGreater, strFile, sep = ",")
    write.table(keggLess, strFileDep, sep = ",")
  }
  if(nrow(ecoGreater)) {
    write.table(ecoGreater, strFile, sep = ",", col.names = !file.exists(strFile), append = T)
    write.table(ecoLess, strFileDep, sep = ",", col.names = !file.exists(strFile), append = T)
  }
  if(nrow(goGreater)) {
    write.table(goGreater, strFile, sep = ",", col.names = !file.exists(strFile), append = T)
    write.table(goLess, strFileDep, sep = ",", col.names = !file.exists(strFile), append = T)
  }
  
  return(1)
}