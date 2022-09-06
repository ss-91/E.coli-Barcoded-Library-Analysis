plotKEGG <- function(resultsDE, strPathwayID, strSpecies, strFileName){
  ## plot a specific KEGG pathway overlay an intensity map
  # INPUT: 
  # resultsDE - a dataframe holding the DE genes and the enrichment values
  # strPathwayID - KEGG pathway to plot (e.g., "00290")
  # strSpecies - species ID ("eco" - E coli, "hsa" - human)
  # strFileName - a string to be inserted to the file names
  
  ## Function definitions
  colorLimits = 2;
  
  
  ## Covert to Entrez IDs, remove extra fields from DF 
  symbol2entrez <- mapIds(org.EcK12.eg.db, rownames(resultsDE), "ENTREZID", keytype="ALIAS")
  symbol2entrez.df <- as.data.frame(symbol2entrez) # make a data frame version of the map
  
  onlyEntrez <- data.frame(gene = symbol2entrez.df$symbol2entrez)
  resultsDE.2 <- cbind(onlyEntrez, resultsDE$log2FoldChange) # bind the Entrez ids and log2FoldChange
  resultsDE.2 <- resultsDE.2[!is.na(resultsDE.2$gene),] # remove NA
  resultsDE.3 <- resultsDE.2[!duplicated(resultsDE.2[,c('gene')]),] # remove duplicates (two genes with the same Entrez Id)
  rownames(resultsDE.3) <- resultsDE.3$gene # replace row names by gene name (Entrez id)
  colnames(resultsDE.3) <- c("gene","log2FoldChange") # change column names to something sensible
  keeps <- "log2FoldChange"
  allChange <- resultsDE.3[ , keeps, drop = FALSE] # weird syntax required to keep "cnts.norm.4" as a dataframe 
  
  #genes<-row.names(allChange)[which(allChange$log2FoldChange>=2)]
  
  ## plot the pathview
  pv.out <- pathview(gene.data = allChange, pathway.id = strPathwayID,species = strSpecies, out.suffix = strFileName, 
                     kegg.native = T, limit=list(colorLimits,colorLimits),same.layer = F, na.col="white",low="blue",mid="gray",high="red");
  
  # set "kegg.native = T" to get a pdf that can be view in Graphviz
  # set "same.layer = T" to put all graphics on same layer (slows the script)
} 
