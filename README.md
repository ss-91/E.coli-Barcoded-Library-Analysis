# E.-coli-Barcoded-Library-Analysis
Basic pipeline for analyzing a barcoded screen

1. Prepare fastq files
	a. Download all relevant fastq files to the working directory
	b. Remove fastq files that are not needed

2. Converte raw reads to barcode counts (THIS CODE RUNS OVERNIGHT)
	* HEAVY computation - the code runs on parallel CPU and typically runs OVERNIGHT on the lab Server (avoid using it if only practicing)
	a. Edit the fastq2barcodeCounts.m to process relevant fastq files
	b. Edit the path to the ASKA_lookup_map.mat on dropbox (a hash pointing barcodes->genes)
	c. Run the script which saves the “dataset” variable (file name typically include the quality score cutoff, e.g., dataset_QS10).
	d. Manually move all fastq file into a subfolder called "fastq"

3. Prepare COUNTS files and plot basic statistics
	a. Edit the basicScreenStatistics.m script (e.g., order replicates and conditions as you wish)
	b. Manually copy/paste columns of the “COUNTS” variables into csv files:
		* Save all biological replicates of a single condition into a separate csv file
		* The first line should be titles ("ID" for gene names and then the replicate names)
		* Paste gene names from Matlab in the left most column
		* Paste the counts in the columns after.
		* Save as comma separated (Excel->save as->File Format: CSV UTF-8)
		* Name each file by its condition).
	c. Manually move all the csv file into a subfolder called "COUNTS" 
	
4. Plot basic statistics of the run 
	a. Edit the Counts_violinPlot.R script to upload the correct COUNTS files
	b. Run the script to /Users/mitchela2/Dropbox (UMass Medical School)/Mitchell Lab/Code repository/ASKA barcode screen/README.txtplot violin and bar plots.
	c. Optional: Manually save the violin plot from the R environment (save in SVG format for compatibility)

5. Identify hits and pathway enrichment
	Note - BarcodeScreen.R is a wrapper script that runs barcodeAnalysis.R on each condition (condition.csv vs control.csv). Enrichment is calculated by GO, Kegg and EcoCYC. While GO and Kegg are downloaded with every run. EcoCYC relies on a static version Amir downloaded on Dec 2019. The run required a local copy of "EcoCyc Pathways Entrez GMT.csv" in the working directory

	a. Edit the BarcodeScreen.R wrapper script. Save and "source" the changes.
	b. Run the script BarcodeScreen.R. Output files are saved to working directory as csv tables
	c. Manually move all the output csv files into a subfolder called "R OUTPUT" 
	d. Manually aggregate results from "hits" files (gene name and log enrichment) to a new excel table called "ScreenHits.xlsx". Save each condition in I separate sheet
	e. For each column of hits add gene description. 
 		- Upload the gene list to EcoCyc as a SmartTable.
		- Add a column to the table by "ADD PROPERTY COLUMN" -> Product. 
		- Export the new online table from EcoCyc: EXPORT -> to SpreadSheet File ... -> "common names"
		- Open the exported table in textedit, copy the "Product" to your hits excel file (double check the gene order and gene number matches).

	f. Manually aggregate results from "Enrichment" files to a new excel table called "PathwayEnrichmennt.xlsx" 

6. Compact GO results with REVIGO (require user interaction on-the-fly)
	a. Edit the plotPathwayEnrichment.m script (uses the excel file PathwayEnrichmennt.xlsx
	b. Run it. It will print to the screen a list of common GO terms.
	c. Copy this list and paste in the REVIGO website (http://revigo.irb.hr/)
	d. In REVIGO:
		- How large would you like the resulting list to be? Select "Medium"
		- Select a database with GO term sizes: Select "E. coli"
	e. Scroll down to see results (result page might look empty in the beginning)
	f. Export the results table "Export results to text table (CSV)"
	g. Sort the results so you can copy only the main categories.
	h. Past the main GO categories into a new Sheet on the PathwayEnrichmennt.xlsx files. Name the sheet REVIGO
	i. Rerun the plotPathwayEnrichment.m to get the "spot" plot for GO pathway enrichment.

6. Hierarchal Clustering of hits (Matlab)
	This script is highly dependent on the screen details and should be developed separately for each screen
