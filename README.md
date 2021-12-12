# Phylomeet: Community Phylogenetics
 
This is code I wrote for a (study)[https://onlinelibrary.wiley.com/doi/10.1002/ece3.3146] I led on how communities of arthropods differ in their evolutionary lineages, depending on which species of host tree they live on. 

## arcot_data
arcot is short for arthropod cottonwood data. The data for this study came from our coauthor Gina Wimp. She gathered the data for her dissertation in 2000-2003. 

### arcot_gw.xlsx 
the original spreadsheet provided by Gina.

### arcot_treenames.csv
contains the names assigned to the cottonwood tree hosts included in this study

### arcot.txt 
the data exported as a tab-delimited text file.

### v1
Exploring this project during Genes to Environments course at Northern Arizona University

### v2
Refining the approach from the class

## Figures
Contains figures output

## LMM
Results from linear mixed models of the results from the analysis.

## Results_ComDist
Output from analysis: community distance, a measure of the difference among communities in the analysis

## Results_NRI
Output from analysis: net relatedness index, a measure of how closely related individuals are within each community

## Results_PD
Output from analysis: phylogenetic distance, a measure of how phylogenetically close a community is

## Trees
Phylogenetic tree framework constructed in Mesquite and exported as a Nexus file, then converted to newick/phylip format. Random variation added to branch lengths here (done in dataprep)

## arcot_analysis.R
Script that conducts main phylogenetic analyses for the project. Pacakge picante provides most functions to calculate the metrics. Results files exported to Results folders.

## arcot_dataprep.R
Script that pulls in raw community data and reformats it in preparation for analysis. It also pulls in phylogenetic tree files, adds random variation 10x, and writes these variations to use in analysis.

## figures_Manuscript.R
Code for producing figures in Figures folder

## phylos_Misof_res.R
Code for adding variation to branch length in the phylogenies

## phynames
dictionary of full taxon names and the codes used for them


