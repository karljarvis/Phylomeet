##########
# Import community arthropod data and phylogenies for community phylogenetic analyses
# Karl Jarvis June 13, 2014
##########


# load picante package
library(picante)

dir = "~/Projects/Phylomeet/Analysis/"
# dir = "/scratch/kj375/Phylomeet/"

# Set directory for input data: community file and phylogeny file
dataDir = paste0(dir, "arcot_data/")

# Folder for phylogenies
treeDir = paste0(dir, "Trees/")

# Set folder for output of PD results
PDdir = paste0(dir, "Results_PD/")

# Set folder for output of NRI results
NRIdir = paste0(dir, "Results_NRI/")

# Set folder for output of ComDist results
comDistDir = paste0(dir, "Results_ComDist/")

# Set directory for figures
figDir <- paste0(dir, "Figures/")


##########
# Load, pool, and organize community data into data structures
##########

# load community data
crosstypes = c("fr","fo","na")
years = factor(c(2000:2003))
types = rep(crosstypes, 4)
yrs = rep(years, each=3)

fullCom = read.table(paste(dataDir, 'arcot.txt', sep=''))
crosstype = rep(rep(c('fr','fo','na','na'), each=10), times=4)
year = rep(2000:2003, each=40)
comdf = cbind(year, crosstype, fullCom)

# pool occurrences by year and by crosstype
CtYrMaker = function(com) 
{ 	
	fr2000 = apply(com[com$year==2000 & com$crosstype=='fr',3:ncol(com)],2,sum)
	fo2000 = apply(com[com$year==2000 & com$crosstype=='fo',3:ncol(com)],2,sum)
	na2000 = apply(com[com$year==2000 & com$crosstype=='na',3:ncol(com)],2,sum)
	fr2001 = apply(com[com$year==2001 & com$crosstype=='fr',3:ncol(com)],2,sum)
	fo2001 = apply(com[com$year==2001 & com$crosstype=='fo',3:ncol(com)],2,sum)
	na2001 = apply(com[com$year==2001 & com$crosstype=='na',3:ncol(com)],2,sum)
	fr2002 = apply(com[com$year==2002 & com$crosstype=='fr',3:ncol(com)],2,sum)
	fo2002 = apply(com[com$year==2002 & com$crosstype=='fo',3:ncol(com)],2,sum)
	na2002 = apply(com[com$year==2002 & com$crosstype=='na',3:ncol(com)],2,sum)
	fr2003 = apply(com[com$year==2003 & com$crosstype=='fr',3:ncol(com)],2,sum)
	fo2003 = apply(com[com$year==2003 & com$crosstype=='fo',3:ncol(com)],2,sum)
	na2003 = apply(com[com$year==2003 & com$crosstype=='na',3:ncol(com)],2,sum)
	out = rbind(fr2000,fo2000,na2000,fr2001,fo2001,na2001,fr2002,fo2002,na2002,fr2003,fo2003,na2003)
	out
}
CtYr = CtYrMaker(com=comdf)

# pool occurrences by crosstype
CtMaker = function(com)
{
	fr = apply(com[com$crosstype=='fr',3:ncol(com)],2,sum)
	fo = apply(com[com$crosstype=='fo',3:ncol(com)],2,sum)
	na = apply(com[com$crosstype=='na',3:ncol(com)],2,sum)
	out = rbind(fr, fo, na)
	out		
}
Ct = CtMaker(com=comdf)	
	
# list of communities to use in this analysis
comList = lapply(list(fullCom, CtYr, Ct), as.matrix)
names(comList) = c('indiv', 'pooled', 'crosstype')
CtYrNames <- c("fr2000","fo2000","na2000","fr2001","fo2001","na2001","fr2002","fo2002","na2002","fr2003","fo2003","na2003")
spaces <- c(0.2,0.2,0.2,1,0.2,0.2,1,0.2,0.2,1,0.2,0.2)

##########
# Load in phylogenies, prune, and create distance matrices for analysis 
##########

# Read in phylogenies
phyleNames = list.files(treeDir, pattern='[.]tre')
phyNames = gsub('.tre','',phyleNames)
phyList = vector("list", length(phyNames))
names(phyList) = gsub('.tre','',phyleNames)
for(i in 1:length(phyleNames))
{
  phyList[[i]] = read.tree(paste0(treeDir, phyleNames[i]))
}

# Prune phylogenies by community
phyPrune = vector('list',length=length(comList))
names(phyPrune) = names(comList)
for (i in 1:length(comList)) 
{	
	phyPrune[[i]] = vector('list', length=length(phyNames))
	names(phyPrune[[i]]) = phyNames
	for (j in 1:length(phyNames)) 
	{ 	
		phyPrune[[i]][[j]] = prune.sample(comList[[i]], phyList[[j]]) 
	}
}

# Create distance matrix for each topology
phyDist = vector('list', length(comList))
names(phyDist) = names(comList)
for (i in 1:length(comList)) 
{	
	phyDist[[i]] = vector('list', length(phyNames))
	names(phyDist[[i]]) = phyNames
	for (j in 1:length(phyNames)) 
	{ 	
		phyDist[[i]][[j]] = cophenetic(phyPrune[[i]][[j]]) 	
	}
}
