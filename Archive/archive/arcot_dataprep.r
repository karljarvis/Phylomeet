#############################################################
### Community Phylogenetic Analyses  
###	Arthropod communities on cottonwood hosts
###	Dataset from Wimp et al. 2004
###
### Code by Karl Jarvis
#############################################################

# load picante package
# install.packages('picante', dependencies = TRUE, repos = 'http://R.research.att.com/')
# install.packages('geiger', dependencies = TRUE, repos = 'http://R.research.att.com/')
library(picante)


# Picante manual
# vignette('picante-intro')

# Set directory for input data: community file and phylogeny file
dataDir = '~/Projects/Phylomeet/Analysis/arcot_data/'
treeDir = paste0(dataDir,'trees')

#############################################################
# Load in communities, take subsets and format for analysis #
#############################################################

# load community data
fullCom = read.table(paste(dataDir, 'arcot.txt', sep=''))
crosstype = rep(rep(c('fr','fo','na','na'), each=10), times=4)
year = c(rep(2000:2003, each=40))
comdf = cbind(year, crosstype, fullCom)

# factors 
comboNames = c('fr2000','fo2000','na2000','fr2001','fo2001','na2001','fr2002','fo2002','na2002','fr2003','fo2003','na2003')
spaces = c(0.2,0.2,0.2,1,0.2,0.2,1,0.2,0.2,1,0.2,0.2)

###############################################
# pool data
	# byCT: sum of occcurrences in all years by crosstype
	frSum = apply(comdf[comdf$crosstype=='fr',3:ncol(comdf)],2,sum)
	foSum = apply(comdf[comdf$crosstype=='fo',3:ncol(comdf)],2,sum)
	naSum = apply(comdf[comdf$crosstype=='na',3:ncol(comdf)],2,sum)
	byCT = rbind(frSum,foSum,naSum)
	
	# byYear: sum of occcurrences by year on all crosstypes
	sum2000 = apply(comdf[comdf$year==2000,3:ncol(comdf)],2,sum)
	sum2001 = apply(comdf[comdf$year==2001,3:ncol(comdf)],2,sum)
	sum2002 = apply(comdf[comdf$year==2002,3:ncol(comdf)],2,sum)
	sum2003 = apply(comdf[comdf$year==2003,3:ncol(comdf)],2,sum)
	byYear = rbind(sum2000,sum2001,sum2002,sum2003)
	
	# byCtYr: sum of occurrences by year and by crosstype
	byCtYrMaker = function(com) 
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
		byCtYr = rbind(fr2000,fo2000,na2000,fr2001,fo2001,na2001,fr2002,fo2002,na2002,fr2003,fo2003,na2003)
		byCtYr
	}
	
	byCtMaker = function(com)
	{
		fr = apply(com[com$crosstype=='fr',3:ncol(com)],2,sum)
		fo = apply(com[com$crosstype=='fo',3:ncol(com)],2,sum)
		na = apply(com[com$crosstype=='na',3:ncol(com)],2,sum)
		byH = rbind(fr, fo, na)
		byH		
	}
	
#####################################################
# Subset the data
	byCtYr = byCtYrMaker(com=comdf)

#####################################################
	# Full community
	# sumAll: sum of all occurrences in all years on all crosstypes
	sumAll = apply(comdf[,4:ncol(comdf)],2,sum)
	
	# CtYr: pooled by crosstype and Year
	CtYr = byCtYrMaker(com=comdf)
	
	# Ct: pooled by crosstype only
	Ct = byCtMaker(com=comdf)	

#####################################################
# 	# Common taxa only
# 	fullsumAll = rbind(fullCom,sumAll)
# 
# 	# sp2: exclude species that have only one occurrence
# 	sp2 = fullsumAll[1:160,sumAll >= 2] 
# 
# 	# sp2CtYr: sp2, pooled by crosstype and year
# 	sp2CtYr = byCtYrMaker(com=cbind(year, crosstype, sp2))
# 
# 	# sp3: exclude species that have only one occurrence
# 	sp3 = fullsumAll[1:160,sumAll >= 3] 
# 
# 	# sp3CtYr: com3sp, pooled by crosstype and year
# 	sp3CtYr = byCtYrMaker(com=cbind(year, crosstype, sp3))
# 
# 	# yr2: all except species that are only present in one year
# 	y = byYear
# 	for(i in 1:ncol(y)) 
# 	{ 
# 		for(j in 1:nrow(y)) 
# 		{ 
# 			if (byYear[j,i] > 0) {y[j,i] = 1} 
# 		}	
# 	}
# 	nyear = apply(y,2,sum)
# 	fullnyear = rbind(fullCom,nyear)
# 	yr2 = fullnyear[1:160,nyear >= 2]
# 
# 	# yr2CtYr: y2, pooled by crosstype and year
# 	yr2CtYr = byCtYrMaker(com=cbind(year, crosstype, yr2))
# 		
# 	# yr2sp2: intersection of sp2 and yr2
# 	sp2year = rbind(sp2,nyear)
# 	yr2sp2 = sp2year[1:160,sp2year[161,] >= 2]
# 
# 	# yr2sp2CtYr: yr2sp2, pooled by crosstype and year
# 	yr2sp2CtYr = byCtYrMaker(com=cbind(year, crosstype, yr2sp2))
# 
# 	# yr2sp3: intersection of sp3 and yr2
# 	sp3year = rbind(sp3,nyear)
# 	yr2sp3 = sp3year[1:160,sp3year[161,] >= 2]
# 
# 	# yr2sp3CtYr: yr2sp3, pooled by crosstype and year
# 	yr2sp3CtYr = byCtYrMaker(com=cbind(year, crosstype, yr2sp3))
# 	
# 	# yr2sp3H: yr2sp3, pooled by crosstype
# 	yr2sp3Ct = byCtMaker(com=cbind(year, crosstype, yr2sp3))
# 	

#####################################################
# Select which community subsets to use

# list of full community and subsets	
	 # comList = lapply(list(fullCom, CtYr, sp2, sp3, yr2, yr2sp2, yr2sp3, sp2CtYr, sp3CtYr, yr2CtYr, yr2sp2CtYr, yr2sp3CtYr), as.matrix)
	 # names(comList) = c('fullCom', 'CtYr', 'sp2', 'sp3', 'yr2', 'yr2sp2', 'yr2sp3', 'sp2CtYr', 'sp3CtYr', 'yr2CtYr', 'yr2sp2CtYr', 'yr2sp3CtYr')
	 
# list of communities to use in this analysis
# 	 comList = lapply(list(fullCom, yr2sp3, CtYr, yr2sp3CtYr, Ct, yr2sp3Ct), as.matrix)
# 	 names(comList) = c('full_indiv', 'common_indiv', 'full_pooled', 'common_pooled', 'full_crosstype', 'common_crosstype')
	 
# list of communities to use in this analysis
  comList = lapply(list(fullCom, CtYr, Ct), as.matrix)
  names(comList) = c('full_indiv', 'full_pooled', 'full_crosstype')


#####################################################
# Load in phylogenies, prune, and create 
# distance matrices for analysis 
#####################################################
# Read in phylogenies
phyleNames = list.files(dataDir, pattern='[.]tre')
phyNames = gsub('.tre','',phyleNames)
phyList = vector("list", length(phyNames))
names(phyList) = gsub('.tre','',phyleNames)
for(i in 1:length(phyleNames))
{
  phyList[[i]] = read.tree(paste0(dataDir, '/trees/', phyleNames[i]))
}

###############################################
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
	phyDist = vector('list',length=length(comList))
	names(phyDist) = names(comList)
	for (i in 1:length(comList)) 
	{	
		phyDist[[i]] = vector('list', length=length(phyNames))
		names(phyDist[[i]]) = phyNames
		for (j in 1:length(phyNames)) 
		{ 	
			phyDist[[i]][[j]] = cophenetic(phyPrune[[i]][[j]]) 	
		}
	}

###############################################
# Plotting

# Plot phylogenies
	# par(mfrow=c(4,2))
	# for(i in 1:length(phyNames))
	# {
		# plot(phyList[[i]], cex=0.5, show.tip.label=F)
	# }



	# plot(phy$full_pooled$ultra, cex=0.6)		

# Plotting phylogenies
	# First five sets of phylogenies to get overview
	# par(mfrow=c(4,3))
	# for (c in 1:9) 
	# {	
		# plot(phy[[c]][[1]], show.tip.label=F) 
	# }

	# # large versions of full phylogenies
	# par(mfrow=c(1,2))
	# plot(phy$fullCom$ultra, cex=0.5)
	# plot(phy$fullCom$equal, cex=0.5)

		