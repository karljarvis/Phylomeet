##########
# Perform community phylogenetic analyses of arthropod communities on cottonwood hosts
# Karl Jarvis June 13, 2014
##########

source("/Users/kjj/Projects/Phylomeet/Analysis/arcot_dataprep_v1.R")
# source("/scratch/kj375/Phylomeet/arcot_dataprep_v1.R")

##########
# Faith's Phylogenetic Diversity (PD)
##########

PD = vector("list", length(comList))
names(PD) = names(comList)
for (i in 1:length(comList)) 
{	
  PD[[i]] = vector("list", length(phyNames))
	for (j in 1:length(phyList)) 
	{ 	
		PD[[i]][[j]] = pd(comList[[i]], phyList[[j]], include.root=TRUE)
		write.csv(PD[[i]][[j]], paste0(PDdir, "PD_", names(comList)[i], "_", phyNames[j], ".csv"))	
    print(j)
	}
}


##########
# Net Relatedness Index (NRI) 
##########

NRI_abund = NRI_pres = vector("list", length(comList))
for (i in 1:length(comList)) 
{	
  NRI_abund[[i]] = NRI_pres[[i]] = vector("list", length(phyNames))
	for (j in 1:length(phyList)) 
	{ 	
		NRI_abund[[i]][[j]] = ses.mpd(comList[[i]], phyDist[[i]][[j]], null.model="richness", abundance.weighted=T, runs = 999)
		NRI_pres[[i]][[j]] = ses.mpd(comList[[i]], phyDist[[i]][[j]], null.model="richness", abundance.weighted=F, runs = 999)
		write.csv(NRI_abund[[i]][[j]], paste0(NRIdir, "NRI_abund_", names(comList)[i],"_", phyNames[j], ".csv"))
		write.csv(NRI_pres[[i]][[j]], paste0(NRIdir, "NRI_pres_", names(comList)[i],"_", phyNames[j], ".csv"))
    print(j)
	}
}


##########
# Phylogenetic Beta Diversity: Community Distance (ComDist)
##########

comDist_abund = comDist_pres = vector("list", length(comList))
names(comDist_abund) = names(comDist_pres) = names(comList)
for (i in 1:length(comList)) 
{	
	comDist_abund[[i]] = comDist_pres[[i]] = vector("list", length(phyNames)) 
	names(comDist_abund[[i]]) = names(comDist_pres[[i]]) = phyNames
	for (j in 1:length(phyList)) 
	{ 	
		comDist_abund[[i]][[j]] = comdist(comList[[i]], phyDist[[i]][[j]], abundance.weighted=T)
		comDist_pres[[i]][[j]] = comdist(comList[[i]], phyDist[[i]][[j]], abundance.weighted=F)
		write.csv(as.matrix(comDist_abund[[i]][[j]]), paste0(comDistDir, "comDist_abund_", names(comList)[i], "_", phyNames[j], ".csv"))
		write.csv(as.matrix(comDist_pres[[i]][[j]]), paste0(comDistDir, "comDist_pres_", names(comList)[i], "_", phyNames[j], ".csv"))
    print(j)
	}
}

