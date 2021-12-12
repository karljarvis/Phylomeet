### This code is to simulate hybrid communities to compare to actual hybrid communities. We want to see if actual hybrid communities are a mix of species that you might expect if you were to throw two genomes together, or whether there are unique, emergent properties to the communities found on hybrids that are fundamentally different than those of simulated communities that blend the characteristics of the pure hosts.

###################################################################
# Load in communities: Fremont, narrowleaf, and F1 for each year
###################################################################

	setwd("~/Documents/Phylomeet/Analysis/Simulated_Hybrids")
	data <- read.csv("~/Documents/Phylomeet/Analysis/Simulated_Hybrids/arcot_gw.csv", row.names=1)
	type <- rep(c(rep("Fremont",10),rep("F1",10),rep("Narrowleaf",20)),4)
	year <- c(rep("2000",40),rep("2001",40),rep("2002",40),rep("2003",40))
	host <- rep(c(1:10,1:10,1:20),4)
	com <- cbind(type,year,host,data)

# Find distribution of abundance for each species in each community
	for(j in 0:4)
	{
		dev.new()
		par(mfrow=c(4,4))
		for(i in 1:16) { hist(com[com$type == "Fremont",j*16+3+i]) }
	}
		dev.new()
		par(mfrow=c(4,4))
		hist(com[com$type == "Fremont",84])
		hist(com[com$type == "Fremont",85])
		hist(com[com$type == "Fremont",86])

# Find distribution of abundance for each species that simulated hybrids should have

# Simulate one hybrid community

# Loop simulations to get many samples
