# getting set up
install.packages("picante", dependencies = TRUE)
library(picante)

# Picante manual
vignette("picante-intro")

# sample dataset
data(phylocom)
names(phylocom)

# short names for the datasets
phy <- phylocom$phylo
comm <- phylocom$sample
traits <- phylocom$traits


### phylogeny I'm playing with ###
phy <- read.tree("phy.tree")
comm <- data.matrix(read.table("comm.txt"))

# phylogenies: see package "ape"
phy # gives overview
plot(phy) # gives you a tree

# community data: format as package "vegan"
comm
class(comm)
colnames(comm)
rownames(comm)

# trait data
head(traits)
traitA <- df2vec(traits,"traitA") #df2vec is dataframe to vector conversion, a function of picante
traitA

# Visualizing trees and data
prunedphy <- prune.sample(comm,phy) # cuts down the phylogeny to include only what's in the community
prunedphy

par(mfrow=c(1,2))
plot(phy)
plot(prunedphy)

# plot presence of species in each of 6 communities
par(mfrow = c(2,3))
for (i in row.names(comm)) 
	{
		plot(prunedphy, show.tip.label = FALSE, main = i) # I can't figure out why "main = i" is here. I see no change in the plots when I change or delete it.
		tiplabels(tip = which(comm[i, ] >0), pch = 20, cex = 2) # tiplabels is an 'ape' command
	}

# plot trait values by using different colors
par(mfrow = c(2,2))
for (i in names(traits)) 
	{
		plot(phy, show.tip.label = FALSE, main = i)
		tiplabels(pch = 22, col = traits[, i] + 1, bg = traits[, i] + 1, cex = 1.5) # col is for text/symbol color, and bg = background color. You'd only need bg here when you have an empty symbol, like the empty square (symbol pch=22) 
	}

# Faith's Phylogenetic Diversity
pd.result <- pd(comm, phy, include.root = TRUE)
pd.result #lower PD = more clumped because they capture less of the phylogenetic diversity of the tree. SR = Species Richness

# Webb et al metrics. Require a distance matrix as input rather than a phylogeny object. A phylo object can be converted to a dist matrix via "cophenetic" function. MPD and MNTD can calculate trait diversity measures if you have a trait distance matrix.
# MPD = Mean Phylogenetic Distance
# MNTD = Mean Neares Taxon Distance
# ses.mpd = -NRI (Net Relatedness Index)
# ses.mntd = -NTI (Nearest Taxon Index)
phydist <- cophenetic(phy)
ses.mpd.result <- ses.mpd(comm,phydist, null.model = 'taxa.labels', abundance.weighted = FALSE, runs = 99) #need way more runs for a good analysis
ses.mpd.result
# mpd.obs.z is -NRI

Phylogenetic Beta Diversity
comdist.result <- comdist(comm, phydist)
comdist.result
library(cluster)
comdist.clusters <- hclust(comdist.result)
plot(comdist.clusters)

