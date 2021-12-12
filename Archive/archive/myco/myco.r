# getting set up
#install.packages("picante", dependencies = TRUE)
library(picante)

# Picante manual
vignette("picante-intro")

# dataset
setwd("/Users/karljarvis/Documents/_NAU/Classes/G2E/Phylomeet/Analysis/myco")
phy <- read.tree("mycultra.phy")
plot(phy)
comm <- data.matrix(read.table("myccom.txt"))

# prune tree - we don't have the need for this
#comm.phy <- prune.sample(comm,phy) # cuts down the phylogeny to include only what's in the community
#comm.phy

#################################
# Faith's Phylogenetic Diversity
#################################
# regular PD
pd.phy <- pd(comm, phy, include.root = TRUE); pd.phy 
write.table(pd.phy,"phylodivul.txt")

# standardized PD, richness null model
ses.pd.phy.richness <- ses.pd(comm, phy, null.model="richness", runs=999)
write.table(ses.pd.phy.richness,"PD_stand_richness.txt")

# standardized PD, Gotelli null model
#ses.pd.phy.indswap <- ses.pd(comm, phy, null.model="independentswap", runs=999, iterations=1000)

#################################
# NRI 
#################################
# abundance not weighted
phydist <- cophenetic(phy)
ses.mpd.myco.uw <- ses.mpd(comm, phydist, null.model = 'taxa.labels', abundance.weighted = FALSE, runs = 999); ses.mpd.myco.uw
write.table(ses.mpd.myco.uw,"NRI_unweighted.txt")

# abundance weighted
phydist <- cophenetic(phy)
ses.mpd.myco.w <- ses.mpd(comm, phydist, null.model = 'taxa.labels', abundance.weighted = TRUE, runs = 999); ses.mpd.myco.w
write.table(ses.mpd.myco.w,"NRI_weighted.txt")

# abundance not weighted
phydist <- cophenetic(phy)
ses.mpd.myco.uw <- ses.mpd(comm, phydist, null.model = 'richness', abundance.weighted = FALSE, runs = 999); ses.mpd.myco.uw
write.table(ses.mpd.myco.uw,"NRI_unweighted_rich.txt")

# NRI - abundance weighted
phydist <- cophenetic(phy)
ses.mpd.myco.w <- ses.mpd(comm, phydist, null.model = 'richness', abundance.weighted = TRUE, runs = 999); ses.mpd.myco.w
write.table(ses.mpd.myco.w,"NRI_weighted_rich.txt")

#################################
# NTI 
#################################
# abundance not weighted
phydist <- cophenetic(phy)
ses.mntd.myco.uw <- ses.mntd(comm, phydist, null.model = 'taxa.labels', abundance.weighted = FALSE, runs = 999); ses.mntd.myco.uw
write.table(ses.mntd.myco.uw,"NTI_unweighted.txt")

# abundance weighted
phydist <- cophenetic(phy)
ses.mntd.myco.w <- ses.mntd(comm, phydist, null.model = 'taxa.labels', abundance.weighted = TRUE, runs = 999); ses.mntd.myco.w
write.table(ses.mntd.myco.w,"NTI_weighted.txt")

#################################
# Phylogenetic Beta Diversity
#################################
phydist <- cophenetic(phy)

# not abundance weighted
comdist.result.uw <- comdist(comm, phydist, abundance.weighted=FALSE)
comdist.result.uw
library(cluster)
comdist.clusters.uw <- hclust(comdist.result.uw)
plot(comdist.clusters.uw, cex=0.5)

# abundance weighted
comdist.result.w <- comdist(comm, phydist, abundance.weighted=TRUE)
comdist.result.w
library(cluster)
comdist.clusters.w <- hclust(comdist.result.w)
plot(comdist.clusters.w, cex=0.5)
write(comdist.result.w, "phylobeta.txt", ncol=nrow(comm))

#################################
# Phylogenetic Beta Diversity NT
#################################
phydist <- cophenetic(phy)

# not abundance weighted
comdistnt.result.uw <- comdistnt(comm, phydist, abundance.weighted=FALSE)
comdistnt.result.uw
library(cluster)
comdistnt.clusters.uw <- hclust(comdistnt.result.uw)
plot(comdistnt.clusters.uw, cex=0.5)

# abundance weighted
comdist.result.w <- comdist(comm, phydist, abundance.weighted=TRUE)
comdist.result.w
library(cluster)
comdist.clusters.w <- hclust(comdist.result.w)
plot(comdist.clusters.w, cex=0.5)
write(comdist.result.w, "phylobeta.txt", ncol=nrow(comm))


#################################
# Phylosor - fraction of branch length shared between two communities
#################################
# this makes less sense for our data, since it's more branch-length heavy.
phylosor.result <- phylosor(comm,phy)
phylosor.clusters <- hclust(phylosor.result)
plot(phylosor.clusters, cex=0.5)

#################################
# Unifrac -unweighted
#################################
unifrac.result <- unifrac(comm,phy)