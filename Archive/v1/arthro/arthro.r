# getting set up
#install.packages("picante", dependencies = TRUE)
library(picante)

# Picante manual
vignette("picante-intro")

# dataset
setwd("/Users/kjj/Documents/_NAU/Classes/G2E/Phylomeet/Analysis/arthro/")
phy <- read.nexus("arthro.nexus")
phy <- read.nexus("arthroultra.trees")
plot(phy, cex=0.25)
comm <- data.matrix(read.table("arthrocomm.txt"))
trait <- read.table("arthrotrait.txt")

# prune tree - we don't have the need for this
#comm.phy <- prune.sample(comm,phy) # cuts down the phylogeny to include only what's in the community
#comm.phy

#################################
# Faith's Phylogenetic Diversity
#################################
# regular PD
pd.phy <- pd(comm, phy, include.root = TRUE); pd.phy 
write.table(pd.phy,"phylodivult.txt")

# standardized PD, richness null model
ses.pd.phy.richness <- ses.pd(comm, phy, null.model="richness", runs=999)
write.table(ses.pd.phy.richness,"PD_stand_richnessult.txt")

#################################
# NRI 
#################################
# NRI - abundance weighted
phydist <- cophenetic(phy)
ses.mpd.arthro.w <- ses.mpd(comm, phydist, null.model = 'richness', abundance.weighted = TRUE, runs = 999); ses.mpd.arthro.w
write.table(ses.mpd.arthro.w,"NRI_weighted_rich.txt")

#################################
# NTI 
#################################
# abundance weighted
phydist <- cophenetic(phy)
ses.mntd.arthro.w <- ses.mntd(comm, phydist, null.model = 'richness', abundance.weighted = TRUE, runs = 999); ses.mntd.arthro.w
write.table(ses.mntd.arthro.w,"NTI_weighted.txt")

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
# nearest taxon across communities
comdistnt.result.w <- comdistnt(comm, phydist, abundance.weighted=TRUE)
comdistnt.result.w
library(cluster)
comdistnt.clusters.w <- hclust(comdistnt.result.w)
plot(comdistnt.clusters.w, cex=0.5)
write(comdistnt.result.w, "phylobeta.txt", ncol=nrow(comm))

#################################
# K stats
multiPhylosignal(trait,phy)
