# getting set up
#install.packages("picante", dependencies = TRUE)
library(picante)

# Picante manual
vignette("picante-intro")

# dataset
setwd("/Users/kjj/Documents/_NAU/Classes/G2E/Phylomeet/Analysis/arthro_cot_bl")
fullphy <- read.nexus("arthro_cot_bltest.trees")
plot(fullphy, cex=0.25)

# Community with all arthropods
#comm <- data.matrix(read.csv("arthro_cot_2003.csv"))

# Community with insects only
comm <- data.matrix(read.table("arthro_cot_2003_insect.txt"))

#trait <- read.table("arthrotrait.txt")
#write.table(comm[,36:111], "acom.txt")

# prune tree
phy <- prune.sample(comm,fullphy) # cuts down the phylogeny to include only what's in the community
plot(phy, cex=0.4)

#################################
# Faith's Phylogenetic Diversity
#################################
# regular PD
out.pd <- pd(comm, phy, include.root = TRUE) 
write.csv(out.pd,"pd_u.csv")

# standardized PD, richness null model
out.ses.pd <- ses.pd(comm, phy, null.model="richness", runs=999)
write.csv(out.ses.pd,"pds_u.csv")

#################################
# NRI 
#################################
# NRI - abundance weighted
phydist <- cophenetic(phy)
out.ses.mpd <- ses.mpd(comm, phydist, null.model = "richness", abundance.weighted = TRUE, runs = 999)
write.csv(out.ses.mpd,"nri_ultra_rich.csv")

#################################
# NTI 
#################################
# abundance weighted
phydist <- cophenetic(phy)
out.ses.mntd <- ses.mntd(comm, phydist, null.model = 'richness', abundance.weighted = TRUE, runs = 999); write.csv(out.ses.mntd,"nti_ultra_rich.csv")

#################################
# Phylogenetic Beta Diversity
#################################
phydist <- cophenetic(phy)

# not abundance weighted
out.comdist.uw <- comdist(comm, phydist, abundance.weighted=FALSE)
write(out.comdist.uw, "comdist.uw.txt", ncol=nrow(comm))
library(cluster)
out.comdist.uw.clusters <- hclust(out.comdist.uw)
plot(out.comdist.uw.clusters, cex=1)

# abundance weighted
out.comdist <- comdist(comm, phydist, abundance.weighted=TRUE)
write(out.comdist, "comdist.txt", ncol=nrow(comm))
library(cluster)
out.comdist.clusters <- hclust(out.comdist)
plot(out.comdist.clusters, cex=1)

#################################
# nearest taxon across communities

# not abundance weighted
phydist <- cophenetic(phy)
out.comdistnt.uw <- comdistnt(comm, phydist, abundance.weighted=FALSE)
library(cluster)
out.comdistnt.clusters <- hclust(out.comdistnt)
plot(out.comdistnt.clusters, cex=1)
write(out.comdistnt, "comdistnt.txt", ncol=nrow(comm))

# abundance weighted
phydist <- cophenetic(phy)
out.comdistnt <- comdistnt(comm, phydist, abundance.weighted=TRUE)
library(cluster)
out.comdistnt.clusters <- hclust(out.comdistnt)
plot(out.comdistnt.clusters, cex=1)
write(out.comdistnt, "comdistnt.txt", ncol=nrow(comm))

#################################
# K stats
multiPhylosignal(trait,phy)
