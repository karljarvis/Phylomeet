# getting set up
#install.packages("picante", dependencies = TRUE)
library(picante)

# Picante manual
vignette("picante-intro")

# dataset
setwd("/Users/kjj/Documents/_NAU/Classes/G2E/Phylomeet/Analysis/arthro_cot")
fullphy <- read.nexus("arthro_cot_ultra.trees")
plot(fullphy, cex=0.25)

# Community with all arthropods
#comm <- data.matrix(read.csv("arthro_cot_2003.csv"))

# Community with insects only
comm <- data.matrix(read.table("arthro_cot_2003.txt"))

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

# NRI - abundance weighted
# independent swap randomization
phydist <- cophenetic(phy)
out.ses.mpd <- ses.mpd(comm, phydist, cstSor = TRUE, null.model = "independentswap", abundance.weighted = TRUE, runs = 999)
write.csv(out.ses.mpd,"nri_ultra_indswap.csv")

# NRI - abundance weighted
# trial swap randomization
phydist <- cophenetic(phy)
out.ses.mpd <- ses.mpd(comm, phydist, cstSor = FALSE, null.model = "trialswap", abundance.weighted = TRUE, runs = 999)
write.csv(out.ses.mpd,"nri_ultra_triswap.csv")
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
#out.comdist.uw <- comdist(comm, phydist, abundance.weighted=FALSE)
#
#library(cluster)
#phy.comdist.uw.clusters <- hclust(phy.comdist)
#plot(phy.comdist.uw.clusters, cex=0.5)

# abundance weighted
out.comdist <- comdist(comm, phydist, abundance.weighted=TRUE)
write(out.comdist, "comdistall.txt", ncol=nrow(comm))

library(cluster)
out.comdist.clusters <- hclust(out.comdist)
plot(out.comdist.clusters, cex=1)

#################################
# nearest taxon across communities
phydist <- cophenetic(phy)
out.comdistnt <- comdistnt(comm, phydist, abundance.weighted=TRUE)
library(cluster)
out.comdistnt.clusters <- hclust(out.comdistnt)
plot(out.comdistnt.clusters, cex=1)
write(out.comdistnt, "comdistntall.txt", ncol=nrow(comm))

#################################
# comm.phylo.cor: Correlations between species co-occurrence and phylogenetic distances
# species.dist metrics:
# cij: Schoener's index of co-occurrence
# jaccard: Jaccard index of co-occurrence
# checkerboard: Checkerboard index of co-occurrence
# doij: DOij index of co-occurrence

out.cor.cij.rich <- comm.phylo.cor(comm, phy, metric = "cij", null.model = "richness", runs = 999)

out.cor.cij.freq <- comm.phylo.cor(comm, phy, metric = "cij", null.model = "frequency", runs = 999)

out.cor.jac.rich <- comm.phylo.cor(comm,phy, metric = "jaccard", null.model = "richness", runs = 999) # lots of errors, may be because of polytomies

out.cor.doij.rich <- comm.phylo.cor(comm,phy, metric = "doij", null.model = "richness", runs = 999)

out.cor.doij.rich <- comm.phylo.cor(comm,phy, metric = "checkerboard", null.model = "richness", runs = 999)