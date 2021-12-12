# getting set up
#install.packages("picante", dependencies = TRUE)
library(picante)

# Picante manual
vignette("picante-intro")

# dataset
setwd("/Users/kjj/Documents/_NAU/Classes/G2E/Phylomeet/Analysis/arthrocot")

# Cottonwood community with all arthropods from 2002
comm <- data.matrix(read.table("arthrocot_2002.txt"))

# load nexus files for phylogenies
phy_equal <- read.nexus("/Users/kjj/Documents/_NAU/Classes/G2E/Phylomeet/Analysis/test/ultra.trees")

# plot
plot(phy_equal)
title(main = "Equal Branch Lengths")

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
# equal
phydist_equal <- cophenetic(phy_equal)
out.ses.mpd <- ses.mpd(comm, phydist_equal, null.model = "richness", abundance.weighted = TRUE, runs = 999)
write.csv(out.ses.mpd,"nri_equal_rich.csv")

# graduated
phydist_grad <- cophenetic(phy_grad)
out.ses.mpd <- ses.mpd(comm, phydist_grad, null.model = "richness", abundance.weighted = TRUE, runs = 999)
write.csv(out.ses.mpd,"nri_grad_rich.csv")

# randomized graduated
phydist_gradrand <- cophenetic(phy_gradrand)
out.ses.mpd <- ses.mpd(comm, phydist_gradrand, null.model = "richness", abundance.weighted = TRUE, runs = 999)
write.csv(out.ses.mpd,"nri_gradrand_rich.csv")

# ultrametric
phydist_ultra <- cophenetic(phy_ultra)
out.ses.mpd <- ses.mpd(comm, phydist_ultra, null.model = "richness", abundance.weighted = TRUE, runs = 999)
write.csv(out.ses.mpd,"nri_ultra_rich.csv")

# randomized ultrametric 
phydist_ultrarand <- cophenetic(phy_ultrarand)
out.ses.mpd <- ses.mpd(comm, phydist_ultrarand, null.model = "richness", abundance.weighted = TRUE, runs = 999)
write.csv(out.ses.mpd,"nri_ultrarand_rich.csv")

# NRI - abundance weighted
# independent swap randomization
#phydist <- cophenetic(phy)
#out.ses.mpd <- ses.mpd(comm, phydist, cstSor = TRUE, null.model = "independentswap", abundance.weighted = TRUE, runs = 999, iterations = 10)
#write.csv(out.ses.mpd,"nri_ultra_indswap.csv")

# NRI - abundance weighted
# trial swap randomization
#phydist <- cophenetic(phy)
#out.ses.mpd <- ses.mpd(comm, phydist, cstSor = FALSE, null.model = "trialswap", abundance.weighted = TRUE, runs = 999)
#write.csv(out.ses.mpd,"nri_ultra_triswap.csv")
#################################
# NTI 
#################################
# NTI - abundance weighted
# equal 
phydist_equal <- cophenetic(phy_equal)
out.ses.mntd <- ses.mntd(comm, phydist_equal, null.model = "richness", abundance.weighted = TRUE, runs = 999)
write.csv(out.ses.mntd,"nti_equal_rich.csv")

# graduated
phydist_grad <- cophenetic(phy_grad)
out.ses.mntd <- ses.mntd(comm, phydist_grad, null.model = "richness", abundance.weighted = TRUE, runs = 999)
write.csv(out.ses.mntd,"nti_grad_rich.csv")

# randomized graduated
phydist_gradrand <- cophenetic(phy_gradrand)
out.ses.mntd <- ses.mntd(comm, phydist_gradrand, null.model = "richness", abundance.weighted = TRUE, runs = 999)
write.csv(out.ses.mntd,"nti_gradrand_rich.csv")

# ultrametric
phydist_ultra <- cophenetic(phy_ultra)
out.ses.mntd <- ses.mntd(comm, phydist_ultra, null.model = "richness", abundance.weighted = TRUE, runs = 999)
write.csv(out.ses.mntd,"nti_ultra_rich.csv")

# randomized ultrametric 
phydist_ultrarand <- cophenetic(phy_ultrarand)
out.ses.mntd <- ses.mntd(comm, phydist_ultrarand, null.model = "richness", abundance.weighted = TRUE, runs = 999)
write.csv(out.ses.mntd,"nti_ultrarand_rich.csv")

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