# getting set up
#install.packages("picante", dependencies = TRUE)
library(picante)

# Picante manual
vignette("picante-intro")

# dataset
setwd("/Users/kjj/Documents/Phylomeet/Analysis/arthrocot")

# Cottonwood community with all arthropods from 2003
comm <- data.matrix(read.table("arthrocot_2003.txt"))

# load nexus files for phylogenies
phy.equal <- read.nexus("/Users/kjj/Documents/Phylomeet/Analysis/arthrocot/arthrocot_equal.trees")
phy.grad <- read.nexus("/Users/kjj/Documents/Phylomeet/Analysis/arthrocot/arthrocot_grad.trees")
phy.gradrand <- read.nexus("/Users/kjj/Documents/Phylomeet/Analysis/arthrocot/arthrocot_gradrand.trees")
phy.ultra <- read.nexus("/Users/kjj/Documents/Phylomeet/Analysis/arthrocot/arthrocot_ultra.trees")
phy.ultrarand <- read.nexus("/Users/kjj/Documents/Phylomeet/Analysis/arthrocot/arthrocot_ultrarand.trees")

# Form distance matrix for each topology
phydist.equal <- cophenetic(phy.equal)
phydist.grad <- cophenetic(phy.grad)
phydist.gradrand <- cophenetic(phy.gradrand)
phydist.ultra <- cophenetic(phy.ultra)
phydist.ultrarand <- cophenetic(phy.ultrarand)


# plot
par(mfrow=c(2,3))

plot(phy.equal, cex=0.25)
title(main = "Equal Branch Lengths")

plot(phy.grad, cex=0.25)
title(main = "Graduated Branch Lengths")

plot(phy.gradrand, cex=0.25) 
title(main = "Randomized Graduated Branch Lengths")

plot(phy.ultra, cex=0.25) 
title(main = "Ultrametric Branch Lengths")

plot(phy.ultrarand, cex=0.25) 
title(main = "Randomized Ultrametric Branch Lengths")


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
out.ses.mpd <- ses.mpd(comm, phydist.equal, null.model = "richness", abundance.weighted = TRUE, runs = 999)
write.csv(out.ses.mpd,"nri_equal_rich.csv")

# graduated
out.ses.mpd <- ses.mpd(comm, phydist.grad, null.model = "richness", abundance.weighted = TRUE, runs = 999)
write.csv(out.ses.mpd,"nri_grad_rich.csv")

# randomized graduated
out.ses.mpd <- ses.mpd(comm, phydist.gradrand, null.model = "richness", abundance.weighted = TRUE, runs = 999)
write.csv(out.ses.mpd,"nri_gradrand_rich.csv")

# ultrametric
out.ses.mpd <- ses.mpd(comm, phydist.ultra, null.model = "richness", abundance.weighted = TRUE, runs = 999)
write.csv(out.ses.mpd,"nri_ultra_rich.csv")

# randomized ultrametric 
out.ses.mpd <- ses.mpd(comm, phydist.ultrarand, null.model = "richness", abundance.weighted = TRUE, runs = 999)
write.csv(out.ses.mpd,"nri_ultrarand_rich.csv")

#################################
# NTI 
#################################
# NTI - abundance weighted
# equal 
out.ses.mntd <- ses.mntd(comm, phydist.equal, null.model = "richness", abundance.weighted = TRUE, runs = 999)
write.csv(out.ses.mntd,"nti_equal_rich.csv")

# graduated
out.ses.mntd <- ses.mntd(comm, phydist.grad, null.model = "richness", abundance.weighted = TRUE, runs = 999)
write.csv(out.ses.mntd,"nti_grad_rich.csv")

# randomized graduated
out.ses.mntd <- ses.mntd(comm, phydist.gradrand, null.model = "richness", abundance.weighted = TRUE, runs = 999)
write.csv(out.ses.mntd,"nti_gradrand_rich.csv")

# ultrametric
out.ses.mntd <- ses.mntd(comm, phydist.ultra, null.model = "richness", abundance.weighted = TRUE, runs = 999)
write.csv(out.ses.mntd,"nti_ultra_rich.csv")

# randomized ultrametric 
out.ses.mntd <- ses.mntd(comm, phydist.ultrarand, null.model = "richness", abundance.weighted = TRUE, runs = 999)
write.csv(out.ses.mntd,"nti_ultrarand_rich.csv")


#################################
# Phylogenetic Beta Diversity
#################################
# community distances
# equal, grad, gradrand, ultra, ultrarand
out.comdist.equal <- comdist(comm, phydist.equal, abundance.weighted=TRUE)
out.comdist.grad <- comdist(comm, phydist.grad, abundance.weighted=TRUE)
out.comdist.gradrand <- comdist(comm, phydist.gradrand, abundance.weighted=TRUE)
out.comdist.ultra <- comdist(comm, phydist.ultra, abundance.weighted=TRUE)
out.comdist.ultrarand <- comdist(comm, phydist.ultrarand, abundance.weighted=TRUE)

# write files
write(out.comdist.equal, "comdist.equal.txt", ncol=nrow(comm))
write(out.comdist.grad, "comdist.grad.txt", ncol=nrow(comm))
write(out.comdist.gradrand, "comdist.gradrand.txt", ncol=nrow(comm))
write(out.comdist.ultra, "comdist.ultra.txt", ncol=nrow(comm))
write(out.comdist.ultrarand, "comdist.ultrarand.txt", ncol=nrow(comm))

# plot clusters
library(cluster)

out.comdist.clusters.equal <- hclust(out.comdist.equal)
out.comdist.clusters.grad <- hclust(out.comdist.grad)
out.comdist.clusters.gradrand <- hclust(out.comdist.gradrand)
out.comdist.clusters.ultra <- hclust(out.comdist.ultra)
out.comdist.clusters.ultrarand <- hclust(out.comdist.ultrarand)

par(mfrow=c(2,3))
plot(out.comdist.clusters.equal, cex=1)
plot(out.comdist.clusters.grad, cex=1)
plot(out.comdist.clusters.gradrand, cex=1)
plot(out.comdist.clusters.ultra, cex=1)
plot(out.comdist.clusters.ultrarand, cex=1)


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