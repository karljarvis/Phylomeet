### Community Phylogenetic Analyses with package picante ###
# R version 2.14.1

# getting set up
#install.packages("picante", dependencies = TRUE)
library(picante)

# Picante manual
vignette("picante-intro")

# dataset
setwd("/Users/kjj/Documents/Phylomeet/Analysis/arcot_gw")

# Cottonwood community with all arthropods, 2000-2003
comm <- data.matrix(read.table("arcot_gw.txt"))

# load nexus files for phylogenies
phy.ultra.full <- read.nexus("ultra.tre")
phy.equal.full <- read.nexus("equal.tre")
phy.grad.full <- read.nexus("grad.tre")

# prune phylogenies to fit community data
phy.ultra <- prune.sample(comm,phy.ultra.full)
phy.equal <- prune.sample(comm,phy.equal.full)
phy.grad <- prune.sample(comm,phy.grad.full)

# Form distance matrix for each topology
phydist.ultra <- cophenetic(phy.ultra)
phydist.equal <- cophenetic(phy.equal)
phydist.grad <- cophenetic(phy.grad)

# plot
par(mfrow=c(1,3))
plot(phy.ultra, cex=0.25) 
title(main = "Ultrametric Branch Lengths")

plot(phy.equal, cex=0.25) 
title(main = "Equal Branch Lengths")

plot(phy.grad, cex=0.25) 
title(main = "Graduated Branch Lengths")

#################################
# Faith's PD
#################################
pd.e <- pd(comm, phy.equal,include.root = TRUE); pd.e
write.csv(pd.e, "pd.equal.csv")

pd.u <- pd(comm, phy.ultra,include.root = TRUE); pd.u
write.csv(pd.u, "pd.ultra.csv")

pd.g <- pd(comm, phy.grad,include.root = TRUE); pd.g
write.csv(pd.g, "pd.grad.csv")

# standardized effect size = ses
ses.pd.u <- ses.pd(comm, phy.ultra,include.root = TRUE, ); ses.pd.u
write.csv(ses.pd.u, "ses.pd.u.csv")


#################################
# NRI 
#################################
# NRI - abundance weighted
# equal
out.ses.mpd <- ses.mpd(comm, phydist.equal, null.model = "richness", abundance.weighted = TRUE, runs = 999)
write.csv(out.ses.mpd,"nri_equal_rich.csv")

# ultrametric
out.ses.mpd <- ses.mpd(comm, phydist.ultra, null.model = "richness", abundance.weighted = TRUE, runs = 999)
write.csv(out.ses.mpd,"nri_ultra_rich.csv")

# graduated
out.ses.mpd <- ses.mpd(comm, phydist.grad, null.model = "richness", abundance.weighted = TRUE, runs = 999)
write.csv(out.ses.mpd,"nri_grad_rich.csv")

#################################
# NTI 
#################################
# NTI - abundance weighted
# equal 
out.ses.mntd <- ses.mntd(comm, phydist.equal, null.model = "richness", abundance.weighted = TRUE, runs = 999)
write.csv(out.ses.mntd,"nti_equal_rich.csv")

# ultrametric
out.ses.mntd <- ses.mntd(comm, phydist.ultra, null.model = "richness", abundance.weighted = TRUE, runs = 999)
write.csv(out.ses.mntd,"nti_ultra_rich.csv")

# graduated
out.ses.mntd <- ses.mntd(comm, phydist.grad, null.model = "richness", abundance.weighted = TRUE, runs = 999)
write.csv(out.ses.mntd,"nti_grad_rich.csv")
#################################
# Phylogenetic Beta Diversity
#################################
# community distances
# equal, grad, gradrand, ultra, ultrarand
comdist.equal <- comdist(comm, phydist.equal, abundance.weighted=TRUE)
comdist.ultra <- comdist(comm, phydist.ultra, abundance.weighted=TRUE)

# write files
m.equal <- as.matrix(comdist.equal)
m.ultra <- as.matrix(comdist.ultra)
write.csv(m.equal, "comdist.equal.csv")
write.csv(m.ultra, "comdist.ultra.csv")

# plot clusters
library(cluster)
comdist.clusters.equal <- hclust(comdist.equal)
comdist.clusters.ultra <- hclust(comdist.ultra)
plot(comdist.clusters.equal, cex=.5)
plot(comdist.clusters.ultra, cex=.5)

# phylo-ordination
# NMDS with package "ecodist"
library(ecodist)
nmds.equal.e <- nmds(comdist.equal, mindim=2, maxdim=2, nits=1)
nmds.equal.min <- nmds.min(nmds.equal.e)
plot(nmds.equal.min, col = "blue")

# NMDS with package "vegan"
library(vegan)
nmds.equal.v <- metaMDS(comdist.equal, autotransform=FALSE)
ordiplot(nmds.equal.v, type = "t")

#################################
# nearest taxon across communities
comdistnt.equal <- comdistnt(comm, phydist.equal, abundance.weighted=TRUE)
comdistnt.ultra <- comdistnt(comm, phydist.ultra, abundance.weighted=TRUE)

# write files
m.equal <- as.matrix(comdistnt.equal)
m.ultra <- as.matrix(comdistnt.ultra)
write.csv(m.equal, "comdistnt.equal.csv")
write.csv(m.ultra, "comdistnt.ultra.csv")

# plot clusters
library(cluster)

comdistnt.clusters.equal <- hclust(comdistnt.equal)
comdistnt.clusters.ultra <- hclust(comdistnt.ultra)

plot(comdistnt.clusters.equal, cex=.5)
plot(comdistnt.clusters.ultra, cex=.5)

#################################
# comm.phylo.cor: Correlations between species co-occurrence and phylogenetic distances
# species.dist metrics:
# cij: Schoener's index of co-occurrence
# jaccard: Jaccard index of co-occurrence
# checkerboard: Checkerboard index of co-occurrence
# doij: DOij index of co-occurrence

out.cor.cij.rich <- comm.phylo.cor(comm, phy.equal, metric = "cij", null.model = "richness", runs = 999)
write.csv(out.cor.cij.rich, "out.cor.cij.rich.csv")

out.cor.cij.freq <- comm.phylo.cor(comm, phy.equal, metric = "cij", null.model = "frequency", runs = 999)
write.csv(out.cor.cij.freq, "out.cor.cij.freq.csv")


out.cor.jac.rich <- comm.phylo.cor(comm, phy.equal, metric = "jaccard", null.model = "richness", runs = 999) # lots of errors, may be because of polytomies

out.cor.doij.rich <- comm.phylo.cor(comm,phy, metric = "doij", null.model = "richness", runs = 999)

out.cor.doij.rich <- comm.phylo.cor(comm,phy, metric = "checkerboard", null.model = "richness", runs = 999)