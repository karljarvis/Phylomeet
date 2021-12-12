# getting set up
#install.packages("picante", dependencies = TRUE)
library(picante)

# Picante manual
vignette("picante-intro")




### branch lengths all equal###
setwd("/Users/karljarvis/Documents/_NAU/Classes/G2E/Phylomeet/picante_tests")
phy <- read.tree("phy.tree")
comm <- data.matrix(read.table("comm.txt"))

# Faith's Phylogenetic Diversity
pd.result <- pd(comm, phy, include.root = TRUE); pd.result

# Webb et al metrics: mpd.obs.z = -NRI
phydist <- cophenetic(phy)
ses.mpd.result <- ses.mpd(comm,phydist, null.model = 'taxa.labels', abundance.weighted = FALSE, runs = 999) 
ses.mpd.result

# Webb et al metrics: mpd.obs.z = -NRI
phydist <- cophenetic(phy)
ses.mpd.result <- ses.mpd(comm,phydist, null.model = 'taxa.labels', abundance.weighted = TRUE, runs = 999) 
ses.mpd.result


# Phylogenetic Beta Diversity
comdist.result <- comdist(comm, phydist)
comdist.result
library(cluster)
comdist.clusters <- hclust(comdist.result)




### branch lengths .1 and 10 ###
phy1 <- read.tree("phy1.tree")
comm <- data.matrix(read.table("comm.txt"))

# Faith's Phylogenetic Diversity
pd1.result <- pd(comm, phy1, include.root = TRUE); pd1.result
#lower PD = more clumped because they capture less of the phylogenetic diversity of the tree. SR = Species Richness

# Webb et al metrics: mpd.obs.z = -NRI
phydist1 <- cophenetic(phy1)
ses1.mpd.result <- ses.mpd(comm,phydist1, null.model = 'taxa.labels', abundance.weighted = FALSE, runs = 999) 
ses1.mpd.result

phydist1 <- cophenetic(phy1)
ses1a.mpd.result <- ses.mpd(comm,phydist1, null.model = 'taxa.labels', abundance.weighted = TRUE, runs = 999) 
ses1a.mpd.result


# Phylogenetic Beta Diversity
comdist1.result <- comdist(comm, phydist1)
comdist1.result
library(cluster)
comdist1.clusters <- hclust(comdist1.result)



### branch lengths  .5 and 2 ###
phy2 <- read.tree("phy2.tree")
comm <- data.matrix(read.table("comm.txt"))

# Faith's Phylogenetic Diversity
pd2.result <- pd(comm, phy2, include.root = TRUE); pd2.result

# Webb et al metrics: mpd.obs.z = -NRI
phydist2 <- cophenetic(phy2)
ses2.mpd.result <- ses.mpd(comm,phydist2, null.model = 'taxa.labels', abundance.weighted = FALSE, runs = 999)
ses2.mpd.result

phydist2 <- cophenetic(phy2)
ses2.mpd.result <- ses.mpd(comm,phydist2, null.model = 'taxa.labels', abundance.weighted = TRUE, runs = 999)
ses2.mpd.result

# Phylogenetic Beta Diversity
comdist2.result <- comdist(comm, phydist2)
comdist2.result
library(cluster)
comdist2.clusters <- hclust(comdist2.result)


# save the stuff
par(mfrow=c(1,3))
plot(phy);plot(phy1);plot(phy2)

par(mfrow=c(1,3))
plot(comdist.clusters, main="phy cluster dendrogram"); plot(comdist1.clusters,main="phy1 cluster dendrogram"); plot(comdist2.clusters, main="phy2 cluster dendrogram")


write.table(c(pd.result, pd1.result, pd2.result), "pdresults.txt")
write.table(c(ses.mpd.result, ses1.mpd.result, ses2.mpd.result), "mpdresults.txt")
write.table(c(comdist.result, comdist1.result, comdist2.result), "comdistresults.txt")
