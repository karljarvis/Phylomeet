##########
# Calculate statistics and plot figures of  community phylogenetic analyses of arthropod communities on cottonwood hosts
# Karl Jarvis June 13, 2014
##########

source("/Users/kjj/Projects/Phylomeet/Analysis/arcot_dataprep_v1.R")
library(ggplot2)

colors = c("black","red","blue")
typelabels = c("Fremont", "Hybrid", "Narrowleaf")
types = rep(crosstypes, 4)
yrs = rep(years, each=3)

#######################
# Figure 1: Phylogeny
phyLongNames <- phyShortNames <- phyList$ultra0
inNames <- read.csv("~/Projects/Phylomeet/Manuscript/Tables&Figures/phylogeny/phynames.csv", colClasses="character")
abbr = inNames[,1]

# phylogeny with order names included
longNames <- apply(inNames[,2:5], 1, function(x) paste(x, collapse="_"))
for(i in 0:9) 
{ 
  longNames <- gsub(i, "sp.", longNames) 
}
longNames <- gsub("sp.sp.", "sp.", longNames)
longNames <- gsub("__", "_", longNames)
longNames <- gsub(" ", "", longNames)

# phylogeny with order names excluded	
shortNames <- apply(inNames[,3:5],1,function(x) paste(x, collapse="_"))
noGenus <- c(188,189,193,199)

for(i in noGenus)
{ 
  shortNames[noGenus] <- paste0(inNames[noGenus,2], "_sp.") 
}

for(i in 0:9)
{ 
  shortNames <- gsub(i, "sp.", shortNames) 
}

shortNames <- gsub("__", "_", shortNames)
shortNames <- gsub("sp.sp.", "sp.", shortNames)

phyLongNames$tip.label <- longNames[match(phyLongNames$tip.label, abbr)]
phyShortNames$tip.label <- shortNames[match(phyShortNames$tip.label, abbr)]

# plotting
pdf("Fig1_phylo_with_orders.pdf", width=8.5, height=11)
plot(phyLongNames, cex=0.4, root.edge=T, no.margin=T)
nodelabels(1:phyLongNames$Nnode, bg="black", col="white", cex=0.5, frame="circle")
box()
dev.off()

pdf("Fig1_phylo.pdf", width=8.5, height=11)
plot(phyShortNames, cex=0.4, edge.width=1.5, root.edge=T, no.margin=T)
nodelabels(1:phyShortNames$Nnode, bg="black", col="white", cex=0.5, frame="circle")
box()
dev.off()

tiff("Fig1_phylo.tiff", width=8.5, height=11, units="in", res=300)
plot(phyShortNames, cex=0.3, edge.width=1.5, root.edge=T, no.margin=T)
nodelabels(1:phyShortNames$Nnode, bg="black", col="white", cex=0.3, frame="circle")
box()
dev.off()

# plot(phyfig, cex=0.5, edge.width=1.5, root.edge=T, no.margin=T, type="fan")



##########
# Figure 2: Faith's Phylogenetic Diversity (PD)
##########

####
# PD individual analyses

# putting individual tree PDs into long format
PDindiv = matrix(NA, 0, 5)
for (i in 1:length(phyNames)) 
{   
  pd = read.csv(paste0(PDdir, "PD_indiv_", phyNames[i], ".csv"), row.names=1)
  df = data.frame(phylo=phyNames[i], crosstype, year, pd)
  PDindiv = rbind(PDindiv, df)
}
PDindiv$year = factor(PDindiv$year)

# plot boxplots
pdf(file=paste0(figDir, "PD_indiv_boxplot.pdf"), width=10, height=8)
ggplot(PDindiv, aes(year,PD)) + 
  geom_boxplot(aes(fill=crosstype)) +
  scale_fill_manual(name = "Cross Type", values = colors, labels = typelabels) +
  facet_wrap(~ phylo) +
  ggtitle("Phylogenetic Diversity\nAnalyses of Communities on Individual Host Trees")
dev.off()

# violin plots
pdf(file=paste0(figDir, "PD_indiv_violin.pdf"), width=10, height=8)
ggplot(PDindiv, aes(year,PD)) + 
  geom_violin(aes(fill=crosstype)) +
  scale_fill_manual(name = "Cross Type", values = colors, labels = typelabels) +
  facet_wrap(~ phylo) +
  ggtitle("Phylogenetic Diversity\nAnalyses of Communities on Individual Host Trees")
dev.off()

pdf(file=paste0(figDir, "PD_indiv_violin_ultraequal.pdf"), width=6, height=4)
ggplot(PDindiv[PDindiv$phylo == "ultra0"|PDindiv$phylo == "equal0",], aes(year,PD)) + 
  geom_violin(aes(fill=crosstype)) +
  scale_fill_manual(name = "Cross Type", values = colors, labels = typelabels) +
  facet_wrap(~ phylo) +
  ggtitle("Phylogenetic Diversity\nAnalyses of Communities on Individual Host Trees")
dev.off()


##########
# Plot PD pooled by crosstype within years
PDpool = matrix(NA, 0, 6)
for (i in 1:length(phyNames)) 
{   
  pd = read.csv(paste0(PDdir, "PD_pooled_", phyNames[i], ".csv"), row.names=1)
  df = data.frame(phylo=phyNames[i], crosstype=types, year=yrs, pd)
  PDpool = rbind(PDpool, df)
}
PDpool$year = factor(PDpool$year)

# plot line plots
pdf(file=paste0(figDir, "PD_pooled_line.pdf"), width=10, height=8)
ggplot(PDpool, aes(year, PD, group=crosstype, color=crosstype)) + 
  geom_line() +
  geom_point() +
  scale_color_manual(name = "Cross Type", values = colors, labels = typelabels) +
  facet_wrap(~ phylo) +
  ggtitle("Phylogenetic Diversity\nAnalyses of Communities Pooled by Host Type")
dev.off()

pdf(file=paste0(figDir, "PD_pooled_line_ultraequal.pdf"), width=6, height=4)
ggplot(PDpool[PDpool$phylo == "ultra0"|PDpool$phylo == "equal0",], aes(year, PD, group=crosstype, color=crosstype)) + 
  geom_line() +
  geom_point() +
  scale_color_manual(name = "Cross Type", values = colors, labels = typelabels) +
  facet_wrap(~ phylo) +
  ggtitle("Phylogenetic Diversity\nAnalyses of Communities Pooled by Host Type")
dev.off()


##########
# Net Relatedness Index (NRI) 
##########

# putting individual tree PDs into long format
NRIindivA = NRIindivP = matrix(NA, 0, 4)
for (i in 1:length(phyNames)) 
{   
  nria = read.csv(paste0(NRIdir, "NRI_abund_indiv_", phyNames[i], ".csv"), row.names=1)
  nrip = read.csv(paste0(NRIdir, "NRI_pres_indiv_", phyNames[i], ".csv"), row.names=1)
  dfa = data.frame(phylo=phyNames[i], crosstype, year, NRI=-nria$mpd.obs.z)
  dfp = data.frame(phylo=phyNames[i], crosstype, year, NRI=-nrip$mpd.obs.z)
  NRIindivA = rbind(NRIindivA, dfa)
  NRIindivP = rbind(NRIindivP, dfp)
}
NRIindivA$year = NRIindivP$year = factor(NRIindivA$year)

# Violin plots
pdf(file=paste0(figDir, "NRI_indiv_abund_violin.pdf"), width=10, height=8)
ggplot(NRIindivA, aes(year,NRI)) + 
  geom_violin(aes(fill=crosstype)) +
  scale_fill_manual(name = "Cross Type", values = colors, labels = typelabels) +
  facet_wrap(~ phylo) +
  ggtitle("Abundance-Weighted Net Relatedness Index\nAnalyses of Communities on Individual Host Trees")
dev.off()

pdf(file=paste0(figDir, "NRI_indiv_pres_violin.pdf"), width=10, height=8)
ggplot(NRIindivP, aes(year,NRI)) + 
  geom_violin(aes(fill=crosstype)) +
  scale_fill_manual(name = "Cross Type", values = colors, labels = typelabels) +
  facet_wrap(~ phylo) +
  ggtitle("Unweighted Net Relatedness Index\nAnalyses of Communities on Individual Host Trees")
dev.off()

# 4 panel figure
NRIindAUE = data.frame(abund="abundance", NRIindivA[NRIindivA$phylo == "ultra0" | NRIindivA$phylo == "equal0",])
NRIindPUE = data.frame(abund="presence", NRIindivP[NRIindivP$phylo == "ultra0" | NRIindivP$phylo == "equal0",])
NRIindUE = rbind(NRIindAUE, NRIindPUE)

pdf(file=paste0(figDir, "NRI_indiv_pa_violin_ultraequal.pdf"), width=10, height=8)
ggplot(NRIindUE, aes(year,NRI)) + 
  geom_violin(aes(fill=crosstype)) +
  scale_fill_manual(name = "Cross Type", values = colors, labels = typelabels) +
  facet_wrap(~ abund + phylo) +
  ggtitle("Net Relatedness Index\nAnalyses of Communities on Individual Host Trees")
dev.off()



##########
# Plot NRI pooled by crosstype within years
NRIpoolA = NRIpoolP = matrix(NA, 0, 4)
for (i in 1:length(phyNames)) 
{   
  nria = read.csv(paste0(NRIdir, "NRI_abund_pooled_", phyNames[i], ".csv"), row.names=1)
  nrip = read.csv(paste0(NRIdir, "NRI_pres_pooled_", phyNames[i], ".csv"), row.names=1)
  dfa = data.frame(phylo=phyNames[i], crosstype=types, year=yrs, NRI=-nria$mpd.obs.z)
  dfp = data.frame(phylo=phyNames[i], crosstype=types, year=yrs, NRI=-nrip$mpd.obs.z)
  NRIpoolA = rbind(NRIpoolA, dfa)
  NRIpoolP = rbind(NRIpoolP, dfp)
}
NRIpoolA$year = NRIpoolP$year = factor(NRIpoolA$year)

# plot line plots
pdf(file=paste0(figDir, "NRI_pooled_abund_line.pdf"), width=10, height=8)
ggplot(NRIpoolA, aes(year, NRI, group=crosstype, color=crosstype)) + 
  geom_line() +
  geom_point() +
  scale_color_manual(name = "Cross Type", values = colors, labels = typelabels) +
  facet_wrap(~ phylo) +
  ggtitle("Abundance-Weighted Net Relatedness Index\nAnalyses of Communities Pooled by Host Type")
dev.off()

pdf(file=paste0(figDir, "NRI_pooled_pres_line.pdf"), width=10, height=8)
ggplot(NRIpoolP, aes(year, NRI, group=crosstype, color=crosstype)) + 
  geom_line() +
  geom_point() +
  scale_color_manual(name = "Cross Type", values = colors, labels = typelabels) +
  facet_wrap(~ phylo) +
  ggtitle("Unweighted Net Relatedness Index\nAnalyses of Communities Pooled by Host Type")
dev.off()

# 4 panel figure
nrip = cbind("pres", NRIpoolP[NRIpoolP$phylo == "ultra0" | NRIpoolP$phylo == "equal0",])
nria = cbind("abund", NRIpoolA[NRIpoolA$phylo == "ultra0" | NRIpoolA$phylo == "equal0",])
names(nrip)[1] = names(nria)[1] = "abund"
n = rbind(nrip, nria)

pdf(file=paste0(figDir, "NRI_pooled_pa_line_ultraequal.pdf"), width=10, height=8)
ggplot(n, aes(year, NRI, group=crosstype, color=crosstype)) + 
  geom_line() +
  geom_point() +
  scale_color_manual(name = "Cross Type", values = colors, labels = typelabels) +
  facet_wrap(~ abund + phylo) +
  ggtitle("Unweighted Net Relatedness Index\nAnalyses of Communities Pooled by Host Type")
dev.off()


##########
# Phylogenetic Beta Diversity: Community Distance (ComDist)
##########

# Load comdist distance matrices into lists
comDistA = comDistP = vector("list",length=length(comList))
names(comDistA) = names(comDistP) = names(comList)
for (i in 1:length(comList)) 
{	
	comDistA[[i]] = comDistP[[i]] = vector("list",length=length(phyNames)) 
	names(comDistA[[i]]) = names(comDistP[[i]]) = phyNames
	for (j in 1:length(phyNames)) 
	{ 	
		comDistA[[i]][[j]] = dist(read.csv(paste0(comDistDir, "ComDist_abund_", names(comList)[i], "_", phyNames[j], ".csv"), row.names=1)) 
		comDistP[[i]][[j]] = dist(read.csv(paste0(comDistDir, "ComDist_pres_", names(comList)[i], "_", phyNames[j], ".csv"), row.names=1)) 
	}
	names(comDistA[[i]]) = names(comDistP[[i]]) = phyNames
}


##########
# Nonmetric MDS 2D
library(MASS)
library(car)
library(devtools)
library(digest)

# individual trees
type = rep(c("Fre","F1","Nar"),4)		# names for labels
cols

typelabels
colors
syms = rep(c(15,16,17),4)

# load data
CDindivA = CDindivP = matrix(NA, 0, 5)
for(i in 1:length(comDistA$indiv))
{
  fitA = isoMDS(comDistA$indiv[[i]], k=2)
  dfA = data.frame(phylo=phyNames[i], crosstype, year, Axis1=fitA$points[,1], Axis2=fitA$points[,2])
  CDindivA = rbind(CDindivA, dfA)

  fitP = isoMDS(comDistP$indiv[[i]], k=2)
  dfP = data.frame(phylo=phyNames[i], crosstype, year, Axis1=fitP$points[,1], Axis2=fitP$points[,2])
  CDindivP = rbind(CDindivP, dfP)
}
CDindivA$year = CDindivP$year = factor(CDindivA$year)

# plot
pdf(file=paste0(figDir, "comDist_indiv_abund_dot.pdf"), width=10, height=10)
ggplot(CDindivA, aes(Axis1, Axis2, color=crosstype)) +
  geom_point() +
  stat_ellipse() +
  scale_color_manual(name = "Cross Type", values = colors, labels = typelabels) +
  facet_wrap(~ phylo) +
  ggtitle("Abundance-Weighted Community Distance\nAnalyses of Communities on Individual Host Trees")
dev.off()

pdf(file=paste0(figDir, "comDist_indiv_pres_dot.pdf"), width=10, height=10)
ggplot(CDindivP, aes(Axis1, Axis2, color=crosstype)) +
  geom_point() +
  stat_ellipse() +
  scale_color_manual(name = "Cross Type", values = colors, labels = typelabels) +
  facet_wrap(~ phylo) +
  ggtitle("Unweighted Community Distance\nAnalyses of Communities on Individual Host Trees")
dev.off()

# 4 panel figure
cdp = cbind("pres", CDindivP[CDindivP$phylo == "ultra0" | CDindivP$phylo == "equal0",])
cda = cbind("abund", CDindivA[CDindivA$phylo == "ultra0" | CDindivA$phylo == "equal0",])
names(cdp)[1] = names(cda)[1] = "abund"
cd = rbind(cdp, cda)

pdf(file=paste0(figDir, "comDist_indiv_pa_dot_ultraequal.pdf"), width=10, height=8)
ggplot(cd, aes(Axis1, Axis2, color=crosstype)) +
  geom_point() +
  stat_ellipse() +
  scale_color_manual(name = "Cross Type", values = colors, labels = typelabels) +
  facet_wrap(~ abund + phylo) +
  ggtitle("Community Distance\nAnalyses of Communities on Individual Host Trees")
dev.off()



#########
# analyses with data pooled by host type within each year

# load data
CDpoolA = CDpoolP = matrix(NA, 0, 5)
for(i in 1:length(comDistA$pooled))
{
  fitA = isoMDS(comDistA$pooled[[i]], k=2)
  dfA = data.frame(phylo=phyNames[i], crosstype=types, year=yrs, Axis1=fitA$points[,1], Axis2=fitA$points[,2])
  CDpoolA = rbind(CDpoolA, dfA)
  
  fitP = isoMDS(comDistP$pooled[[i]], k=2)
  dfP = data.frame(phylo=phyNames[i], crosstype=types, year=yrs, Axis1=fitP$points[,1], Axis2=fitP$points[,2])
  CDpoolP = rbind(CDpoolP, dfP)
}
CDpoolA$year = CDpoolP$year = factor(CDpoolA$year)

# plot
pdf(file=paste0(figDir, "comDist_pool_abund_dot.pdf"), width=10, height=10)
ggplot(CDpoolA, aes(Axis1, Axis2, color=crosstype)) +
  geom_point() +
  stat_ellipse() +
  scale_color_manual(name = "Cross Type", values = colors, labels = typelabels) +
  facet_wrap(~ phylo) +
  ggtitle("Abundance-Weighted Community Distance\nAnalyses of Communities Pooled by Host Type")
dev.off()

pdf(file=paste0(figDir, "comDist_pool_pres_dot.pdf"), width=10, height=10)
ggplot(CDpoolP, aes(Axis1, Axis2, color=crosstype)) +
  geom_point() +
  stat_ellipse() +
  scale_color_manual(name = "Cross Type", values = colors, labels = typelabels) +
  facet_wrap(~ phylo) +
  ggtitle("Unweighted Community Distance\nAnalyses of Communities Pooled by Host Type")
dev.off()

# 4 panel figure
cdpp = cbind("pres", CDpoolP[CDpoolP$phylo == "ultra0" | CDpoolP$phylo == "equal0",])
cdpa = cbind("abund", CDpoolA[CDpoolA$phylo == "ultra0" | CDpoolA$phylo == "equal0",])
names(cdpp)[1] = names(cdpa)[1] = "abund"
cdp = rbind(cdpp, cdpa)

pdf(file=paste0(figDir, "comDist_pooled_pa_dot_ultraequal.pdf"), width=10, height=8)
ggplot(cdp, aes(Axis1, Axis2, color=crosstype)) +
  geom_point() +
  stat_ellipse() +
  scale_color_manual(name = "Cross Type", values = colors, labels = typelabels) +
  facet_wrap(~ abund + phylo) +
  ggtitle("Community Distance\nAnalyses of Communities Pooled by Host Type")
dev.off()



############################
# Graphs
d = comDistA$pooled[[1]]
distmat = matrix(comDistA$pooled[[1]])
hy = t(combn(labels(comDistA$pooled[[1]]), 2))
colnames(hy) = c("hy1", "hy2")
hy1 <- colsplit(string=hy[,1], pattern="20", names=c("host1","year1"))
hy2 <- colsplit(string=hy[,2], pattern="20", names=c("host2","year2"))
df = data.frame(hy, hy1, hy2, dist=as.numeric(distmat))

ggplot(df, aes(host1,dist)) + geom_boxplot()
ggplot(df, aes(host2,dist)) + geom_boxplot()

  scale_fill_manual(name = "Cross Type", values = colors, labels = typelabels) +
  facet_wrap(~ phylo) +
  ggtitle("Phylogenetic Diversity\nAnalyses of Communities on Individual Host Trees")






df$id <- apply(df[,1:4], 1, function(x) {paste(x, collapse="")})

fofo <- df$h1 == "fo" & df$h2 == "fo"
frfr <- df$h1 == "fr" & df$h2 == "fr"
nana <- df$h1 == "na" & df$h2 == "na"
fofr <- df$h1 == "fo" & df$h2 == "fr"
frfo <- df$h1 == "fr" & df$h2 == "fo"
nafr <- df$h1 == "na" & df$h2 == "fr"
frna <- df$h1 == "fr" & df$h2 == "na"
nafo <- df$h1 == "na" & df$h2 == "fo"
fona <- df$h1 == "fo" & df$h2 == "na"

df$x[fofo] <- "within hyb"
df$x[frfr] <- "within fre"
df$x[nana] <- "within nar"
df$x[fofr] <- "between hyb & fre"
df$x[frfo] <- "between hyb & fre"
df$x[nafr] <- "between nar & fre"
df$x[frna] <- "between nar & fre"
df$x[nafo] <- "between hyb & nar"
df$x[fona] <- "between hyb & nar"

####################
# Boxplot of pooled results
w <- df[df$x == "within hyb" | df$x == "within fre" | df$x == "within nar",]
b <- df[df$x == "between hyb & fre" | df$x == "between nar & fre" | df$x == "between hyb & nar",]

pdf(paste0(FIGwd,"Fig4_PB_boxplots.pdf"), width=10, height=6)
par(mfrow=c(1,2), oma=c(0,0,3,0))
par(mar=c(5,5,5,1))
boxplot(dist ~ x, data=w, 
        ylab="Community Distance", 
        axes=F, boxwex=0.8,
        ylim=c(14,33))
mtext("Within Host Types (Different Years)", side=3, line=1, cex=1)
axis(2)
axis(1, at=1:3, labels=NA)
box()
mtext(c("Fremont", "Hybrid", "Narrowleaf"), side=1, at=1:3, line=0.8) 

par(mar=c(5,1,5,5))
boxplot(dist ~ x, data=b, at=c(2,3,1),
        ylab="Community Distance", 
        axes=F, varwidth=T,
        ylim=c(14,33))
mtext("Among Host Types (All Years)", side=3, line=1, cex=1)
axis(4)
axis(1, at=1:3, labels=NA)
box()
mtext(c("Fremont & \nNarrowleaf", "Hybrid & \nFremont","Hybrid & \nNarrowleaf"), 
      side=1, at=1:3, line=1.6)

mtext("Phylobetadiversity", outer=TRUE, side=3, line=0, cex=1, font=2)
dev.off()


