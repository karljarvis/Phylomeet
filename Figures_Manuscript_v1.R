##########
# Calculate statistics and plot figures of  community phylogenetic analyses of arthropod communities on cottonwood hosts
# Karl Jarvis June 13, 2014
##########

library(ggplot2)
library(gridExtra)

source("~/GoogleDriveSUU/Phylomeet/Analysis_Misof_res/arcot_dataprep_v1.R")

colors = c("black","red","blue")
typelabels = c("Fremont", "Hybrid", "Narrowleaf")
types = rep(crosstypes, 4)
yrs = rep(years, each=3)

# Figure 1: Phylogeny -----------------------------------------------------

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
pdf(paste0(figDir, "Fig1_phylo_with_orders.pdf"), width=8.5, height=11)
plot(phyLongNames, cex=0.4, root.edge=T, no.margin=T)
nodelabels(1:phyLongNames$Nnode, bg="black", col="white", cex=0.5, frame="circle")
box()
dev.off()

pdf(paste0(figDir, "Fig1_phylo.pdf"), width=8.5, height=11)
plot(phyShortNames, cex=0.4, edge.width=1.5, root.edge=T, no.margin=T)
nodelabels(1:phyShortNames$Nnode, bg="black", col="white", cex=0.5, frame="circle")
box()
dev.off()

tiff(paste0(figDir, "Fig1_phylo.tiff"), width=8.5, height=11, units="in", res=300)
plot(phyShortNames, cex=0.3, edge.width=1.5, root.edge=T, no.margin=T)
nodelabels(1:phyShortNames$Nnode, bg="black", col="white", cex=0.3, frame="circle")
box()
dev.off()

# plot(phyfig, cex=0.5, edge.width=1.5, root.edge=T, no.margin=T, type="fan")

# Load data for Figures 2-3 and Table 2 --------------------------------

# putting MPD and NRI results for individual and pooled communities into long format
PDi = PDp = MNia = MNip = MNpa = MNpp = data.frame() # MPD + NRI + (Individual or Pooled) + (Abundance or Presence)
for (i in 1:length(phyNames)) 
{   
  pdi = read.csv(paste0(PDdir, "PD_indiv_", phyNames[i], ".csv"), row.names=1)
  pdp = read.csv(paste0(PDdir, "PD_pooled_", phyNames[i], ".csv"), row.names=1)
  ia = read.csv(paste0(NRIdir, "NRI_abund_indiv_", phyNames[i], ".csv"), row.names=1)
  ip = read.csv(paste0(NRIdir, "NRI_pres_indiv_", phyNames[i], ".csv"), row.names=1)
  pa = read.csv(paste0(NRIdir, "NRI_abund_pooled_", phyNames[i], ".csv"), row.names=1)
  pp = read.csv(paste0(NRIdir, "NRI_pres_pooled_", phyNames[i], ".csv"), row.names=1)
  PDi = rbind(PDi, data.frame(phylo=phyNames[i], crosstype, year, pdi))
  PDp = rbind(PDp, data.frame(phylo=phyNames[i], crosstype=types, year=yrs, pdp))
  MNia = rbind(MNia, data.frame(phylo=phyNames[i], crosstype, year, ia[,c(1:4,6:7)]))
  MNip = rbind(MNip, data.frame(phylo=phyNames[i], crosstype, year, ip[,c(1:4,6:7)]))
  MNpa = rbind(MNpa, data.frame(phylo=phyNames[i], crosstype=types, year=yrs, pa[,c(1:4,6:7)]))
  MNpp = rbind(MNpp, data.frame(phylo=phyNames[i], crosstype=types, year=yrs, pp[,c(1:4,6:7)]))
}
names(MNia)[4:9] = names(MNpa)[4:9] = c("SR","MPD","Rand_MPD","Rand_SD","mpd.obs.z","pval")
MNpa$NRI = -MNpa$mpd.obs.z
MNia$NRI = -MNia$mpd.obs.z
MNia$crosstype = factor(MNia$crosstype, levels=c("fo","fr","na"))
MNia$crosstype_fr = factor(MNia$crosstype, levels=c("fr","fo","na"))
MNia$crosstype_na = factor(MNia$crosstype, levels=c("na","fo","fr"))
MNpa$crosstype = factor(MNpa$crosstype, levels=c("fo","fr","na"))
MNpa$crosstype_fr = factor(MNpa$crosstype, levels=c("fr","fo","na"))
MNpa$crosstype_na = factor(MNpa$crosstype, levels=c("na","fo","fr"))
MNia$year = factor(MNia$year)
MNpa$year = factor(MNpa$year)

# subset to only ultra0
PDiu = PDi[PDi$phylo == "ultra0",-1]
PDpu = PDp[PDp$phylo == "ultra0",-1]
MNiau = MNia[MNia$phylo == "ultra0",-1]
MNpau = MNpa[MNpa$phylo == "ultra0",-1]

# Figure 2: PD ------------------------------------------------------------

PDy = c(PDiu$PD,PDpu$PD)
PDylims = c(min(PDy), max(PDy))

levels(PDiu$crosstype) = c("Hybrid","Fremont","Narrowleaf")
PDiu$crosstype = factor(PDiu$crosstype, levels(PDiu$crosstype)[c(2,1,3)])
levels(PDpu$crosstype) = c("Hybrid","Fremont","Narrowleaf")
PDpu$crosstype = factor(PDpu$crosstype, levels(PDpu$crosstype)[c(2,1,3)])

# Plots
PDi.plot = ggplot(PDiu, aes(year, PD)) + 
  geom_boxplot(aes(fill=crosstype)) + ylim(PDylims) +
  scale_fill_manual(name = "Tree Type", values = colors) +
  #   theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank()) +
  theme_bw() +
  ggtitle("Phylogenetic Diversity\nCommunities on Individual Trees")

PDp.plot = ggplot(PDpu, aes(year, PD, group=crosstype, color=crosstype)) + 
  geom_line() + geom_point() + ylim(PDylims) +
  scale_color_manual(name = "Tree Type", values = colors) +
  theme_bw() +
  #   theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
  #         axis.text.x = element_blank(), axis.text.y = element_blank()) +
  ggtitle("Phylogenetic Diversity\nCommunities Pooled by Tree Type")

tiff(file.path(figDir, "Fig3_PD_indiv_pooled.tif"), units="in", res=200, width=10, height=5)
grid.arrange(PDi.plot, PDp.plot, ncol=2)
dev.off()

pdf(file.path(figDir, "PD_vs_SR_plot.pdf"))
plot(PDiu$SR, PDiu$PD)
abline(lm(PD ~ SR, data=PDiu))
dev.off()

sink(file.path(figDir, "PD_vs_SR_model.txt"))
summary(lm(PD ~ SR, data=PDiu))
sink()

library(reshape2)
library(plyr)
PDu = PDue[PDue$phylo == "ultra0",]
PDm = dcast(PDu, crosstype ~ year, mean, value.var=PD)
ddply(.data=PDu, .variables=.(crosstype,year), summarize, mean=mean)

PDu$SR = as.numeric(PDu$SR)
p = data.frame(crosstype=PDu$crosstype, year=PDu$year, PD=PDu$PD)
pc = dcast(p, crosstype ~ year, mean)
rownames(pc) = pc$crosstype
pc = pc[,-1]
lapply(pc, function(x) max(x)-min(x))
pcm = unlist(lapply(pc, function(x) mean(x)))
max(pcm)-min(pcm)
apply(pc,1,function(x) max(x)-min(x))

# Figure 2: Mean Phylogenetic Distance (MPD) ------------------------------


levels(MNiau$crosstype) = c("Hybrid","Fremont","Narrowleaf")
MNiau$crosstype = factor(MNiau$crosstype, levels(MNiau$crosstype)[c(2,1,3)])
levels(MNpau$crosstype) = c("Hybrid","Fremont","Narrowleaf")
MNpau$crosstype = factor(MNpau$crosstype, levels(MNpau$crosstype)[c(2,1,3)])


MNpau$MPDu = MNpau$MPD + MNpau$Rand_SD
MNpau$MPDl = MNpau$MPD - MNpau$Rand_SD
MPDy = c(MNiau$MPD,MNpau$MPDu,MNpau$MPDl)
MPDylims = c(min(MPDy), max(MPDy))

# Plots
MPDiaplot = ggplot(MNiau, aes(year, MPD)) + 
  geom_boxplot(aes(fill=crosstype)) + ylim(MPDylims) +
  scale_fill_manual(name = "Tree Type", values = colors, labels = typelabels) +
#   theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank()) +
  theme_bw() +
  ggtitle("Mean Phylogenetic Distance\nCommunities on Individual Trees")

MPDpaplot = ggplot(MNpau, aes(year, MPD, group=crosstype, color=crosstype)) + 
  geom_line() + geom_point() + ylim(MPDylims) +
  geom_errorbar(aes(ymin=MPDl, ymax=MPDu), width=0.2) +
  scale_color_manual(name = "Tree Type", values = colors, labels = typelabels) +
  theme_bw() +
  #   theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
#         axis.text.x = element_blank(), axis.text.y = element_blank()) +
  ggtitle("Mean Phylogenetic Distance\nCommunities Pooled by Tree Type")

pdf(file=paste0(figDir, "Fig3_MPD_indiv_pooled.pdf"), width=10, height=5)
grid.arrange(MPDiaplot, MPDpaplot, ncol=2)
dev.off()

# Figure 2: Net Relatedness Index (NRI)  --------------------------------

MNpau$NRIu = MNpau$NRI + MNpau$Rand_SD
MNpau$NRIl = MNpau$NRI - MNpau$Rand_SD

NRIy = c(MNiau$NRI, MNpau$NRIu, MNpau$NRIl)
NRIylims = c(min(NRIy), max(NRIy))

# Plots
NRIiaplot = ggplot(MNiau, aes(year,NRI)) + 
  geom_boxplot(aes(fill=crosstype)) + ylim(NRIylims) +
  scale_fill_manual(name = "Tree Type", values = colors, labels = typelabels) +
#   theme(legend.position = "none") +
  theme_bw() +
  ggtitle("Net Relatedness Index\nCommunities on Individual Trees")

NRIpaplot = ggplot(MNpau, aes(year, NRI, group=crosstype, color=crosstype)) + 
  geom_line() + geom_point() + ylim(NRIylims) +
  geom_errorbar(aes(ymin = MNpau$NRIl, ymax = MNpau$NRIu), width=0.2) +
  scale_color_manual(name = "Tree Type", values = colors, labels = typelabels) +
#   theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
  theme_bw() +
  ggtitle("Net Relatedness Index\nCommunities Pooled by Tree Type")

pdf(file=paste0(figDir, "Fig2_NRI_indiv_pooled.pdf"), width=10, height=5)
grid.arrange(NRIiaplot, NRIpaplot, ncol=2)
dev.off() 

tiff(file = file.path(figDir, "Fig2_PD_MPD_NRI.tif"), width = 10, height = 12, units="in", res=200)
grid.arrange(PDi.plot, PDp.plot, MPDiaplot, MPDpaplot, NRIiaplot, NRIpaplot, ncol=2)
dev.off()

# Figure 3: Load data and NMDS --------------------------------------------

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

# load individual data and perform 2D NMDS
library(MASS)
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

# load pooled data
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


# PERMANOVA ---------------------------------------------------------------

library(vegan)

# ultra PERMANOVA
cd.ia.perm = adonis(comDistA$indiv$ultra0 ~ crosstype + year, 
                    data=CDindivA[CDindivA$phylo == "ultra0",])
cd.pa.perm = adonis(comDistA$pooled$ultra0 ~ crosstype + year, 
                    data=CDpoolA[CDpoolA$phylo == "ultra0",])
cd.ip.perm = adonis(comDistP$indiv$ultra0 ~ crosstype + year, 
                    data=CDindivP[CDindivP$phylo == "ultra0",])
cd.pp.perm = adonis(comDistP$pooled$ultra0 ~ crosstype + year, 
                    data=CDpoolP[CDpoolP$phylo == "ultra0",])

# write to file
perm.df = rbind(cbind(num="abund",com="indiv",cd.ia.perm$aov.tab[1:2,c(1,5,6)]),
                cbind(num="pres",com="indiv",cd.ip.perm$aov.tab[1:2,c(1,5,6)]),
                cbind(num="abund",com="pool",cd.pa.perm$aov.tab[1:2,c(1,5,6)]),
                cbind(num="pres",com="pool",cd.pp.perm$aov.tab[1:2,c(1,5,6)]))
write.csv(perm.df, file.path(figDir, "PERMANOVA_ultra.csv"))

sink(file.path(figDir, "CD_adonis.txt"))
print("Community Distance individual trees, abundance-weighted")
cd.ia.perm
print("Community Distance pooled trees, abundance-weighted")
cd.pa.perm
print("Community Distance individual trees, presence only")
cd.ip.perm
print("Community Distance pooled trees, presence only")
cd.pp.perm
sink()

# Figure 3: ComDist (CD) --------------------------------------------------

# plot
cdp = cbind("unweighted", CDindivP[CDindivP$phylo == "ultra0",])
cda = cbind("abundance weighted", CDindivA[CDindivA$phylo == "ultra0",])
names(cdp)[1] = names(cda)[1] = "abund"
cdi = rbind(cdp, cda)

# reorder factors
cdi$crosstype = factor(cdi$crosstype,levels(cdi$crosstype)[c(2,1,3)])
levels(cdi$crosstype) = c("Fremont","Hybrid","Narrowleaf")


CDiplot = ggplot(cdi, aes(Axis1, Axis2, color=crosstype)) +
  geom_point() + 
  stat_ellipse(type = "norm", level = 0.95) +
  scale_color_manual(name = "Tree Type", values = colors) +
  facet_wrap(~ abund) +
  theme_bw() +
  ggtitle("Ordinations of Community Distance Estimates\nCommunities on Individual Trees")

cdpp = cbind("unweighted", CDpoolP[CDpoolP$phylo == "ultra0",])
cdpa = cbind("abundance weighted", CDpoolA[CDpoolA$phylo == "ultra0",])
names(cdpp)[1] = names(cdpa)[1] = "abund"
cdp = rbind(cdpp, cdpa)

# reorder factors
cdp$crosstype = factor(cdp$crosstype,levels(cdp$crosstype)[c(2,1,3)])
levels(cdp$crosstype) = c("Fremont","Hybrid","Narrowleaf")

# plot
CDpplot = ggplot(cdp, aes(Axis1, Axis2, color=crosstype)) +
  geom_point() +
  stat_ellipse(type = "norm", level = 0.95) +
  scale_color_manual(name = "Tree Type", values = colors) +
  facet_wrap(~ abund, scales="fixed") +
  theme_bw() +
  ggtitle("Communities Pooled by Tree Type")

# cd = rbind(cbind(pool="indiv", abund = "taxa not weighted", CDindivP[CDindivP$phylo == "ultra0",]),
#            cbind(pool="pool", abund = "taxa weighted by abundance", CDindivP[CDpoolP$phylo == "ultra0",]),
#            cbind(pool="indiv", abund = "taxa not weighted", CDindivP[CDindivA$phylo == "ultra0",]),
#            cbind(pool="pool", abund = "taxa weighted by abundance", CDindivP[CDpoolA$phylo == "ultra0",]))
# 
# CD.plot = ggplot(cd, aes(Axis1, Axis2, color=crosstype)) +
#   geom_point() +
#   stat_ellipse(type = "norm", level = 0.95) +
#   scale_color_manual(name = "Cross Type", values = colors, labels = typelabels) +
#   facet_wrap( ~ abund + pool, scales="fixed") +
#   theme_bw() +
#   ggtitle("Communities Pooled by Host Type")


tiff(file=paste0(figDir, "Fig3_comDist_indiv_pooled.tif"), units="in", res=200, width=10, height=10)
grid.arrange(CDiplot, CDpplot, ncol=1)
dev.off()

# Table 2: PD, MPD, and NRI results --------------------------------------------

# abundance weighted NRI and MPD results for pooled and individual analyses
library(reshape2)
library(metap)

# Means of individual
imelt = melt(MNiau, id=c("year", "crosstype"), measure.vars = c("SR","MPD","Rand_MPD","Rand_SD","NRI","pval"))
idf = dcast(imelt, year + crosstype ~ variable, mean)
names(idf)[-1:-2] = paste0("indiv_", names(idf)[-1:-2])

ip = imelt[imelt$variable == "pval", c(1,2,4)]
fishersMethod = function(x) pchisq(-2*sum(log(x)),df=2*length(x),lower=FALSE)
ipdf = dcast(ip, year + crosstype ~ 0, fishersMethod)
idf$indiv_pval_fisher = ipdf[,3]

PDmelt = melt(PDiu, id=c("year", "crosstype"), measure.vars = c("SR","PD"))
PDdf = dcast(PDmelt, year + crosstype ~ variable, mean)
idf$indiv_PD = PDdf$PD

idf$crosstype = factor(idf$crosstype,levels(idf$crosstype)[c(2,1,3)])
idf = idf[order(idf$year, idf$crosstype),]


pdf = MNpau[,c(-1,-2)]
pdf$PD = PDpu$PD
names(pdf) = paste0("pool_",names(pdf))

tab = data.frame(idf,pdf)
write.csv(tab, file.path(figDir, "Table2_PD_MPD_NRI.csv"), row.names=F)


# Linear mixed models -----------------------------------------------------

# would be better to include in candidate set:
# type of phylo, ultra or equal, not just specific one
# specific host tree, not just type
library(nlme)
library(AICcmodavg)

PDue = PDi[PDi$phylo == "equal0" | PDi$phylo == "ultra0",]
PDue$phylo = factor(as.character(PDue$phylo))
s = read.csv(file.path(dataDir, "arcot_treenames.csv"), header = T, row.names=1)
names(s) = "treename"
s$tree = paste0(crosstype, "_", s$treename)
PDue$tree = s$tree

PDue$dummy = 1
PD.lmm.1 = lme(
  fixed = PD ~ 1,
  random = ~ 1 | dummy,
  data = PDue,
  method = "ML"
)

PDue$SR = scale(PDue$SR)
PD.lmm.list = list(
  PD.lmm.1 <- PD.lmm.1,
  PD.lmm.2 <- update(PD.lmm.1, fixed = PD ~ crosstype),
  PD.lmm.3 <- update(PD.lmm.2, random = ~ 1 | year),
  PD.lmm.4 <- update(PD.lmm.2, random = ~ 1 | phylo),
  PD.lmm.5 <- update(PD.lmm.2, random = ~ 1 | tree),
  PD.lmm.6 <- update(PD.lmm.2, random = list(~1|phylo, ~1|year)),
  PD.lmm.7 <- update(PD.lmm.2, random = list(~1|tree, ~1|year)),
  PD.lmm.8 <- update(PD.lmm.2, random = list(~1|phylo, ~1|tree)),
  PD.lmm.9 <- update(PD.lmm.2, random = list(~1|phylo, ~1|year, ~1|tree)),
  
  PD.lmm.10 <- update(PD.lmm.1, fixed = PD ~ SR),
  PD.lmm.11 <- update(PD.lmm.10, random = ~ 1 | year),
  PD.lmm.12 <- update(PD.lmm.10, random = ~ 1 | phylo),
  PD.lmm.13 <- update(PD.lmm.10, random = ~ 1 | tree),
  PD.lmm.14 <- update(PD.lmm.10, random = list(~1|phylo, ~1|year)),
  PD.lmm.15 <- update(PD.lmm.10, random = list(~1|tree, ~1|year)),
  PD.lmm.16 <- update(PD.lmm.10, random = list(~1|phylo, ~1|tree)),
  PD.lmm.17 <- update(PD.lmm.10, random = list(~1|phylo, ~1|year, ~1|tree)),
  
  PD.lmm.18 <- update(PD.lmm.1, fixed = PD ~ SR + crosstype),
  PD.lmm.19 <- update(PD.lmm.18, random = ~ 1 | year),
  PD.lmm.20 <- update(PD.lmm.18, random = ~ 1 | phylo),
  PD.lmm.21 <- update(PD.lmm.18, random = ~ 1 | tree),
  PD.lmm.22 <- update(PD.lmm.18, random = list(~1|phylo, ~1|year)),
  PD.lmm.23 <- update(PD.lmm.18, random = list(~1|tree, ~1|year)),
  PD.lmm.24 <- update(PD.lmm.18, random = list(~1|phylo, ~1|tree)),
  PD.lmm.25 <- update(PD.lmm.18, random = list(~1|phylo, ~1|year, ~1|tree))
)

# model selection table
PD.lmm.names = names(PD.lmm.list) = paste0("PD.lmm.", 1:length(PD.lmm.list))
PD.aic.df = aictab(PD.lmm.list, PD.lmm.names)
write.csv(PD.aic.df, file.path(dir, "LMM", "PD_aictab.csv"))
write.csv(summary(PD.lmm.25)$tTable, file.path(dir, "LMM", "PD_fullmod.csv"))

# MPD
MN = MNia[MNia$phylo == "equal0" | MNia$phylo == "ultra0",]
MN$phylo = factor(as.character(MN$phylo))
s = read.csv(file.path(dataDir, "arcot_treenames.csv"), header = T, row.names=1)
names(s) = "treename"
s$tree = paste0(crosstype, "_", s$treename)
MN$tree = s$tree

MN$dummy = 1
MPD.lmm.1 = lme(
  fixed = MPD ~ 1,
  random = ~ 1 | dummy,
  data = MN,
  method = "ML"
)

MN$SR = scale(MN$SR)
MPD.lmm.list = list(
  MPD.lmm.1 <- MPD.lmm.1,
  MPD.lmm.2 <- update(MPD.lmm.1, fixed = MPD ~ crosstype),
  MPD.lmm.3 <- update(MPD.lmm.2, random = ~ 1 | year),
  MPD.lmm.4 <- update(MPD.lmm.2, random = ~ 1 | phylo),
  MPD.lmm.5 <- update(MPD.lmm.2, random = ~ 1 | tree),
  MPD.lmm.6 <- update(MPD.lmm.2, random = list(~1|phylo, ~1|year)),
  MPD.lmm.7 <- update(MPD.lmm.2, random = list(~1|tree, ~1|year)),
  MPD.lmm.8 <- update(MPD.lmm.2, random = list(~1|phylo, ~1|tree)),
  MPD.lmm.9 <- update(MPD.lmm.2, random = list(~1|phylo, ~1|year, ~1|tree))
)

library(lme4)
MPD.lmer.1 = lmer(formula = MPD ~ crosstype + (1|year) + (1|phylo) + (1|tree), REML=F, data=MN)
test = lmer(MPD ~ 1 + (1|phylo), data=MN)

# model selection table
MPD.lmm.names = names(MPD.lmm.list) = paste0("MPD.lmm.", 1:length(MPD.lmm.list))
MPD.aic.df = aictab(MPD.lmm.list, MPD.lmm.names)
write.csv(MPD.aic.df, file.path(dir, "LMM", "MPD_aictab.csv"))

# 
# # NRI
MN$dummy = 1
NRI.lmm.1 = lme(
  fixed = NRI ~ 1,
  random = ~ 1 | dummy,
  data = MN,
  method = "ML"
)

MN$SR = scale(MN$SR)
NRI.lmm.list = list(
  NRI.lmm.1 <- NRI.lmm.1,
  NRI.lmm.2 <- update(NRI.lmm.1, fixed = NRI ~ crosstype),
  NRI.lmm.3 <- update(NRI.lmm.2, random = ~ 1 | year),
  NRI.lmm.4 <- update(NRI.lmm.2, random = ~ 1 | phylo),
  NRI.lmm.5 <- update(NRI.lmm.2, random = ~ 1 | tree),
  NRI.lmm.6 <- update(NRI.lmm.2, random = list(~1|phylo, ~1|year)),
  NRI.lmm.7 <- update(NRI.lmm.2, random = list(~1|tree, ~1|year)),
  NRI.lmm.8 <- update(NRI.lmm.2, random = list(~1|phylo, ~1|tree)),
  NRI.lmm.9 <- update(NRI.lmm.2, random = list(~1|phylo, ~1|year, ~1|tree))
)

# model selection table
NRI.lmm.names = names(NRI.lmm.list) = paste0("NRI.lmm.", 1:length(NRI.lmm.list))
NRI.aic.df = aictab(NRI.lmm.list, NRI.lmm.names)
write.csv(NRI.aic.df, file.path(dir, "LMM", "NRI_aictab.csv"))

# model averaging
modavg(NRI.lmm.list, "crosstype", warn = F)



##################
# LME with pooled data
# MPD
MNp = MNpa[MNpa$phylo == "equal0" | MNpa$phylo == "ultra0",]
MNp$phylo = factor(as.character(MNp$phylo))
MNp$dummy = 1

MPDp.lmm.1 = lme(
  fixed = MPD ~ 1,
  random = ~ 1 | dummy,
  data = MNp,
  method = "ML"
)

MNp$SR = scale(MNp$SR)
MPDp.lmm.list = list(
  MPDp.lmm.1 <- MPDp.lmm.1,
  MPDp.lmm.2 <- update(MPDp.lmm.1, fixed = MPD ~ crosstype),
  MPDp.lmm.3 <- update(MPDp.lmm.2, random = ~ 1 | year),
  MPDp.lmm.4 <- update(MPDp.lmm.2, random = ~ 1 | phylo),
  MPDp.lmm.5 <- update(MPDp.lmm.2, random = list(~1|phylo, ~1|year))
)

# model selection table
MPDp.lmm.names = names(MPDp.lmm.list) = paste0("MPD.lmm.", 1:length(MPDp.lmm.list))
MPDp.aic.df = aictab(MPDp.lmm.list, MPDp.lmm.names)
write.csv(MPDp.aic.df, file.path(dir, "LMM", "MPDp_aictab.csv"))


# NRI
NRIp.lmm.1 = lme(
  fixed = NRI ~ 1,
  random = ~ 1 | dummy,
  data = MNp,
  method = "ML"
)

MNp$SR = scale(MNp$SR)
NRIp.lmm.list = list(
  NRIp.lmm.1 <- NRIp.lmm.1,
  NRIp.lmm.2 <- update(NRIp.lmm.1, fixed = NRI ~ crosstype),
  NRIp.lmm.3 <- update(NRIp.lmm.2, random = ~ 1 | year),
  NRIp.lmm.4 <- update(NRIp.lmm.2, random = ~ 1 | phylo),
  NRIp.lmm.5 <- update(NRIp.lmm.2, random = list(~1|phylo, ~1|year))
)

# model selection table
NRIp.lmm.names = names(NRIp.lmm.list) = paste0("NRI.lmm.", 1:length(NRIp.lmm.list))
NRIp.aic.df = aictab(NRIp.lmm.list, NRIp.lmm.names)
write.csv(NRIp.aic.df, file.path(dir, "LMM", "NRIp_aictab.csv"))




































# 
# NRI.lmm.4.fo.REML = lme(NRI~crosstype_na, random=list(~1|phylo, ~1|year), data=MNia, method="REML")
# summary(NRI.lmm.4.na.REML)
# intervals(NRI.lmm.4.na.REML)
# 
# NRI.lmm.4.fr.REML = lme(NRI~crosstype_fr, random=list(~1|phylo, ~1|year), data=MNia, method="REML")
# summary(NRI.lmm.4.fr.REML)
# 
# MNia$crosstype = factor(MNia$crosstype, levels = c("na","fo","fr"))
# NRI.lmm.4.na.REML = lme(NRI~crosstype, random=list(~1|phylo, ~1|year), data=MNia, method="REML")
# summary(NRI.lmm.4.na.REML)




library(lme4)
# Mixed model with all phylogenies as random effects
PD.m2n = lme(fixed = PD ~ crosstype, random = list(~1 | year, ~1 | phylo), data=PDi)
PD.s2n = summary(PD.m2n)

PD.m2 = lmer(formula = PD ~ crosstype + SR + (1 | year) + (1 | phylo) + (1|tree), PDue)
PD.s2 = summary(PD.m2)
(PDci.m2 = round(confint.merMod(PD.m2), 2))
write.csv(PDci.m2, file.path(dir, "LMM", "PD_CI.csv"))
write.csv(round(PD.s2$coefficients, 2), file.path(dir, "LMM", "PD_coef.csv"))

MPD.m2 = lmer(formula = MPD ~ crosstype + (1 | year) + (1 | phylo) + (1|tree), MN)
MPD.s2 = summary(MPD.m2)
(MPDci.m2 = round(confint.merMod(MPD.m2), 2))
write.csv(MPDci.m2, file.path(dir, "LMM", "MPD_CI.csv"))
write.csv(round(MPD.s2$coefficients,2), file.path(dir, "LMM", "MPD_coef.csv"))

NRI.m2 = lmer(formula = NRI ~ crosstype + (1 | year) + (1 | phylo) + (1|tree), MN)
NRI.s2 = summary(NRI.m2)
(NRIci.m2 = round(confint.merMod(NRI.m2), 2))
write.csv(NRIci.m2, file.path(dir, "LMM", "NRI_CI.csv"))
write.csv(round(NRI.s2$coefficients,2), file.path(dir, "LMM", "NRI_coef.csv"))


# pooled analyses
PDp.m2 = lmer(formula = PD ~ crosstype + SR + (1 | year) + (1 | phylo), PDue)
(PDp.s2 = summary(PD.m2))
(PDpci.m2 = round(confint.merMod(PD.m2), 2))
write.csv(PDci.m2, file.path(dir, "LMM", "PD_CI.csv"))
write.csv(round(PD.s2$coefficients, 2), file.path(dir, "LMM", "PD_coef.csv"))

MPDp.m2 = lmer(formula = MPD ~ crosstype + (1 | year) + (1 | phylo), MNp)
MPDp.s2 = summary(MPDp.m2)
(MPDpci.m2 = round(confint.merMod(MPDp.m2), 2))
write.csv(MPDpci.m2, file.path(dir, "LMM", "MPDp_CI.csv"))
write.csv(round(MPDp.s2$coefficients,2), file.path(dir, "LMM", "MPDp_coef.csv"))

NRIp.m2 = lmer(formula = NRI ~ crosstype + (1 | year) + (1 | phylo), MNp)
NRIp.s2 = summary(NRIp.m2)
(NRIpci.m2 = round(confint.merMod(NRIp.m2), 2))
write.csv(NRIpci.m2, file.path(dir, "LMM", "NRIp_CI.csv"))
write.csv(round(NRIp.s2$coefficients,2), file.path(dir, "LMM", "NRIp_coef.csv"))



#####
# LM results
# sig differences among PD indiv results?
PDiu$year = as.numeric(as.character(PDiu$year))
PD.lm = lm(PD ~ crosstype + year, data=PDiu)
summary(PD.lm)

