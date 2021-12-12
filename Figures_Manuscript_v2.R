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

phyLongNames <- phyShortNames <- phyGenSpNames <- phyFamNames <- phyList$ultra0
inNames <- read.csv(file.path(dir, "phynames.csv"), colClasses="character")
abbr = inNames[,1]

# phylogeny with order names included
longNames <- apply(inNames[,2:5], 1, function(x) paste(x, collapse="_"))
for(i in 0:9) longNames <- gsub(i, "sp.", longNames)
longNames <- gsub("sp.sp.", "sp.", longNames)
longNames <- gsub("__", "_", longNames)
longNames <- gsub(" ", "", longNames)

# phylogeny with order names excluded	
shortNames <- apply(inNames[,3:5],1,function(x) paste(x, collapse="_"))
noGenus <- c(188,189,193,199)
for(i in noGenus) shortNames[noGenus] <- paste0(inNames[noGenus,2], "_sp.") 
for(i in 0:9) shortNames <- gsub(i, "sp.", shortNames) 
shortNames <- gsub("__", "_", shortNames)
shortNames <- gsub("sp.sp.", "sp.", shortNames)

# phylogeny with family names excluded	
GenSpNames <- apply(inNames[,4:5],1,function(x) paste(x, collapse="_"))
for (i in 0:9) GenSpNames = gsub(paste0("_",i), paste0("sp.",i), GenSpNames)

# phylogeny with only family names	
FamNames = colsplit(shortNames, "_", c("fam","gen","sp"))[,1]

phyLongNames$tip.label <- longNames[match(phyLongNames$tip.label, abbr)]
phyShortNames$tip.label <- shortNames[match(phyShortNames$tip.label, abbr)]
phyGenSpNames$tip.label <- GenSpNames[match(phyGenSpNames$tip.label, abbr)]
phyFamNames$tip.label <- FamNames[match(phyFamNames$tip.label, abbr)]

# plotting
pdf(paste0(figDir, "Fig1_phylo_with_orders.pdf"), width=8.5, height=11)
plot(phyLongNames, cex=0.4, root.edge=T, no.margin=T)
nodelabels(1:phyLongNames$Nnode, bg="black", col="white", cex=0.5, frame="circle")
box()
dev.off()

pdf(paste0(figDir, "Fig1_phylo.pdf"), width=5, height=7)
plot(phyShortNames, cex=0.4, edge.width=1.5, root.edge=T, no.margin=T)
nodelabels(1:phyShortNames$Nnode, bg="black", col="white", cex=0.5, frame="rect")
box()
dev.off()

pdf(paste0(figDir, "Fig1_phylo_Genus_species.pdf"), width=5, height=7)
plot(phyGenSpNames, cex=0.4, edge.width=1.5, root.edge=T, no.margin=T)
nodelabels(1:phyGenSpNames$Nnode, bg="black", col="white", cex=0.5, frame="rect")
box()
dev.off()

pdf(paste0(figDir, "Fig1_phylo_Family.pdf"), width=5, height=7)
plot(phyFamNames, cex=0.4, edge.width=1.5, root.edge=T, no.margin=T)
nodelabels(1:phyFamNames$Nnode, bg="black", col="white", cex=0.5, frame="rect")
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
# PDu = PDue[PDue$phylo == "ultra0",]
# PDm = dcast(PDu, crosstype ~ year, mean, value.var=PD)
# ddply(.data=PDu, .variables=.(crosstype,year), summarize, mean=mean)
# 
# PDu$SR = as.numeric(PDu$SR)
# p = data.frame(crosstype=PDu$crosstype, year=PDu$year, PD=PDu$PD)
# pc = dcast(p, crosstype ~ year, mean)
# rownames(pc) = pc$crosstype
# pc = pc[,-1]
# lapply(pc, function(x) max(x)-min(x))
# pcm = unlist(lapply(pc, function(x) mean(x)))
# max(pcm)-min(pcm)
# apply(pc,1,function(x) max(x)-min(x))

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

pdf(file=paste0(figDir, "Fig2_MPD_indiv_pooled.pdf"), width=10, height=5)
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

tiff(file = file.path(figDir, "Fig2_PD_MPD_NRI.tif"), width = 10, height = 12, units="in", res=200)
grid.arrange(PDi.plot, PDp.plot, MPDiaplot, MPDpaplot, NRIiaplot, NRIpaplot, ncol=2)
dev.off()

pdf(file = file.path(figDir, "Fig2_PD_MPD_NRI.pdf"), width = 7, height = 7)
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
  ggtitle("Ordinations of community distance estimates \namong communities on individual trees")

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
  ggtitle("Ordinations of community distance estimates \namong communities pooled by tree type")

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

pdf(file=paste0(figDir, "Fig3_comDist_indiv_pooled.pdf"), width=7, height=7)
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

library(lme4)
# Individual trees
# Mixed model with all phylogenies as random effects
# PD
PDi.ue = PDi[PDi$phylo == "equal0" | PDi$phylo == "ultra0",]
PDi.ue$phylo = factor(as.character(PDi.ue$phylo))
s = read.csv(file.path(dataDir, "arcot_treenames.csv"), header = T, row.names=1)
names(s) = "treename"
s$tree = paste0(crosstype, "_", s$treename)
PDi.ue$tree = s$tree
PDi.ue$dummy = 1
PD.mod = lmer(formula = PD ~ crosstype + SR + (1 | year) + (1 | phylo) + (1|tree), PDi.ue)
PD.coef = data.frame(coef(summary(PD.mod)))
PD.coef$pval = 2 * (1 - pnorm(abs(PD.coef$t.value)))
PD.coef = cbind(PD.coef, confint.merMod(PD.mod)[5:8,])


# MPD
MN = MNia[MNia$phylo == "equal0" | MNia$phylo == "ultra0",]
MN$phylo = factor(as.character(MN$phylo))
s = read.csv(file.path(dataDir, "arcot_treenames.csv"), header = T, row.names=1)
names(s) = "treename"
s$tree = paste0(crosstype, "_", s$treename)
MN$tree = s$tree
MN$dummy = 1

MPD.mod = lmer(formula = MPD ~ crosstype + (1 | year) + (1 | phylo) + (1|tree), MN)
MPD.coef = data.frame(coef(summary(MPD.mod)))
MPD.coef$pval = 2 * (1 - pnorm(abs(MPD.coef$t.value)))
MPD.coef = cbind(MPD.coef, confint.merMod(MPD.mod)[5:7,])

# NRI
NRI.mod = lmer(formula = NRI ~ crosstype + (1 | year) + (1 | phylo) + (1|tree), MN)
NRI.coef = data.frame(coef(summary(NRI.mod)))
NRI.coef$pval = 2 * (1 - pnorm(abs(NRI.coef$t.value)))
NRI.coef = cbind(NRI.coef, confint.merMod(NRI.mod)[5:7,])

# Table of LME of individual analyses
indiv.df = rbind(cbind(metric='PD', PD.coef),
                 cbind(metric='MPD', MPD.coef),
                 cbind(metric='NRI', NRI.coef))


# pooled analyses
# PD
PDp.ue = PDp[PDp$phylo == "equal0" | PDp$phylo == "ultra0",]
PDp.ue$phylo = factor(as.character(PDp.ue$phylo))
PDp.ue$dummy = 1

PDp.mod = lmer(formula = PD ~ crosstype + SR + (1 | year) + (1 | phylo), PDp.ue)
PDp.coef = data.frame(coef(summary(PDp.mod)))
PDp.coef$pval = 2 * (1 - pnorm(abs(PDp.coef$t.value)))
PDp.coef = cbind(PDp.coef, confint.merMod(PDp.mod)[4:7,])

# MPD
MNp.ue = MNpa[MNpa$phylo == "equal0" | MNpa$phylo == "ultra0",]
MNp.ue$phylo = factor(as.character(MNp.ue$phylo))
MNp.ue$dummy = 1

MPDp.mod = lmer(formula = MPD ~ crosstype + (1 | year) + (1 | phylo), MNp.ue)
MPDp.coef = data.frame(coef(summary(MPDp.mod)))
MPDp.coef$pval = 2 * (1 - pnorm(abs(MPDp.coef$t.value)))
MPDp.coef = cbind(MPDp.coef, confint.merMod(MPDp.mod)[4:6,])

# NRI
NRIp.mod = lmer(formula = NRI ~ crosstype + (1 | year) + (1 | phylo), MNp)
NRIp.coef = data.frame(coef(summary(NRIp.mod)))
NRIp.coef$pval = 2 * (1 - pnorm(abs(NRIp.coef$t.value)))
NRIp.coef = cbind(NRIp.coef, confint.merMod(NRIp.mod)[4:6,])

# Table of LME of pooled analyses
pooled.df = rbind(cbind(metric='PD', PDp.coef),
                 cbind(metric='MPD', MPDp.coef),
                 cbind(metric='NRI', NRIp.coef))

LMM.df = rbind(cbind(level='indiv', indiv.df),
                cbind(level='pooled', pooled.df))
rownames(LMM.df) = NULL
LMM.df$t.value = NULL
LMM.df[,c(3,4,6,7)] = lapply(LMM.df[,c(3,4,6,7)], function(x) round(x, 2))
LMM.df[,5] = round(LMM.df[,5],3)
write.csv(LMM.df, file.path(dir, "LMM", "Table_2_LMM.csv"))



# Results: community abundance --------------------------------------------

# total individuals present in dataset
totabund = sum(comList[[1]][grep("00", rownames(comList[[1]])), ],
     comList[[1]][grep("01", rownames(comList[[1]])), ],
     comList[[1]][grep("02", rownames(comList[[1]])), ],
     comList[[1]][grep("03", rownames(comList[[1]])), ])

# % of species present per pooled community
mean(apply(comList[[2]], 1, function(x) length(x[x>0]))/ncol(comList[[2]]))
# 0.319

# mean abundance per species, pooled communities
mean(apply(comList[[2]], 2, function(x) mean(x[x>0])))
# 12.1

# of species present in pooled community, min and max abundance per species
mean(apply(comList[[2]], 2, function(x) min(x[x>0])))
# 2.0
mean(apply(comList[[2]], 2, function(x) max(x[x>0])))
# 42.5



# % of species present per non-pooled community
mean(apply(comList[[1]], 2, function(x) length(x[x>0])))/ncol(comList[[2]])
# 0.075

# mean abundance per species, non-pooled communities
mean(apply(comList[[1]], 2, function(x) mean(x[x>0])))
# 2.6

# of species present in non-pooled community, min and max abundance per species
mean(apply(comList[[1]], 2, function(x) min(x[x>0])))
# 1.1
mean(apply(comList[[1]], 2, function(x) max(x[x>0])))
# 15.5
