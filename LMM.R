####
# Linear mixed models 
# Objective: test how tree type influences phylogenetic diversity (PD), mean PD (MPD), and net relatedness index (NRI)

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

###
# nlme worked better for building the models, and lme4 worked better for model selection on those models

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
