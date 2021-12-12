
library(ape)
library(pez)
source("/Users/kjj/Projects/Phylomeet/Analysis/arcot_dataprep_v1.R")

# community
com = fullCom
nspp <- ncol(com)
nsite <- nrow(com)

# phylogeny
phy = phyPrune$indiv$ultra0

# variance-covariance matrix, i.e. branch length separating species
Vphy <- vcv(phy)
Vphy <- Vphy/(det(Vphy)^(1/nspp))

# probability of occurrence of each species in each community
com.v = unlist(as.data.frame(t(com)))
prob = log(com.v + 1)/max(log(com.v + 1))

# presence of species in each community
pres = com.v
pres[pres > 0] = 1

# factors
site <- factor(rownames(com))
species <- factor(colnames(com))
env = 

# phylogeny set as covariance structure random effect
r.intercept.spp.indep <- list(1, sp = species, covar = diag(nspp))
r.intercept.spp.phy <- list(1, sp = species, covar = Vphy)
r.slope.spp.indep <- list(env, sp = species, covar = diag(nspp))
r.slope.spp.phy <- list(env, sp = species, covar = Vphy)
r.site <- list(1, site = site, covar = diag(nsite))
rnd.effects <- list(r.intercept.spp.indep,
                    r.intercept.spp.phy, r.slope.spp.indep,
                    r.slope.spp.phy, r.site)

# run the model
# test how species presence/abundance is predicted by 
model = communityPGLMM(formula = pres ~ site, family = "binomial", sp = species, site = site, random.effects = rnd.effects, REML = TRUE, verbose = FALSE)

communityPGLMM.binary.LRT(model, re.number = 1)
communityPGLMM.binary.LRT(model, re.number = 2)

