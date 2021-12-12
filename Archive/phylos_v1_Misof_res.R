##########
# Import phylogeny, edit branch lengths, and export phylogenies with different branch lengths
# Karl Jarvis June 13, 2014
##########
library(ape)

dir = "~/Projects/Phylomeet/Analysis/"
dataDir = paste0(dir, "arcot_data/")
treeDir = paste0(dir, "Trees_Misof_res/")

# number of randomized branch lengths
nrand = 9

##########
# Import phylogeny with ultrametricized branch lengths (from Mesquite)
##########

ultra = read.nexus(paste0(treeDir,"ultra.nex"))
ultra = ladderize(ultra, FALSE)
write.tree(ultra, paste0(treeDir,"ultra0.tre"))
ultras = list(ultra)
names(ultras) = "ultra"
edgeLength = length(ultra$edge.length)

# create more trees with randomly altered branch lengths
for(i in 1:nrand)
{
  ultras[[i+1]] = ultra
  ultras[[i+1]]$edge.length = ultra$edge.length + runif(edgeLength, min=-0.99, max=0.99)
  names(ultras)[[i+1]] = paste0("urand", i)
  write.tree(ultras[[i+1]], paste0(treeDir, "ultra", i, ".tre"))
}

##########
# Change branch lenths so each edge length is set equal 
equal = ultra
equal$edge.length = rep(1, length=edgeLength)
equals = list(equal)
names(equals) = "equal"
write.tree(equal, paste0(treeDir, "equal0.tre"))

# Create trees with randomly altered branch lengths
for(i in 1:nrand)
{
  equals[[i+1]] = equal
  equals[[i+1]]$edge.length = runif(edgeLength, min=0.01, max=2)
  names(equals)[[i+1]] = paste0("erand", i)
  write.tree(equals[[i+1]], paste0(treeDir, "equal", i, ".tre"))
}

##########
# Combine all phylogenies into one list
phyloList = as.list(c(ultras, equals))
phyNames = names(phyloList)

plot(ultra, cex=0.5)
