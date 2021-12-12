# 
DataDir = '~/Projects/Phylomeet/Analysis/arcot_data/'
TreeDir = paste0(DataDir,'trees')

# Phylogeny with ultrametricized branch lengths
# read in phylogeny 
ultra = read.nexus(paste0(DataDir,'ultra.tre'))

# u1 = read.tree(paste0(DataDir,'ultra1.tre'))
# u1l = ladderize(u1)
# plot(u1, type='fan')

write.tree(ultra, paste0(DataDir,'ultra1.tre'))
ultras = list(ultra)
names(ultras) = 'ultra'
edgeLength = length(ultra$edge.length)

# number of randomized branch lengths
nrand = 9

# create nrand more trees with randomly altered branch lengths
for(i in 1:nrand)
{
  ultras[[i+1]] = ultra
  ultras[[i+1]]$edge.length = ultra$edge.length + runif(edgeLength, min=-0.99, max=0.99)
  names(ultras)[[i+1]] = paste0('urand', i)
  write.tree(ultras[[i+1]], paste0(DataDir, 'ultra', i+1, '.tre'))
}

# Phylogeny with each edge length set equal 
# read in phylogeny 
equal = ultra
equal$edge.length = rep(1, length=edgeLength)
equals = list(equal)
names(equals) = 'equal'
write.tree(equal, paste0(DataDir, 'equal1.tre'))

# create nrand more trees with randomly altered branch lengths
for(i in 1:nrand)
{
  equals[[i+1]] = equal
  equals[[i+1]]$edge.length = runif(edgeLength, min=0.01, max=2)
  names(equals)[[i+1]] = paste0('erand', i)
  write.tree(equals[[i+1]], paste0(DataDir, 'equal', i+1, '.tre'))
}

# Combine all phylogenies into one list
phyloList = as.list(c(ultras, equals))
phyNames = names(phylos)