data(iris)
iris <- iris[seq(1, 150, by=3),]
iris.md <- distance(iris[,1:4], "mahal")

# Minimum-stress 2-dimensional nonmetric multidimensional scaling configuration
# Uses small number of separate ordinations (5) to increase speed of example.
# Use more for final analysis.
iris.nmds <- nmds(iris.md, mindim=2, maxdim=2, nits=3)
iris.nmin <- nmds.min(iris.nmds)

# Plot NMDS result with symbols denoting species
plot(iris.nmin, pch=as.numeric(iris[,5]))

# Fit vectors for the main variables to the NMDS configuration
iris.vf <- vf(iris.nmin, iris[,1:4], nperm=10)
plot(iris.vf, col="blue")
