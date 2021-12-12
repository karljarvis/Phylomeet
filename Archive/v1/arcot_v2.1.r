#############################################################
### Community Phylogenetic Analyses  
###	Arthropod communities on cottonwood hosts
###	Dataset from Wimp et al. 2004
###
### Code by Karl Jarvis
#############################################################

# Gery notes for next steps:
# (3) pooling each host across years
# (2) NRI with new topology
# (1) phylobeta diversity


# load package
# install.packages("picante", dependencies = TRUE, repos = "http://R.research.att.com/")
library(picante)

# Picante manual
# vignette("picante-intro")

# Set directory for input data
DATAwd <- "~/Projects/Phylomeet/Analysis/arcot_data/"

#############################################################
# Load in communities, take subsets and format for analysis #
#############################################################

# load community data
fullcom <- read.table(paste(DATAwd, "arcot.txt", sep=""))
host <- rep(c(rep("fremont",10), rep("f1",10), rep("backcross",10), rep("narrowleaf",10)),4)
genhost <- rep(c(rep("fr",10), rep("fo",10), rep("na",20)),4)
year <- c(rep(2000,40),rep(2001,40),rep(2002,40),rep(2003,40))
com <- cbind(year, genhost, fullcom)

# pool data
	# byHost: sum of occcurrences in all years by host
	sumFr <- apply(com[com$genhost=="fr",3:ncol(com)],2,sum)
	sumFo <- apply(com[com$genhost=="fo",3:ncol(com)],2,sum)
	sumNa <- apply(com[com$genhost=="na",3:ncol(com)],2,sum)
	byHost <- rbind(sumFr,sumFo,sumNa)
	
	# byYear: sum of occcurrences by year on all hosts
	sum2000 <- apply(com[com$year==2000,3:ncol(com)],2,sum)
	sum2001 <- apply(com[com$year==2001,3:ncol(com)],2,sum)
	sum2002 <- apply(com[com$year==2002,3:ncol(com)],2,sum)
	sum2003 <- apply(com[com$year==2003,3:ncol(com)],2,sum)
	byYear <- rbind(sum2000,sum2001,sum2002,sum2003)
	
	# byHY: sum of occurrences by year and by host
	byHYmaker <- function(com) 
	{ 	
		fr2000 <- apply(com[com$year==2000 & com$genhost=="fr",3:ncol(com)],2,sum)
		fo2000 <- apply(com[com$year==2000 & com$genhost=="fo",3:ncol(com)],2,sum)
		na2000 <- apply(com[com$year==2000 & com$genhost=="na",3:ncol(com)],2,sum)
		fr2001 <- apply(com[com$year==2001 & com$genhost=="fr",3:ncol(com)],2,sum)
		fo2001 <- apply(com[com$year==2001 & com$genhost=="fo",3:ncol(com)],2,sum)
		na2001 <- apply(com[com$year==2001 & com$genhost=="na",3:ncol(com)],2,sum)
		fr2002 <- apply(com[com$year==2002 & com$genhost=="fr",3:ncol(com)],2,sum)
		fo2002 <- apply(com[com$year==2002 & com$genhost=="fo",3:ncol(com)],2,sum)
		na2002 <- apply(com[com$year==2002 & com$genhost=="na",3:ncol(com)],2,sum)
		fr2003 <- apply(com[com$year==2003 & com$genhost=="fr",3:ncol(com)],2,sum)
		fo2003 <- apply(com[com$year==2003 & com$genhost=="fo",3:ncol(com)],2,sum)
		na2003 <- apply(com[com$year==2003 & com$genhost=="na",3:ncol(com)],2,sum)
		byHY <- rbind(fr2000,fo2000,na2000,fr2001,fo2001,na2001,fr2002,fo2002,na2002,fr2003,fo2003,na2003)
		byHY
	}
	
	byHmaker <- function(com)
	{
		fr <- apply(com[com$genhost=="fr",3:ncol(com)],2,sum)
		fo <- apply(com[com$genhost=="fo",3:ncol(com)],2,sum)
		na <- apply(com[com$genhost=="na",3:ncol(com)],2,sum)
		byH <- rbind(fr, fo, na)
		byH		
	}
	
	byHY <- byHYmaker(com=com)

	# sumAll: sum of all occurrences in all years on all hosts
	sumAll <- apply(com[,4:ncol(com)],2,sum)
	
	# HY: pooled by Host and Year
	HY <- byHYmaker(com=com)
	
	# H: pooled by Host only
	H <- byHmaker(com=com)	

# subset data based on commonness
fullsumAll <- rbind(fullcom,sumAll)

	# sp2: exclude species that have only one occurrence
	sp2 <- fullsumAll[1:160,sumAll >= 2] 

	# sp2HY: com2sp, pooled by host and year
	sp2HY <- byHYmaker(com=cbind(year, genhost, sp2))

	# sp3: exclude species that have only one occurrence
	sp3 <- fullsumAll[1:160,sumAll >= 3] 

	# sp3HY: com3sp, pooled by host and year
	sp3HY <- byHYmaker(com=cbind(year, genhost, sp3))

	# yr2: all except species that are only present in one year
	y <- byYear
	for(i in 1:ncol(y)) 
	{ 
		for(j in 1:nrow(y)) 
		{ 
			if (byYear[j,i] > 0) {y[j,i] <- 1} 
		}	
	}
	nyear <- apply(y,2,sum)
	fullnyear <- rbind(fullcom,nyear)
	yr2 <- fullnyear[1:160,nyear >= 2]

	# yr2HY: y2, pooled by host and year
	yr2HY <- byHYmaker(com=cbind(year, genhost, yr2))
		
	# yr2sp2: intersection of sp2 and yr2
	sp2year <- rbind(sp2,nyear)
	yr2sp2 <- sp2year[1:160,sp2year[161,] >= 2]

	# yr2sp2HY: yr2sp2, pooled by host and year
	yr2sp2HY <- byHYmaker(com=cbind(year, genhost, yr2sp2))

	# yr2sp3: intersection of sp3 and yr2
	sp3year <- rbind(sp3,nyear)
	yr2sp3 <- sp3year[1:160,sp3year[161,] >= 2]

	# yr2sp3HY: yr2sp3, pooled by host and year
	yr2sp3HY <- byHYmaker(com=cbind(year, genhost, yr2sp3))
	
	# yr2sp3H: yr2sp3, pooled by host
	yr2sp3H <- byHmaker(com=cbind(year, genhost, yr2sp3))
	

# list of full community and subsets	
	 # coms <- lapply(list(fullcom, HY, sp2, sp3, yr2, yr2sp2, yr2sp3, sp2HY, sp3HY, yr2HY, yr2sp2HY, yr2sp3HY), as.matrix)
	 # names(coms) <- c("fullcom", "HY", "sp2", "sp3", "yr2", "yr2sp2", "yr2sp3", "sp2HY", "sp3HY", "yr2HY", "yr2sp2HY", "yr2sp3HY")
	 

# list of the full community and the smallest subset	
	 # coms <- lapply(list(fullcom, HY, yr2sp3, yr2sp3HY), as.matrix)
	 # names(coms) <- c("fullcom", "HY", "yr2sp3", "yr2sp3HY")

# list of the smallest subset	
	 coms <- lapply(list(yr2sp3, yr2sp3HY, yr2sp3H), as.matrix)
	 names(coms) <- c("indiv", "pooled", "host")
	 

#########################################################################
# Load in phylogenies, prune, and create distance matrices for analysis #
#########################################################################

# Phylogeny with ultrametricized branch lengths
	# read in phylogeny 
	ultra <- read.nexus(paste(DATAwd,"phylos/","ultra",".tre",sep=""))
	ultras <- list(ultra)
	names(ultras) <- "ultra"
	lel <- length(ultra$edge.length)

	# number of randomized branch lengths
	nrand = 9

	# create nrand more trees with randomly altered branch lengths
	for(i in 1:nrand)
	{
		ultras[[i+1]] <- ultra
		ultras[[i+1]]$edge.length <- ultra$edge.length + runif(lel, min=-0.99, max=0.99)
		names(ultras)[[i+1]] <- paste0("urand", i)
	}

# Phylogeny with each edge length set equal 
	# read in phylogeny 
	equal <- ultra
	equal$edge.length <- rep(1, length=lel)
	equals <- list(equal)
	names(equals) <- "equal"

	# create nrand more trees with randomly altered branch lengths
	for(i in 1:nrand)
	{
		equals[[i+1]] <- equal
		equals[[i+1]]$edge.length <- runif(lel, min=0.01, max=2)
		names(equals)[[i+1]] <- paste0("erand", i)
	}
	
# Combine all phylogenies into one list
	phylos <- as.list(c(ultras, equals))
	phyNames <- names(phylos)

# Plot phylogenies
	dev.new()
	par(mfrow=c(4,2))
	for(i in 1:length(phylos))
	{
		plot(phylos[[i]], cex=0.5, show.tip.label=F)
	}
	
# Prune phylogenies by community
	phy <- vector("list",length=length(coms))
	names(phy) <- names(coms)
	for (c in 1:length(coms)) 
	{	
		phy[[c]] <- vector("list", length=length(phyNames))
		names(phy[[c]]) <- phyNames
		for (p in 1:length(phylos)) 
		{ 	
			phy[[c]][[p]] <- prune.sample(coms[[c]], phylos[[p]]) 
		}
	}

	plot(phy$pooled$ultra, cex=0.6)		

# Create distance matrix for each topology
	phydist <- vector("list",length=length(coms))
	names(phydist) <- names(coms)
	for (c in 1:length(coms)) 
	{	phydist[[c]] <- vector("list", length=length(phyNames))
		names(phydist[[c]]) <- phyNames
		for (p in 1:length(phylos)) 
		{ 	phydist[[c]][[p]] <- cophenetic(phy[[c]][[p]]) 	}
	}

# Plotting phylogenies
	# First five sets of phylogenies to get overview
	# par(mfrow=c(4,3))
	# for (c in 1:9) 
	# {	
		# plot(phy[[c]][[1]], show.tip.label=F) 
	# }

	# # large versions of full phylogenies
	# par(mfrow=c(1,2))
	# plot(phy$fullcom$ultra, cex=0.5)
	# plot(phy$fullcom$equal, cex=0.5)

		
#################################
# Faith's PD
#################################
# Set folder for output
	PDwd <- "~/Projects/Phylomeet/Analysis/Results_PD/"
	setwd(PDwd)

# Make a list to contain PD results
	PD <- vector("list",length=length(coms))
	names(PD) <- names(coms)

# Find PD for all communities and branch lengths
	for (c in 1:length(coms)) 
	{	PD[[c]] <- vector("list",length=length(phyNames))
		for (p in 1:length(phylos)) 
		{ 	
			PD[[c]][[p]] <- pd(coms[[c]], phy[[c]][[p]], include.root=TRUE)
			write.csv(PD[[c]][[p]], paste0(PDwd, "PD_", names(coms)[c], "_", phyNames[p], ".csv"))	# write PD to file
		}
	}
	
# Read in PD from file
	PD <- vector("list",length=length(coms))
	names(PD) <- names(coms)
	for (c in 1:length(coms)) 
	{	
		# create list to hold each set of results within each community 
		PD[[c]] <- vector("list",length=length(phyNames))
		names(PD[[c]]) <- phyNames
		for (p in 1:length(phylos)) 
		{ 	
			# read in results
			PD[[c]][[p]] <- read.csv(paste0(PDwd, "PD_", names(coms)[c], "_", phyNames[p], ".csv"), row.names=1)	
		}
	}

# Analyze PD
years <- 2000:2003
hosts <- c("fr","fo","na")
HYnames <- c("fr2000","fo2000","na2000","fr2001","fo2001","na2001","fr2002","fo2002","na2002","fr2003","fo2003","na2003")
dir.create("Figures/")
setwd(paste0(PDwd,"Figures/"))
spaces <- c(0.2,0.2,0.2,1,0.2,0.2,1,0.2,0.2,1,0.2,0.2)

# Individual
	PDindiv <- vector("list",length=length(phyNames))
	names(PDindiv) <- phyNames
	for (p in 1:length(phylos)) 
	{ 	
		PDmean <- PDse <- c()
		PDhy <- cbind(genhost, year, PD[[1]][[p]])
		for(y in 1:length(years))
		{
			for(h in 1:length(hosts))
			{
				PDmean <- c(PDmean, mean(PDhy[PDhy[,1] == hosts[h] & PDhy[,2] == years[y],4]))
				PDse <- c(PDse, sd(PDhy[PDhy[,1] == hosts[h] & PDhy[,2] == years[y],4])/sqrt(length(PDhy[PDhy[,1] == hosts[h] & PDhy[,2] == years[y],4])))
			}
		}
		PDindiv[[p]] <- cbind(PDmean,PDse)
		rownames(PDindiv[[p]]) <- HYnames
	}


	pdf( file=paste0( "PD_", names(PD[1]), ".pdf" ), width=10, height=10 )
	par(mfrow=c(4,5))

	for (p in 1:10) 
	{ 	
		barplot(PDindiv[[p]][,1], col = c("black", "red", "blue"), space=spaces, main="Ultra")
	}

	for (p in 11:20) 
	{ 	
		barplot(PDindiv[[p]][,1], col = c("black", "red", "blue"), space=spaces, main="Equal")
	}		
	dev.off()




# Pooled and by host
for(c in 2:length(coms))
{
	pdf( file=paste0( "PD_", names(PD[c]), ".pdf" ), width=10, height=10 )
	par(mfrow=c(4,5))
	for (p in 1:10) 
	{ 	
		barplot(PD[[c]][[p]]$PD, space = c(0.2,0.2,0.2,1,0.2,0.2,1,0.2,0.2,1,0.2,0.2), ylim = c(0,140), col = c("black", "red", "blue"), main="Ultra")
	}

	for (p in 11:20) 
	{ 	
		barplot(PD[[c]][[p]]$PD, space = c(0.2,0.2,0.2,1,0.2,0.2,1,0.2,0.2,1,0.2,0.2), ylim = c(0,140), col = c("black", "red", "blue"), main="Equal")
	}		
	dev.off()
}




for( c in 1:length(ComDist) )
{

	par( mfrow=c(2,5) )
	for( p in 1:10 )
	{
		fit <- isoMDS(PD[[c]][[p]], k=2)
		x <- fit$points[,1]
		y <- fit$points[,2]
		plot( x, y, col=cols, pch=syms, cex=2, main = paste( names( PD[c]), names( PD[[c]][p] ), sep=" " ) )
	}
	dev.off()
}





#################################
# NRI 
#################################
# NRI - abundance weighted
# takes 10 seconds at 99 runs for each branch length setting

# Set folder for output
NRIwd <- "~/Projects/Phylomeet/Analysis/Results_NRI/"

# Make a list to contain NRI results
NRI <- vector("list",length=length(coms))

# Find NRI for all communities and branch lengths
for (c in 1:length(coms)) 
{	NRI[[c]] <- vector("list",length=length(phyNames))
	for (p in 1:length(phylos)) 
	{ 	NRI[[c]][[p]] <- ses.mpd(coms[[c]], phydist[[c]][[p]], null.model = "richness", abundance.weighted = TRUE, runs = 999)

		# write results to file
		write.csv(NRI[[c]][[p]], paste0(NRIwd, "NRI_", names(coms)[c],"_", phyNames[p], ".csv"))
		
		# read in results from file
		NRI[[c]][[p]] <- read.csv(paste0(NRIwd, "NRI_", names(coms)[c],"_", phyNames[p], ".csv"))
		
	}
}


#################################
# Phylogenetic Beta Diversity
#################################
# set folder for output
ComDistwd <- "~/Projects/Phylomeet/Analysis/Results_ComDist/"

# make a list to contain comdist results
ComDist <- vector("list",length=length(coms))

# find comdist for all communities and branch lengths
	# loop through community subsets
	for (c in 1:length(coms)) 
	{	
		# list to contain results for phylogeny
		ComDist[[c]] <- vector("list",length=length(phyNames)) 
		names(ComDist[[c]]) <- phyNames
		
		# loop through phylogenies within each community subset
		for (p in 1:length(phylos)) 
			{ 	
			ComDist[[c]][[p]] <- comdist(coms[[c]], phydist[[c]][[p]], abundance.weighted=T)
			write.csv( as.matrix( ComDist[[c]][[p]]), paste0( ComDistwd, "ComDist_", names(coms)[c], "_", phyNames[p], ".csv"))
			}

		# alternatively, import pre-run results
		for (p in 1:length(phylos)) 
		{ 	
			ComDist[[c]][[p]] <- dist( read.csv( paste0( ComDistwd, "ComDist_", names(coms)[c], "_", phyNames[p], ".csv" ), row.names=1 ) ) 
		}
	}
	

# Nonmetric MDS 2D
	library(MASS)
	library(car)
	type <- rep(c("Fre","F1","Nar"),4)		# names for labels
	cols <- rep(c("black","red","blue"),4)					# colors for labels
	syms <- rep(c(15,16,17),4)
	names(ComDist) <- names(coms)


# plot into pdf files: ultra
for( c in 1:length(ComDist) )
{
	pdf( file=paste0( names( ComDist[c]), names( ComDist[[c]][1] ), ".pdf" ), width=11, height=5.5 )
	par( mfrow=c(2,5) )
	for( p in 1:10 )
	{
		fit <- isoMDS(ComDist[[c]][[p]], k=2)
		x <- fit$points[,1]
		y <- fit$points[,2]
		plot( x, y, col=cols, pch=syms, cex=2, main = paste( names( ComDist[c]), names( ComDist[[c]][p] ), sep=" " ) )
	}
	dev.off()
}

# plot into pdf files: equal
for( c in 1:length(ComDist) )
{
	pdf( file=paste0( names( ComDist[c]), names( ComDist[[c]][11] ), ".pdf" ), width=11, height=5.5 )
	par( mfrow=c(2,5) )
	for( p in 11:20 )
	{
		fit <- isoMDS(ComDist[[c]][[p]], k=2)
		x <- fit$points[,1]
		y <- fit$points[,2]
		plot( x,y, col=cols, pch=syms, cex=2, main = paste( names( ComDist[c]), names( ComDist[[c]][p] ), sep=" " ) )
	}
	dev.off()
}


	




	# # plot to quartz windows
	# for( c in 1:length(ComDist) )
	# {
		# dev.new()
		# par( mfrow=c(2,5) )
		# for( p in 1:10 )
		# {
			# fit <- isoMDS(ComDist[[c]][[p]], k=2)
			# x <- fit$points[,1]
			# y <- fit$points[,2]
			# plot( x,y, col=cols, pch=16, cex=2, main = paste( names( ComDist[c]), names( ComDist[[c]][p] ), sep=" " ) )
		# }
	
		# dev.new()
		# par(mfrow=c(2,5))
		# for( p in 11:20 )
		# {
			# fit <- isoMDS(ComDist[[c]][[p]], k=2)
			# x <- fit$points[,1]
			# y <- fit$points[,2]
			# plot( x,y, col=cols, pch=16, cex=2, main = paste( names( ComDist[c]), names( ComDist[[c]][p] ), sep=" " ) )
		# }
	# }






	# par(mfrow=c(2,4))
	# for( i in 1:length(ComDist[[2]]))
	# {
		# fit <- isoMDS(ComDist[[2]][[i]], k=2)
		# x <- fit$points[,1]
		# y <- fit$points[,2]
		# plot( x,y, col=cols, main = paste0( "Pooled Hosts ", names( ComDist[[2]][i] ) ) )
	# }

			# scatterplot(x,y, groups=type, by.group=T, reg=F, ellipse=F, smooth=F, level=0.80, xlab="Axis 1", ylab="Axis 2", main="Nonmetric MDS of Phylogenetic Community Distance of \nArthopod Communities on Fremont, F1, and Narrowleaf", legend.plot=FALSE, pch=16, cex=0)
			# text(x, y, labels = rownames(byHY), cex=1, col=cols)
			# legend("topleft",legend=c("Fre","Hyb","Nar"),col=c(2,1,3), pch=19)
	
	
# 3D Analysis 
	library(MASS)
	fit <- isoMDS(ComDist$yr2sp3HY$ultra, k=3) 
	x <- fit$points[,1]
	y <- fit$points[,2]
	z <- fit$points[,3]
	cols <- rep(c(2,1,3),4)
	library(rgl)
	plot3d(x,y,z, type="s", col=cols, main="Phylobetadiversity 3D NMDS")
	
	# Nonmetric MDS - individual
	fit <- isoMDS(ComDist$com3spHY$ultra, k=2)

	# 2x2 panel by year, no ellipses
	cols <- rep(c(rep(1,10),rep(2,10),rep(3,20)),4)
	x <- fit$points[,1]
	y <- fit$points[,2]
	par(mfrow=c(2,2))
	for(i in 0:3) 
	{	plot(x[40*i+(1:40)], y[40*i+(1:40)], xlab="Axis 1", ylab="Axis 2", main="Nonmetric MDS", type="p", pch=3, col=cols[1:40])	
	}

	# 2x2 panel by year with ellipses
	require(car)
	par(mfrow=c(2,2))
	for(i in 0:3) 
	{
	dataEllipse(x[40*i+(1:10)],y[40*i+(1:10)],levels=0.90, xlim=c(-7,7), ylim=c(-7,7), col="blue", xlab="Axis 1", ylab="Axis 2", main=paste("Phylobetadiversity",2000+i,sep=" "))
	dataEllipse(x[40*i+(11:20)],y[40*i+(11:20)],levels=0.90, xlim=c(-7,7), ylim=c(-7,7), add=TRUE, col="red")
	dataEllipse(x[40*i+(21:40)],y[40*i+(21:40)],levels=0.90, xlim=c(-7,7), ylim=c(-7,7), add=TRUE, col="green")
	legend("topleft",legend=c("Fre","Hyb","Nar"),col=c(1:3), pch=c(1,1,1))
	}

	# all in one plot
	grpnames <- c("Fre 2000","Fre 2001","Fre 2002","Fre 2003", "Hyb 2000", "Hyb 2001", "Hyb 2002", "Hyb 2003", "Nar 2000", "Nar 2001", "Nar 2002", "Nar 2003")
	cols <- c("blue1","blue2","blue3","blue4", "red1","red2","red3","red4","green1","green2","green3","green4")
	plot(x[1:10], y[1:10], xlim=c(-8,8), ylim=c(-8,8), col="blue1", xlab="Axis 1", ylab="Axis 2", main="Phylobetadiversity of Host Tree Communities", pch=15)
	for(i in 0:3)
	{	
	points(x[40*i+1:10], y[40*i+1:10], col=paste("blue",i+1,sep=""), pch=15)
	points(x[40*i+11:20],y[40*i+11:20], col=paste("red",i+1,sep=""), pch=17)
	points(x[40*i+21:40],y[40*i+21:40], col=paste("green",i+1,sep=""), pch=19)
	legend("topleft", legend=grpnames, col=cols, pch=c(rep(15,4),rep(17,4),rep(19,4)))
	}	

# 3D analysis - individual host trees
	library(MASS)
	fit <- isoMDS(comdist.equal, k=3)
	x <- fit$points[,1]
	y <- fit$points[,2]
	z <- fit$points[,3]
	
	type <- rep(c("Fre","F1","Nar","Nar"),4)
	cols <- rep(c(2,1,3,3),4)
	rownames(byTY)

	library(rgl)
	plot3d(x,y,z, type="s", size=1.5, col=cols, main="Phylobetadiversity 3D NMDS")


	# MANOVA
	type <- rep(c(rep("Fre",10),rep("Hyb",10),rep("Nar",20)),4)
	fit <- manova(cbind(x,y) ~ factor(type))
	(sumfit<-summary(fit,test="Wilks"))
	
	fit2000 <- manova(cbind(x[1:40],y[1:40]) ~ factor(type[1:40]))
	(sumfit2000<-summary(fit2000,test="Wilks"))

	fit2001 <- manova(cbind(x[41:80],y[41:80]) ~ factor(type[1:40]))
	(sumfit2001<-summary(fit2001,test="Wilks"))

	fit2002 <- manova(cbind(x[81:120],y[81:120]) ~ factor(type[1:40]))
	(sumfit2002<-summary(fit2002,test="Wilks"))

	fit2003 <- manova(cbind(x[121:160],y[121:160]) ~ factor(type[1:40]))
	(sumfit2003<-summary(fit2003,test="Wilks"))

	xy <- cbind(x,y)
	(m2000 <- colMeans(xy[1:40,]))
	(m2001 <- colMeans(xy[41:80,]))
	(m2002 <- colMeans(xy[81:120,]))
	(m2003 <- colMeans(xy[121:160,]))
	(w <- sumfit$SS$Residuals)
	
	# MVN - 2000-2003
	Skres<-fit$residuals
	S<-cov(Skres)
	MHX<-mahalanobis(Skres, center=rep(0,2), cov=S)
	Qtiles<-qchisq(ppoints(length(type)),df=2)
	qqplot(Qtiles, MHX,main = expression("Phylobetadiversity Data: Q-Q plot of " * ~D^2 * " for residuals vs. quantiles of" * ~ chi[2]^2), xlab=expression(chi[2]^2 * " quantiles"), 
	ylab="Mahalanobis Distance for model residuals")
	abline(0, 1)
	source('~/Documents/Classes/STA572_Multivar_Stats/CvM.R', chdir = TRUE)
	CvM.test(MHX,df=2,nsim=10000)

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
