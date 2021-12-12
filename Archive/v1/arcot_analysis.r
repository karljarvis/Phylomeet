#########################################################################
# Community Phylogenetic Analysis of arthropod communities on cottonwood hosts
#########################################################################

#########################################################################
# Set folder locations

# Set folder for output of PD results
PDDir <- "~/Projects/Phylomeet/Analysis/Results_PD/"

# Set folder for output of NRI results
NRIDir <- "~/Projects/Phylomeet/Analysis/Results_NRI/"

# Set folder for output of ComDist results
comDistDir <- "~/Projects/Phylomeet/Analysis/Results_ComDist/"

# Set directory for figures
figDir <- "~/Projects/Phylomeet/Analysis/Figures/"

library(reshape2)

#########################################################################
# Faith's Phylogenetic Diversity (PD)
#########################################################################
# Make a list to contain PD results
	PD <- vector("list",length=length(coms))
	names(PD) <- names(coms)

#########################################################################
# Calculate PD for all communities and branch lengths
	for (c in 1:length(coms)) 
	{	PD[[c]] <- vector("list",length=length(phyNames))
		for (p in 1:length(phylos)) 
		{ 	
			PD[[c]][[p]] <- pd(coms[[c]], phy[[c]][[p]], include.root=TRUE)
			write.csv(PD[[c]][[p]], paste0(PDDir, "PD_", names(coms)[c], "_", phyNames[p], ".csv"))	# write PD to file
		}
	}
	
#########################################################################
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
			PD[[c]][[p]] <- read.csv(paste0(PDDir, "PD_", names(coms)[c], "_", phyNames[p], ".csv"), row.names=1)	
		}
	}


#########################################################################
# Summarize and plot PD results

#########################################################################
# PD: Individual
	# PD: Summary statistics
source('~/Projects/Phylomeet/Analysis/IndivSummaryMaker.R', chdir = TRUE)
PD_full_indiv <- IndivSummaryMaker(result=PD, community=1, column=3)
PD_common_indiv <- IndivSummaryMaker(result=PD, community=2, column=3)

	# Plot full
	pdf( file=paste0( figDir, "PD_", names(PD[1]), ".pdf" ), width=10, height=10 )
	par(mfrow=c(4,5))
	for (p in 1:20) 
	{ 	barplot(PD_full_indiv[[p]][,3], col = c("black", "red", "blue"), space=spaces, main=names(PD_full_indiv)[[p]], ylim=c(0,120))	}
	dev.off()
	
	# Plot common
	pdf( file=paste0( figDir, "PD_", names(PD[2]), ".pdf" ), width=10, height=10 )
	par(mfrow=c(4,5))
	for (p in 1:20) 
	{ 	
		barplot(PD_common_indiv[[p]][,3], col = c("black", "red", "blue"), space=spaces, main=names(PD_common_indiv)[[p]], ylim=c(0,120))
	}
	dev.off()


# single representative plot
PDplot <- data.frame(crosstype,year=substr(year,3,4),PD[[1]][[1]])
ctnew <- gsub("fo","hy",crosstype)
boxplot(PD ~ ctnew + year, PDplot, at=c(1,2,3,5,6,7,9,10,11,13,14,15), col=c("blue","green3","yellow"))
boxplot(SR ~ ctnew + year, PDplot, at=c(1,2,3,5,6,7,9,10,11,13,14,15) ,col=c("blue","green3","yellow"))



#########################################################################
# PD: Pooled by crosstype within years
	# PD: full community data pooled
	pdf( file=paste0( figDir, "PD_", names(PD[3]), ".pdf" ), width=10, height=10 )
	par(mfrow=c(4,5))
	for (p in 1:20) 
	{ 	
		barplot(PD[[3]][[p]][,1], space = c(0.2,0.2,0.2,1,0.2,0.2,1,0.2,0.2,1,0.2,0.2), col = c("black", "red", "blue"), main=names(PD$full_pooled)[[p]], ylim=c(0,250))
	}
	dev.off()

	# PD: common community data pooled
	pdf( file=paste0( figDir, "PD_", names(PD[4]), ".pdf" ), width=10, height=10 )
	par(mfrow=c(4,5))
	for (p in 1:20) 
	{ 	
		barplot(PD[[4]][[p]][,1], space = c(0.2,0.2,0.2,1,0.2,0.2,1,0.2,0.2,1,0.2,0.2), col = c("black", "red", "blue"), main=names(PD$common_pooled)[[p]], ylim=c(0,250))
	}
	dev.off()

barplot(PD[[3]][[1]][,1], space = c(0.2,0.2,0.2,1,0.2,0.2,1,0.2,0.2,1,0.2,0.2), col=c("blue","green3","yellow"))
barplot(PD[[3]][[1]][,2], space = c(0.2,0.2,0.2,1,0.2,0.2,1,0.2,0.2,1,0.2,0.2), col=c("blue","green3","yellow"))

#########################################################################
# PD: Pooled by crosstype across all years
	# PD: full community data pooled by crosstype
	pdf( file=paste0( figDir, "PD_", names(PD[5]), ".pdf" ), width=10, height=10 )
	par(mfrow=c(4,5))
	for (p in 1:20) 
	{ 	
		barplot(PD[[5]][[p]][,1], col = c("black", "red", "blue"), main=names(PD$full_crosstype)[[p]], ylim=c(0,300))
	}
	dev.off()

	# PD: common community data pooled by crosstype
	pdf( file=paste0( figDir, "PD_", names(PD[6]), ".pdf" ), width=10, height=10 )
	par(mfrow=c(4,5))
	for (p in 1:20) 
	{ 	
		barplot(PD[[6]][[p]][,1], col = c("black", "red", "blue"), main=names(PD$full_crosstype)[[p]], ylim=c(0,300))
	}
	dev.off()




#########################################################################
# Net Relatedness Index (NRI) 
#########################################################################
# NRI - abundance weighted


#########################################################################
# Make a list to contain NRI results
NRI_abund <- NRI_pres <- vector("list",length=length(coms))

# Estimate NRI for all communities and branch lengths
	for (c in 1:length(coms)) 
	{	NRI_abund[[c]] <- vector("list",length=length(phyNames))
		NRI_pres[[c]] <- vector("list",length=length(phyNames))
		for (p in 1:length(phylos)) 
		{ 	
			NRI_abund[[c]][[p]] <- ses.mpd(coms[[c]], phydist[[c]][[p]], null.model = "richness", abundance.weighted = TRUE, runs = 9)
			NRI_pres[[c]][[p]] <- ses.mpd(coms[[c]], phydist[[c]][[p]], null.model = "richness", abundance.weighted = FALSE, runs = 9)
	
			# write results to file
			write.csv(NRI_abund[[c]][[p]], paste0(NRIDir, "NRI_abund_", names(coms)[c],"_", phyNames[p], ".csv"))
			write.csv(NRI_pres[[c]][[p]], paste0(NRIDir, "NRI_pres_", names(coms)[c],"_", phyNames[p], ".csv"))
		}
	}

#########################################################################
# Read in NRI from file
	NRI_abund <- NRI_pres <- vector("list",length=length(coms))
	names(NRI_abund) <- names(NRI_pres) <- names(coms)
	for (c in 1:length(coms)) 
	{	
		# create list to hold each set of results within each community 
		NRI_abund[[c]] <- NRI_pres[[c]] <- vector("list",length=length(phyNames))
		names(NRI_abund[[c]]) <- names(NRI_pres[[c]]) <- phyNames
		for (p in 1:length(phylos)) 
		{ 	
			# read in results
			NRI_abund[[c]][[p]] <- read.csv(paste0(NRIDir, "NRI_abund_", names(coms)[c], "_", phyNames[p], ".csv"), row.names=1)	
			NRI_pres[[c]][[p]] <- read.csv(paste0(NRIDir, "NRI_pres_", names(coms)[c], "_", phyNames[p], ".csv"), row.names=1)		
		}
	}

#########################################################################
# Summarize and Plot NRI


#################################	
# single representative plot: individual
ctnew <- gsub("fo","hy",crosstype)
NRIplot <- data.frame(ctnew,year=substr(year,3,4),NRI=NRI_abund[[1]][[1]][,6])
boxplot(NRI ~ ctnew + year, NRIplot, at=c(1,2,3,5,6,7,9,10,11,13,14,15), col=c("blue","green3","yellow"))

#################################	
# single representative plot: pooled
barplot(NRI_abund$full_pooled[[1]][,6], space = c(0.2,0.2,0.2,1,0.2,0.2,1,0.2,0.2,1,0.2,0.2), col=c("blue","green3","yellow"))




#########################################################################
# NRI: Individual
	# NRI: Summary statistics of individual NRI results
	source('~/Projects/Phylomeet/Analysis/IndivSummaryMaker.R', chdir = TRUE)
	NRIindiv <- vector("list",length=4)
	names(NRIindiv) <- c("NRI_full_indiv_abund", "NRI_full_indiv_pres", "NRI_common_indiv_abund", "NRI_common_indiv_pres")
	NRIindiv[[1]] <- IndivSummaryMaker(result=NRI_abund, community=1, column=8)
	NRIindiv[[2]] <- IndivSummaryMaker(result=NRI_pres, community=1, column=8)
	NRIindiv[[3]] <- IndivSummaryMaker(result=NRI_abund, community=2, column=8)
	NRIindiv[[4]] <- IndivSummaryMaker(result=NRI_pres, community=2, column=8)


for(c in 1:length(NRI_abund))
{
	for(p in 1:length(NRI_abund[[1]]))
	{
		df <- data.frame(year, ct=crosstype, nri=NRI_abund[[c]][[p]]$mpd.obs.z)
		m <- acast(df, ct ~ year, mean)
		s <- acast(df, year ~ ct, sd)
	}
}



		barplot(-m, col = c("black", "red", "blue"), beside=TRUE)



	# NRI: Plot individual means
	for(i in 1:4)
	{
		pdf( file=paste0( figDir, names(NRIindiv[i]), ".pdf"), width=10, height=10 )
		par(mfrow=c(4,5))
		for (p in 1:20) 
		{ 	barplot(-NRIindiv[[i]][[p]][,3], col = c("black", "red", "blue"), space=spaces, main=names(NRIindiv[[i]])[[p]], ylim=c(-2,2)) 	}
		dev.off()
	}	
	


	# NRI: Pooled by crosstype within years
	for(i in 3:4)
	{
		pdf( file=paste0( figDir, "NRI_", names(NRI_abund[i]), "_abund.pdf" ), width=10, height=10 )
		par(mfrow=c(4,5))
		for (p in 1:20) 
		{ 	barplot(-NRI_abund[[i]][[p]][,6], space = c(0.2,0.2,0.2,1,0.2,0.2,1,0.2,0.2,1,0.2,0.2), col = c("black", "red", "blue"), main=names(NRI_abund[[i]][p]), ylim=c(-2,2))	}
		dev.off()
		
		pdf( file=paste0( figDir, "NRI_", names(NRI_pres[i]), "_pres.pdf" ), width=10, height=10 )
		par(mfrow=c(4,5))
		for (p in 1:20) 
		{ 	barplot(-NRI_pres[[i]][[p]][,6], space = c(0.2,0.2,0.2,1,0.2,0.2,1,0.2,0.2,1,0.2,0.2), col = c("black", "red", "blue"), main=names(NRI_pres[[i]][p]), ylim=c(-2,2))	}
		dev.off()
	}

	# NRI: Pooled by cross type
	for(i in 5:6)
	{
		pdf( file=paste0( figDir, "NRI_", names(NRI_abund[i]), "_abund.pdf" ), width=10, height=10 )
		par(mfrow=c(4,5))
		for (p in 1:20) 
		{ 	barplot(-NRI_abund[[i]][[p]][,6], col = c("black", "red", "blue"), main=names(NRI_abund[[i]][p]), ylim=c(-2,2))	}
		dev.off()
		
		pdf( file=paste0( figDir, "NRI_", names(NRI_abund[i]), "_pres.pdf" ), width=10, height=10 )
		par(mfrow=c(4,5))
		for (p in 1:20) 
		{ 	barplot(-NRI_pres[[i]][[p]][,6], col = c("black", "red", "blue"), main=names(NRI_pres[[i]][p]), ylim=c(-2,2))	}
		dev.off()
	}

	

#################################
# Phylogenetic Beta Diversity: Community Distance (ComDist)
#################################
# make a list to contain comdist results
ComDist_abund <- ComDist_pres <- vector("list",length=length(coms))
names(ComDist_abund) <- names(ComDist_pres) <- names(coms)

# calculate comdist for all communities and branch lengths
	# loop through community subsets
	for (c in 1:length(coms)) 
	{	
		# list to contain results for phylogeny
		ComDist_abund[[c]] <- ComDist_pres[[c]] <- vector("list",length=length(phyNames)) 
		names(ComDist_abund[[c]]) <- names(ComDist_pres[[c]]) <- phyNames
		
		# loop through phylogenies within each community subset
		for (p in 1:length(phylos)) 
		{ 	
			ComDist_abund[[c]][[p]] <- comdist( coms[[c]], phydist[[c]][[p]], abundance.weighted=TRUE )
			ComDist_pres[[c]][[p]] <- comdist( coms[[c]], phydist[[c]][[p]], abundance.weighted=FALSE )
			write.csv( as.matrix( ComDist_abund[[c]][[p]]), paste0( comDistDir, "ComDist_abund_", names(coms)[c], "_", phyNames[p], ".csv"))
			write.csv( as.matrix( ComDist_pres[[c]][[p]]), paste0( comDistDir, "ComDist_pres_", names(coms)[c], "_", phyNames[p], ".csv"))
		}
	}

# import results
	for (c in 1:length(coms)) 
	{	
		# list to contain results for phylogeny
		ComDist_abund[[c]] <- ComDist_pres[[c]] <- vector("list",length=length(phyNames)) 
		names(ComDist_abund[[c]]) <- names(ComDist_pres[[c]]) <- phyNames
		
	# import results
	for (p in 1:length(phylos)) 
		{ 	
			ComDist_abund[[c]][[p]] <- dist( read.csv( paste0( comDistDir, "ComDist_abund_", names(coms)[c], "_", phyNames[p], ".csv" ), row.names=1 ) ) 
			ComDist_pres[[c]][[p]] <- dist( read.csv( paste0( comDistDir, "ComDist_pres_", names(coms)[c], "_", phyNames[p], ".csv" ), row.names=1 ) ) 
		}
		names(ComDist_abund[[c]]) <- names(ComDist_pres[[c]]) <- phyNames
	}
		ComDist <- list(ComDist_abund, ComDist_pres)
		names(ComDist) <- c("ComDist_abund", "ComDist_pres")

############################
# graph of comdist results vs host type
cd <- ComDist_abund$full_pooled$ultra
ty <- t(combn(CtYrNames,2))
c12 <- colsplit(string=ty[,1], pattern="20", names=c("t1","y1"))
c34 <- colsplit(string=ty[,2], pattern="20", names=c("t2","y2"))
m <- data.frame(c12, c34, dist=as.numeric(cd))

fofo <- m[m$t1 == "fo" & m$t2 == "fo",]
fofr <- rbind( m[m$t1 == "fo" & m$t2 == "fr",] , m[m$t1 == "fr" & m$t2 == "fo",] )
frfr <- m[m$t1 == "fr" & m$t2 == "fr",]




fofo <- m[grepl("fo",m[,1]),]

fo <- m[unique(c(grep("fo", m[,1]), grep("fo", m[,2]))),]
plot(fo[,3])

m[grep("fo", m[,2]),]
plot(m[,3],)


df <- data.frame(a1=c("x1","y","x2","y"),a2=c("x","x","y","y"), a3=c(1,2,3,4))

df[df$a1 == "x" & df$a2 == "x",]
df[df$a1 == df[grep("x",df$a1),] & df$a2 == "x",]
df[grep("x",df$a1) & df$a2 == "x",]
df[unique(c(grep("x",df$a1), grep("x",df$a2))),]

df[grep("x",df$a1),]
df[grep("y",df$a1),]


#############################
# Nonmetric MDS 2D
	library(MASS)
	library(car)
	type <- rep(c("Fre","F1","Nar"),4)		# names for labels
	cols <- rep(c("black","red","blue"),4)	# colors for labels
	syms <- rep(c(15,16,17),4)

	# plot into pdf files
	for( a in 1:length(ComDist))
	{
		for( c in 1:length(ComDist_abund) )
		{
			pdf( file=paste0( figDir, names(ComDist[a]), "_", names(ComDist[[a]][c]), ".pdf" ), width=10, height=10 )
			par( mfrow=c(4,5) )
			for( p in 1:20 )
			{
				fit <- isoMDS(ComDist[[a]][[c]][[p]], k=2)
				x <- fit$points[,1]
				y <- fit$points[,2]
				plot( x, y, col=cols, pch=syms, cex=2, main = paste0( names( ComDist[[a]][c]), " ", names( ComDist[[a]][[c]][p] ) ) )
			}
			dev.off()
		}		
	}


	
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

