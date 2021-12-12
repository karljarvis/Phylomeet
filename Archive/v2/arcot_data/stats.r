setwd("/Users/kjj/Documents/Phylomeet/Analysis/arcot_common2_2y")

##################
# Faith's PD
##################
pd.u <- read.csv("pd.u.csv")

fr00.pd.u <- as.numeric(pd.u[1:10,2])
fo00.pd.u <- as.numeric(pd.u[11:20,2])
na00.pd.u <- as.numeric(pd.u[21:40,2])
fr01.pd.u <- as.numeric(pd.u[41:50,2])
fo01.pd.u <- as.numeric(pd.u[51:60,2])
na01.pd.u <- as.numeric(pd.u[61:80,2])
fr02.pd.u <- as.numeric(pd.u[81:90,2])
fo02.pd.u <- as.numeric(pd.u[91:100,2])
na02.pd.u <- as.numeric(pd.u[101:120,2])
fr03.pd.u <- as.numeric(pd.u[121:130,2])
fo03.pd.u <- as.numeric(pd.u[131:140,2])
na03.pd.u <- as.numeric(pd.u[141:160,2])

# standard errors
se <- function(x) sd(x)/sqrt(length(x))

fr00.se.u <- se(fr00.pd.u)
fo00.pd.u <- as.numeric(pd.u[11:20,2])
na00.pd.u <- as.numeric(pd.u[21:40,2])
fr01.pd.u <- as.numeric(pd.u[41:50,2])
fo01.pd.u <- as.numeric(pd.u[51:60,2])
na01.pd.u <- as.numeric(pd.u[61:80,2])
fr02.pd.u <- as.numeric(pd.u[81:90,2])
fo02.pd.u <- as.numeric(pd.u[91:100,2])
na02.pd.u <- as.numeric(pd.u[101:120,2])
fr03.pd.u <- as.numeric(pd.u[121:130,2])
fo03.pd.u <- as.numeric(pd.u[131:140,2])
na03.pd.u <- as.numeric(pd.u[141:160,2])

# bar graph
library(gplots)
pd.u.means <- c(mean(fr00.pd.u), mean(fo00.pd.u), mean(na00.pd.u), mean(fr01.pd.u), mean(fo01.pd.u), mean(na01.pd.u), mean(fr02.pd.u), mean(fo02.pd.u), mean(na02.pd.u), mean(fr03.pd.u), mean(fo03.pd.u), mean(na03.pd.u))

barplot(pd.u.means, space = c(0.2,0.2,0.2,1,0.2,0.2,1,0.2,0.2,1,0.2,0.2), ylim = c(0,120), col = c("blue", "darkgreen", "gold"))
plotCI(pd.u.means, uiw = )


# histograms of nri
par(mfrow=c(4,3))
hist(fr00.pd.u)
hist(fo00.pd.u)
hist(na00.pd.u)
hist(fr01.pd.u)
hist(fo01.pd.u)
hist(na01.pd.u)
hist(fr02.pd.u)
hist(fo02.pd.u)
hist(na02.pd.u)
hist(fr03.pd.u)
hist(fo03.pd.u)
hist(na03.pd.u)

################
# PD QUESTIONS #
################

######################################################
# Do the host types differ across the whole dataset?
######################################################
pd.u.k <- kruskal.test(list(fr00.pd.u, fo00.pd.u, na00.pd.u, fr01.pd.u, fo01.pd.u, na01.pd.u, fr02.pd.u, fo02.pd.u, na02.pd.u, fr03.pd.u, fo03.pd.u, na03.pd.u)); pd.u.k
# TS = 74, p = e-11
pd.u.k <- kruskal.test(list(fr01.pd.u, fo01.pd.u, na01.pd.u, fr02.pd.u, fo02.pd.u, na02.pd.u, fr03.pd.u, fo03.pd.u, na03.pd.u)); pd.u.k
# TS=65, p = e-11
# Yup.



######################################################
# Does PD differ by host type?
######################################################
# mann-whitney test
frfo00.pd.u.mw <- wilcox.test(fr00.pd.u, fo00.pd.u); frfo00.pd.u.mw
nafo00.pd.u.mw <- wilcox.test(na00.pd.u, fo00.pd.u); nafo00.pd.u.mw
frfo01.pd.u.mw <- wilcox.test(fr01.pd.u, fo01.pd.u); frfo01.pd.u.mw
nafo01.pd.u.mw <- wilcox.test(na01.pd.u, fo01.pd.u); nafo01.pd.u.mw
frfo02.pd.u.mw <- wilcox.test(fr02.pd.u, fo02.pd.u); frfo02.pd.u.mw
nafo02.pd.u.mw <- wilcox.test(na02.pd.u, fo02.pd.u); nafo02.pd.u.mw
frfo03.pd.u.mw <- wilcox.test(fr03.pd.u, fo03.pd.u); frfo03.pd.u.mw
nafo03.pd.u.mw <- wilcox.test(na03.pd.u, fo03.pd.u); nafo03.pd.u.mw

# output
pd.ph.mw <- c(frfo00.pd.u.mw$p.value, nafo00.pd.u.mw$p.value, frfo01.pd.u.mw$p.value, nafo01.pd.u.mw$p.value, frfo02.pd.u.mw$p.value, nafo02.pd.u.mw$p.value, frfo03.pd.u.mw$p.value, nafo03.pd.u.mw$p.value)
ph.names <- c("frfo00", "nafo00", "frfo01", "nafo01", "frfo02", "nafo02", "frfo03", "nafo03")
pd.ph.mat <- cbind(ph.names, pd.ph.mw)
write.csv(pd.ph.mat, "pd.ph.csv")


#############
# nri ultra #
#############
nri.u <- read.csv("nri_ultra_rich.csv")
nri.u <- as.matrix(nri.u)

fr00.nri.u <- -as.numeric(nri.u[1:10,7])
fo00.nri.u <- -as.numeric(nri.u[11:20,7])
na00.nri.u <- -as.numeric(nri.u[21:40,7])
fr01.nri.u <- -as.numeric(nri.u[41:50,7])
fo01.nri.u <- -as.numeric(nri.u[51:60,7])
na01.nri.u <- -as.numeric(nri.u[61:80,7])
fr02.nri.u <- -as.numeric(nri.u[81:90,7])
fo02.nri.u <- -as.numeric(nri.u[91:100,7])
na02.nri.u <- -as.numeric(nri.u[101:120,7])
fr03.nri.u <- -as.numeric(nri.u[121:130,7])
fo03.nri.u <- -as.numeric(nri.u[131:140,7])
na03.nri.u <- -as.numeric(nri.u[141:160,7])

fr00m.nri.u <- mean(fr00.nri.u)
fo00m.nri.u <- mean(fo00.nri.u)
na00m.nri.u <- mean(na00.nri.u)
fr01m.nri.u <- mean(fr01.nri.u)
fo01m.nri.u <- mean(fo01.nri.u)
na01m.nri.u <- mean(na01.nri.u)
fr02m.nri.u <- mean(fr02.nri.u)
fo02m.nri.u <- mean(fo02.nri.u)
na02m.nri.u <- mean(na02.nri.u)
fr03m.nri.u <- mean(fr03.nri.u)
fo03m.nri.u <- mean(fo03.nri.u)
na03m.nri.u <- mean(na03.nri.u)


# bar graph of 2001
nri.01.u <- c(mean(fr01.nri.u), mean(fo01.nri.u), mean(na01.nri.u))
barplot(nri.01.u, col=c("blue", "darkgreen", "gold"), axes = T, ylim=c(0,1.3))

# bar graph of mean of all
nri.aves.u <- c(mean(fr00m.nri.u,fr01m.nri.u,fr02m.nri.u,fr03m.nri.u), mean(fo00m.nri.u, fo01m.nri.u, fo02m.nri.u, fo03m.nri.u), mean(na00m.nri.u, na01m.nri.u, na02m.nri.u, na03m.nri.u))
barplot(nri.aves.u, col=c("blue", "darkgreen", "gold"), ylim = c(0,1.6))

# bar graph of all
nri.00_03.u <- c(mean(fr00.nri.u), mean(fo00.nri.u), mean(na00.nri.u), mean(fr01.nri.u), mean(fo01.nri.u), mean(na01.nri.u), mean(fr02.nri.u), mean(fo02.nri.u), mean(na02.nri.u), mean(fr03.nri.u), mean(fo03.nri.u), mean(na03.nri.u))
barplot(nri.00_03.u, space = c(0.2,0.2,0.2,1,0.2,0.2,1,0.2,0.2,1,0.2,0.2), col=c("blue", "darkgreen", "gold"), ylab = "NRI")

# histograms of nri.u
par(mfrow=c(4,3))
hist(fr00.nri.u)
hist(fo00.nri.u)
hist(na00.nri.u)
hist(fr01.nri.u)
hist(fo01.nri.u)
hist(na01.nri.u)
hist(fr02.nri.u)
hist(fo02.nri.u)
hist(na02.nri.u)
hist(fr03.nri.u)
hist(fo03.nri.u)
hist(na03.nri.u)

#################
# NRI QUESTIONS #
#################

######################################################
# Is there non-neutral phylogenetic structure?
######################################################
# wilcoxon test on all host type groups
fr00.nri.u.w <- wilcox.test(fr00.nri.u); fr00.nri.u.w
fo00.nri.u.w <- wilcox.test(fo00.nri.u); fo00.nri.u.w
na00.nri.u.w <- wilcox.test(na00.nri.u); na00.nri.u.w
fr01.nri.u.w <- wilcox.test(fr01.nri.u); fr01.nri.u.w
fo01.nri.u.w <- wilcox.test(fo01.nri.u); fo01.nri.u.w
na01.nri.u.w <- wilcox.test(na01.nri.u); na01.nri.u.w
fr02.nri.u.w <- wilcox.test(fr02.nri.u); fr02.nri.u.w
fo02.nri.u.w <- wilcox.test(fo02.nri.u); fo02.nri.u.w
na02.nri.u.w <- wilcox.test(na02.nri.u); na02.nri.u.w
fr03.nri.u.w <- wilcox.test(fr03.nri.u); fr03.nri.u.w
fo03.nri.u.w <- wilcox.test(fo03.nri.u); fo03.nri.u.w
na03.nri.u.w <- wilcox.test(na03.nri.u); na03.nri.u.w

# extract p-values
nri.u.w <- as.numeric(c(fr00.nri.u.w$p.value, fo00.nri.u.w$p.value, na00.nri.u.w$p.value, fr01.nri.u.w$p.value, fo01.nri.u.w$p.value, na01.nri.u.w$p.value, fr02.nri.u.w$p.value, fo02.nri.u.w$p.value, na02.nri.u.w$p.value, fr03.nri.u.w$p.value, fo03.nri.u.w$p.value, na03.nri.u.w$p.value))

# make character vector of group names, bind to p-values, export to csv
comm.labels <- c("fr00", "fo00", "na00", "fr01", "fo01", "na01", "fr02", "fo02", "na02", "fr03", "fo03", "na03")
nri.u.w.out <- cbind(comm.labels, nri.u.w)
write.csv(nri.u.w.out, "nri.u.w.out.csv")

######################################################
# Do the host types differ across the whole dataset?
######################################################
# kruskal-wallis test
nri.u.k <- kruskal.test(list(fr00.nri.u, fo00.nri.u, na00.nri.u, fr01.nri.u, fo01.nri.u, na01.nri.u, fr02.nri.u, fo02.nri.u, na02.nri.u, fr03.nri.u, fo03.nri.u, na03.nri.u))
nri.u.k

# resulting p-value <<< 0.0001
# Chi2 test statistic: 63.3
# So: yes they most definitely differ

# kruskal-wallis test without 2000
nri.u.k <- kruskal.test(list(fr01.nri.u, fo01.nri.u, na01.nri.u, fr02.nri.u, fo02.nri.u, na02.nri.u, fr03.nri.u, fo03.nri.u, na03.nri.u)); nri.u.k
# same result, TS of 53 instead, p-value still <<<0.0001


######################################################
# Do the host types differ within each year?
######################################################
# kruskal-wallis test
nri00.u.k <- kruskal.test(list(fr00.nri.u, fo00.nri.u, na00.nri.u)); nri00.u.k
nri01.u.k <- kruskal.test(list(fr01.nri.u, fo01.nri.u, na01.nri.u)); nri01.u.k
nri02.u.k <- kruskal.test(list(fr02.nri.u, fo02.nri.u, na02.nri.u)); nri02.u.k
nri03.u.k <- kruskal.test(list(fr03.nri.u, fo03.nri.u, na03.nri.u)); nri03.u.k

# output
nri.htwy <- c(nri00.u.k$p.value, nri01.u.k$p.value, nri02.u.k$p.value, nri03.u.k$p.value)
year.labels <- c("2000", "2001", "2002", "2003")
nri.htwy.m <- cbind(year.labels, nri.htwy)
write.csv(nri.htwy.m, "nri.htwy.csv")

######################################################
# Do the parental species differ from hybrids (within each year)? 
######################################################
# mann-whitney test
frfo00.nri.u.mw <- wilcox.test(fr00.nri.u, fo00.nri.u); frfo00.nri.u.mw
nafo00.nri.u.mw <- wilcox.test(na00.nri.u, fo00.nri.u); nafo00.nri.u.mw
frfo01.nri.u.mw <- wilcox.test(fr01.nri.u, fo01.nri.u); frfo01.nri.u.mw
nafo01.nri.u.mw <- wilcox.test(na01.nri.u, fo01.nri.u); nafo01.nri.u.mw
frfo02.nri.u.mw <- wilcox.test(fr02.nri.u, fo02.nri.u); frfo02.nri.u.mw
nafo02.nri.u.mw <- wilcox.test(na02.nri.u, fo02.nri.u); nafo02.nri.u.mw
frfo03.nri.u.mw <- wilcox.test(fr03.nri.u, fo03.nri.u); frfo03.nri.u.mw
nafo03.nri.u.mw <- wilcox.test(na03.nri.u, fo03.nri.u); nafo03.nri.u.mw

# output
nri.ph.mw <- c(frfo00.nri.u.mw$p.value, nafo00.nri.u.mw$p.value, frfo01.nri.u.mw$p.value, nafo01.nri.u.mw$p.value, frfo02.nri.u.mw$p.value, nafo02.nri.u.mw$p.value, frfo03.nri.u.mw$p.value, nafo03.nri.u.mw$p.value)
ph.names <- c("frfo00", "nafo00", "frfo01", "nafo01", "frfo02", "nafo02", "frfo03", "nafo03")
nri.ph.mat <- cbind(ph.names, nri.ph.mw)
write.csv(nri.ph.mat, "nri.ph.csv")

######################################################
# Does each host type differ from year to year? 
######################################################
# kruskal-wallis test
nri.fr.u.k <- kruskal.test(list(fr00.nri.u, fr01.nri.u, fr02.nri.u, fr03.nri.u))
nri.fo.u.k <- kruskal.test(list(fo00.nri.u, fo01.nri.u, fo02.nri.u, fo03.nri.u))
nri.na.u.k <- kruskal.test(list(na00.nri.u, na01.nri.u, na02.nri.u, na03.nri.u))

# output
nri.htd <- c(nri.fr.u.k$p.value, nri.fo.u.k$p.value, nri.na.u.k$p.value)
htd.names <- c("fr", "fo", "na")
htd.mat <- cbind(nri.htd,htd.names)
write.csv(htd.mat, "htdy.csv")

# year to year comparison
nri.fr0001.u.mw <- wilcox.test(fr00.nri.u, fr01.nri.u)
nri.fr0102.u.mw <- wilcox.test(fr01.nri.u, fr02.nri.u)
nri.fr0203.u.mw <- wilcox.test(fr02.nri.u, fr03.nri.u)
nri.fo0001.u.mw <- wilcox.test(fo00.nri.u, fo01.nri.u)
nri.fo0102.u.mw <- wilcox.test(fo01.nri.u, fo02.nri.u)
nri.fo0203.u.mw <- wilcox.test(fo02.nri.u, fo03.nri.u)
nri.na0001.u.mw <- wilcox.test(na00.nri.u, na01.nri.u)
nri.na0102.u.mw <- wilcox.test(na01.nri.u, na02.nri.u)
nri.na0203.u.mw <- wilcox.test(na02.nri.u, na03.nri.u)

# output
nri.htd.fr <- c(nri.fr0001.u.mw$p.value, nri.fr0102.u.mw$p.value, nri.fr0203.u.mw$p.value)
nri.htd.fo <- c(nri.fo0001.u.mw$p.value, nri.fo0102.u.mw$p.value, nri.fo0203.u.mw$p.value)
nri.htd.na <- c(nri.na0001.u.mw$p.value, nri.na0102.u.mw$p.value, nri.na0203.u.mw$p.value)
htd.namesffn <- c("2000-01", "2001-02", "2002-03")
htd.ffn.mat <- cbind(htd.namesffn, nri.htd.fr, nri.htd.fo, nri.htd.na)
write.csv(htd.ffn.mat, "htdy.ffn.csv")

##################
# 2001 nti data
##################
nti <- read.csv("nti_ultra_rich.csv")
nti <- as.matrix(nti)

fr01.nti <- as.numeric(nti[41:50,7])
fo01.nti <- as.numeric(nti[51:60,7])
na01.nti <- as.numeric(nti[61:80,7])

# bar graph
nti01.nti <- c(mean(fr01.nti), mean(fo01.nti), mean(na01.nti))
barplot(nti01.nti)

# histograms of nti
par(mfrow=c(1,3))
hist(fr01.nti)
hist(fo01.nti)
hist(na01.nti)

# wilcoxon test
fr01.nti.w <- wilcox.test(fr01.nti); fr01.nti.w
fo01.nti.w <- wilcox.test(fo01.nti); fo01.nti.w
na01.nti.w <- wilcox.test(na01.nti); na01.nti.w

##################
# 2002 nri data
##################
nri <- read.csv("nri_ultra_rich.csv")
nri <- as.matrix(nri)

fr02.nri <- as.numeric(nri[41:50,7])
fo02.nri <- as.numeric(nri[51:60,7])
na02.nri <- as.numeric(nri[61:80,7])

# bar graph
nri02.nri <- c(mean(fr02.nri), mean(fo02.nri), mean(na02.nri))
barplot(nri02.nri)

# histograms of nri
par(mfrow=c(1,3))
hist(fr02.nri)
hist(fo02.nri)
hist(na02.nri)

# wilcoxon test
fr02.nri.w <- wilcox.test(fr02.nri); fr02.nri.w
fo02.nri.w <- wilcox.test(fo02.nri); fo02.nri.w
na02.nri.w <- wilcox.test(na02.nri); na02.nri.w

##################
# 2002 nti data
##################
nti <- read.csv("nti_ultra_rich.csv")
nti <- as.matrix(nti)

fr02.nti <- as.numeric(nti[41:50,7])
fo02.nti <- as.numeric(nti[51:60,7])
na02.nti <- as.numeric(nti[61:80,7])

# bar graph
nti02.nti <- c(mean(fr02.nti), mean(fo02.nti), mean(na02.nti))
barplot(nti02.nti)

# histograms of nti
par(mfrow=c(1,3))
hist(fr02.nti)
hist(fo02.nti)
hist(na02.nti)

# wilcoxon test
fr02.nti.w <- wilcox.test(fr02.nti); fr02.nti.w
fo02.nti.w <- wilcox.test(fo02.nti); fo02.nti.w
na02.nti.w <- wilcox.test(na02.nti); na02.nti.w

##################
# 2003 nri data
##################
nri <- read.csv("nri_ultra_rich.csv")
nri <- as.matrix(nri)

fr03.nri <- as.numeric(nri[41:50,7])
fo03.nri <- as.numeric(nri[51:60,7])
na03.nri <- as.numeric(nri[61:80,7])

# bar graph
nri03.nri <- c(mean(fr03.nri), mean(fo03.nri), mean(na03.nri))
barplot(nri03.nri)

# histograms of nri
par(mfrow=c(1,3))
hist(fr03.nri)
hist(fo03.nri)
hist(na03.nri)

# wilcoxon test
fr03.nri.w <- wilcox.test(fr03.nri); fr03.nri.w
fo03.nri.w <- wilcox.test(fo03.nri); fo03.nri.w
na03.nri.w <- wilcox.test(na03.nri); na03.nri.w

##################
# 2003 nti data
##################
nti <- read.csv("nti_ultra_rich.csv")
nti <- as.matrix(nti)

fr03.nti <- as.numeric(nti[41:50,7])
fo03.nti <- as.numeric(nti[51:60,7])
na03.nti <- as.numeric(nti[61:80,7])

# bar graph
nti03.nti <- c(mean(fr03.nti), mean(fo03.nti), mean(na03.nti))
barplot(nti03.nti)

# histograms of nti
par(mfrow=c(1,3))
hist(fr03.nti)
hist(fo03.nti)
hist(na03.nti)

# wilcoxon test
fr03.nti.w <- wilcox.test(fr03.nti); fr03.nti.w
fo03.nti.w <- wilcox.test(fo03.nti); fo03.nti.w
na03.nti.w <- wilcox.test(na03.nti); na03.nti.w