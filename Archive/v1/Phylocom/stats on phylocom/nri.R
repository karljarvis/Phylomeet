# I want to do an ANOVA of NRI values(a measure of phylogenetic relatedness)

setwd("/Users/karljarvis/Documents/_NAU/Classes/G2E/Community Phylogenetics/R")

#missing data is NaN or NA
arthro_nri <- read.csv("arthro_nri.csv")
attach(arthro_nri)

#check for normal looking boxplots
par(mfrow = c(1,1))
boxplot(NRI ~ host_type)

arthronona <- na.omit(arthro_nri)

par(mfrow = c(3,1))

R <- subset(NRI, host_type == "R")
hist(R, 20, xlim =c(-2,4))
abline(v=mean(R))

S <- subset(arthronona$NRI, host_type == "S")
hist(S, 20, xlim =c(-2,4))
abline(v=mean(S))

D <- subset(NRI, host_type == "D")
hist(D, 20, xlim =c(-2,4))
abline(v=mean(D))




