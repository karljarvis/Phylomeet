setwd("/Users/karljarvis/Documents/_NAU/Classes/G2E/Community Phylogenetics/R")
data <- read.csv("myco_sr.csv")
is.data.frame(data)
attach(data)
S
R
par(mfrow=c(3,1))
hist(S, 20)
hist(R, 20)


boxplot(data)
wilcox.test(S,R, alternative = "greater")

wilcox.test(S,0,alternative = 'greater')