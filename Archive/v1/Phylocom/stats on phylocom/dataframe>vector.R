# the goal of this script is to convert a dataset from 2D to 1D. Basically, I want to take the columns from the dataset and stack them on top of each other with the first column on top, second column beneath the first, and so on.

setwd("/Users/karljarvis/Documents/_NAU/Classes/G2E/Community Phylogenetics")

hosttrees <- read.csv("hosttrees.csv")

#convert the dataframe to a matrix
htm <- as.matrix(hosttrees)

#convert the matrix to a vector
htv <- as.vector(htm)

#export the vector
write.csv(htv, file = "htv.csv")

