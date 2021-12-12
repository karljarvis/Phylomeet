IndivSummaryMaker <- function(result=PD, column=3) 
{
	summary <- vector("list",length=length(phyNames))
	names(summary) <- phyNames
	crosstypes <- c("fr","fo","na")
	years <- 2000:2003
	
	for (i in 1:length(phyNames)) 
	{ 	
		mean <- se <- numeric(0)
		res <- cbind(crosstype, year, result[[1]][[i]])
		for(j in 1:length(crosstypes))
		{
			for(k in 1:length(years))
			{
				mean <- c(mean, mean(res[res[,1] == crosstypes[j] & res[,2] == years[k], column]))
				se <- c(se, sd(res[res[,1] == crosstypes[j] & res[,2] == years[k], column] ) / sqrt( length( res[res[,1] == crosstypes[k] & res[,2] == years[j], column])))
			}
		}
		ct <- rep(c("fr","fo","na"),4)
		yr <- c(rep(2000,3), rep(2001,3), rep(2002,3), rep(2003,3))
		summary[[i]] <- data.frame(crosstype=ct ,year=yr, mean=mean, se=se)
		rownames(summary[[i]]) <- paste0(ct, yr)
	}
	summary
}
