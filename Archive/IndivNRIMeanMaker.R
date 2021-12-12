
IndivNRIMeanMaker <- function(result=NRI_abund, community=1) 
{
	mn <- se <- vector("list",length=length(phyNames))
	names(mn) <- names(se) <- phyNames

	for (p in 1:length(phylos)) 
	{ 	
		df <- data.frame(year, ct=ctyear, nri=result[[community]][[p]]$mpd.obs.z)
		mn[[p]] <- acast(df, ct ~ year, mean)
		# se[[p]] <- acast(df, ct ~ year, function(x) {sd(x)/sqrt(length(x))})
	}
	mn
	# se
}

test <- IndivNRIMeanMaker(NRI_abund,1)