# getting set up
#install.packages("picante", dependencies = TRUE)
library(picante)

# Picante manual
vignette("picante-intro")

# load nexus files
phy_equal <- read.nexus("/Users/kjj/Documents/_NAU/Classes/G2E/Phylomeet/Analysis/arthrocot/arthrocot_equal/arthrocot_equal.trees")
phy_grad <- read.nexus("/Users/kjj/Documents/_NAU/Classes/G2E/Phylomeet/Analysis/arthrocot/arthrocot_grad/arthrocot_grad.trees")
phy_gradrand <- read.nexus("/Users/kjj/Documents/_NAU/Classes/G2E/Phylomeet/Analysis/arthrocot/arthrocot_gradrand/arthrocot_gradrand.trees")
phy_ultra <- read.nexus("/Users/kjj/Documents/_NAU/Classes/G2E/Phylomeet/Analysis/arthrocot/arthrocot_ultra/arthrocot_ultra.trees")
phy_ultrarand <- read.nexus("/Users/kjj/Documents/_NAU/Classes/G2E/Phylomeet/Analysis/arthrocot/arthrocot_ultrarand/arthrocot_ultrarand.trees")

# plot
plot(phy_equal, cex=0.25)
title(main = "Equal Branch Lengths")

plot(phy_grad, cex=0.25)
title(main = "Graduated Branch Lengths")

plot(phy_gradrand, cex=0.25) 
title(main = "Randomized Graduated Branch Lengths")

plot(phy_ultra, cex=0.25) 
title(main = "Ultrametric Branch Lengths")

plot(phy_ultrarand, cex=0.25) 
title(main = "Randomized Ultrametric Branch Lengths")

# plot topologies without tips
par(mfrow=c(2,3))

plot(phy_equal, show.tip.label = FALSE)
title(main = "Equal Branch Lengths")

plot(phy_grad, show.tip.label = FALSE)
title(main = "Graduated Branch Lengths")

plot(phy_gradrand, show.tip.label = FALSE) 
title(main = "Randomized Graduated Branch Lengths")

plot(phy_ultra, show.tip.label = FALSE) 
title(main = "Ultrametric Branch Lengths")

plot(phy_ultrarand, show.tip.label = FALSE) 
title(main = "Randomized Ultrametric Branch Lengths")