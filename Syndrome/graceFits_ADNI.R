## Not run: 
library(grace)
library(fda)
library(mvtnorm)
setwd("/Users/Teresa/Documents/Dissertation/R_Dissertation/Syndrome")
source('graceT.R')
#clustB
dfGrace = graceC3
gT = with(dfGrace,graceT(argvals = time,y  = quant, outcome = type,id = RID,group = factor(clustB),plots = TRUE))
gT1 = with(subset(dfGrace,clustB==1),graceT(argvals = time,y  = quant, outcome = type,id = RID,group = factor(clustB),plots = FALSE))
gT2 = with(subset(dfGrace,clustB==2),graceT(argvals = time,y  = quant, outcome = type,id = RID,group = factor(clustB),plots = FALSE))
gT3 = with(subset(dfGrace,clustB==3),graceT(argvals = time,y  = quant, outcome = type,id = RID,group = factor(clustB),plots = FALSE))

saveRDS(gT1,'g3T1.RDA')
saveRDS(gT2,'g3T2.RDA')
saveRDS(gT3,'g3T3.RDA')

r = data.frame(rbind(g1$fits$av45,g1$fits$fdg,g1$fits$hoc))
c2g1g2sim =ddply(data.frame(rbind(g1$fits$av45$subset,g1$fits$fdg$subset,g1$fits$hoc$subset)),'outcome',function(x){sum(x$smooth.resids^2)})


