## Not run: 
library(grace)
library(fda)
library(mvtnorm)
setwd("~/Documents/Dissertation/R_Dissertation/Syndrome")
source('graceT.R')
#clustB
dfGrace = graceC3

gTSplit = with(dfGrace,graceT(argvals = time,y  = quant, outcome = type,id = RID,group = factor(clustB),plots = TRUE))
TT1 = with(subset(dfGrace,clustB==1),graceT(argvals = time,y  = quant, outcome = type,id = RID,group = factor(clustB),plots = TRUE))
TT2 = with(subset(dfGrace,clustB==2),graceT(argvals = time,y  = quant, outcome = type,id = RID,group = factor(clustB),plots = TRUE))
TT3 = with(subset(dfGrace,clustB==3),graceT(argvals = time,y  = quant, outcome = type,id = RID,group = factor(clustB),plots = TRUE))
#T4 = with(subset(dfGrace,clustB==4),graceT(argvals = time,y  = quant, outcome = type,id = RID,group = factor(clustB),plots = FALSE))
gTAll = with(dfGrace,graceT(argvals = time,y  = quant, outcome = type,id = RID,plots = TRUE))

saveRDS(gTSplit, 'gtsplit.RDA')
saveRDS(T1,'g4T1.RDA')
saveRDS(T2,'g4T2.RDA')
saveRDS(T3,'g4T3.RDA')
saveRDS(T4,'g4T3.RDA')

g1sim =ddply(data.frame(rbind(g1$fits$av45$subset,g1$fits$fdg$subset,g1$fits$hoc$subset)),'outcome',function(x){sum(x$smooth.resids^2)})
g2sim =ddply(data.frame(rbind(g2$fits$av45$subset,g2$fits$fdg$subset,g2$fits$hoc$subset)),'outcome',function(x){sum(x$smooth.resids^2)})


h1sim =ddply(data.frame(rbind(h1$fits$av45$subset,h1$fits$fdg$subset,h1$fits$hoc$subset)),'outcome',function(x){sum(x$smooth.resids^2)})
h2sim =ddply(data.frame(rbind(h2$fits$av45$subset,h2$fits$fdg$subset,h2$fits$hoc$subset)),'outcome',function(x){sum(x$smooth.resids^2)})
h3sim =ddply(data.frame(rbind(h3$fits$av45$subset,h3$fits$fdg$subset,h3$fits$hoc$subset)),'outcome',function(x){sum(x$smooth.resids^2)})

simFunc  = function(df){ddply(data.frame(rbind(df$fits$av45$subset,df$fits$fdg$subset,df$fits$hoc$subset)),'outcome',function(x){sum(x$smooth.resids^2)})}

T1sim = simFunc(T1)
T2sim = simFunc(T2)
T3sim = simFunc(T3)
T4sim = simFunc(T4)

sum(g1sim$V1,g2sim$V1) + 2*7
sum(h1sim$V1,h2sim$V1,h3sim$V1) + 3*7
sum(T1sim$V1,T2sim$V1,T3sim$V1,T4sim$V1) + 4*7

h2sim
h3sim
