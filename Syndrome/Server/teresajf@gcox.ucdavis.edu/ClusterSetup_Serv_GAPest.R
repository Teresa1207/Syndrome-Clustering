### Syndrome Classification    ###
### Model Setup and Simulation ###
# S2: Syndrome 2: Vascular-->HPCV, Amy is idependent. 
# In syndrome 1, Amyloid advances to the 50th percentile prior to HPCV followed by vascular. 
# In syndrome 2, Vascular advances to the 50th percentile prior to HPCV. Amyloid is indepenent. 
#setwd("/Users/Teresa/Documents/Dissertation/R_Dissertation/Syndrome/Server")
#########################################
library(mvtnorm)
library(MASS)
#library(ggplot2)
#library(GGally)
library(reshape2)
library(plyr)
#library(geepack)
library(cluster)
#library(flexclust)
#library(kml3d)
#library(clusterCrit)
`%ni%` <- Negate(`%in%`)
###############################################################################
####### Cluster Design #######
###############################################################################
#setwd("~/Documents/Dissertation/R_Dissertation/Syndrome/Server/ErrorEstServer")
source('param.R')
source('ModelSetup_Serv.R')
source('DataSetup_Serv.R')

res = lapply(c(1:7), function(pdS){
counter = 1
nRun = 500

idx = matrix(0,nrow = nRun, ncol = 1)

#sdRun = sdSim[sdS]
while(counter<(nRun+1)){

modelSet(vs=3)
dataServ(sdS=1,pdS)

setLong = readRDS('setL.RDA') #df.long
setMelt = readRDS('setM.RDA') #df.melt

model <- function(x){
  fit1 = try(lm(qnorm(ve) ~ oTime,data=x))
  data.frame(qL=coef(fit1)[[1]], qS=coef(fit1)[[2]])
}


s4 = setMelt
clust = ddply(s4,.(RID,variable),function(x){model(x)})

clustL = dcast(clust[,c("RID",'variable',"qL","qS")],RID~variable,value.var = c('qL'))
colnames(clustL) = c('RID','L1',"L2","L3")
clustS = dcast(clust[,c("RID",'variable',"qL","qS")],RID~variable,value.var = c('qS'))
colnames(clustS) = c('RID','S1',"S2","S3")
cCast = merge(clustL,clustS,by = 'RID',all = TRUE)

##Transform Betas##
cComp = cCast[complete.cases(cCast),]

obs = 6
cData1 = ddply(cComp,.(RID),summarise,
               b12 = (L1-L2)*(obs-1)+(S1 - S2)*(obs-1)^2/2,
               b13 = (L1-L3)*(obs-1)+(S1 - S3)*(obs-1)^2/2#,
               #b23 = (L2-L3)*(obs-1)+(S2 - S3)*(obs-1)^2/2
)
gsP.Z <- clusGap(cData1[,-1], FUN = kmeans, K.max = 8, B = 200)

idx[counter,] <- maxSE(gsP.Z$Tab[17:24],gsP.Z$Tab[25:32], method = 'Tibs2001SEmax')

counter = counter+1
}

scen = data.frame(clus = idx)
})
saveRDS(res, file = 'resGap1.RDA')


