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
library(flexclust)
library(kml3d)
library(clusterCrit)
`%ni%` <- Negate(`%in%`)
###############################################################################
####### Cluster Design #######
###############################################################################
pScen = 2
sdScen = 3
varScen = 8

#set.seed(12071983)
counter = 1
nRun = 2

Index = matrix(0,nrow = nRun,ncol = 4)
IndexS1 = matrix(0,nrow = nRun,ncol = 4)
IndexS2 = matrix(0,nrow = nRun,ncol = 4)
IndexS3 = matrix(0,nrow = nRun,ncol = 4)
IndexKM = matrix(0,nrow = nRun,ncol = 4)
idx = matrix(0,nrow = nRun, ncol = 4)
cS = 3
while(counter<(nRun+1)){
source('ModelSetup_Serv.R')
source('DataSetup_Serv.R')

set1 = readRDS('set1a.RDA') #Syndrome 1
set2 = readRDS('set2a.RDA') #Syndrome 2
set3 = readRDS('set3a.RDA') #Syndrome 3
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

k1 = data.frame(RID = cData1$RID,clustB =kmeans(cData1[,-1],cS,nstart=50)$cluster)

k1.s = merge(k1,setLong[,c('RID',"Syndrome",'dGrp')],by = 'RID',all.x= TRUE)

Sng = data.frame(RID = cComp$RID,clustSng1 =kmeans(cComp[,c(2,5)],cS,nstart=50)$cluster,
                  clustSng2 =kmeans(cComp[,c(3,6)],cS,nstart=50)$cluster,
                  clustSng3 =kmeans(cComp[,c(4,7)],cS,nstart=50)$cluster)
Sng.m = merge(Sng,setLong[,c('RID',"Syndrome",'dGrp')],by = 'RID',all.x= TRUE)
Sng.m1 = Sng.m[!duplicated(Sng.m$RID),]
k1.s1 = k1.s[!duplicated(k1.s$RID),]

Index[counter,] = with(k1.s1,comPart(clustB,Syndrome))
IndexS1[counter,] = with(Sng.m1,comPart(clustSng1,Syndrome))
IndexS2[counter,] = with(Sng.m1,comPart(clustSng2,Syndrome))
IndexS3[counter,] = with(Sng.m1,comPart(clustSng3,Syndrome))

setMelt$tAge2 = round(setMelt$tAge)
rd = dcast(data = setMelt, formula = RID+oTime~variable, value.var = 've')
rd1 = merge(rd,setMelt[!duplicated(setMelt$RID),c("RID",'tAge2')],by = c('RID'), all.x = TRUE)
rd1$tAge3 = rd1$oTime+rd1$tAge2
r = reshape(rd1,idvar = 'RID',timevar = 'tAge3',direction = 'wide')
repL = length(unique(rd1$tAge3))
cldR <- cld3d(r, timeInData = list(M1 = 1:repL * 5-2, M2 = 1:repL *5-1, M3 = 1:repL*5))
kml3d(cldR)
clustKM = getClusters(cldR,cS)
comp = data.frame(df.long[!duplicated(df.long$RID),c('RID','Syndrome')],KML = getClusters(cldR,3))
IndexKM[counter,] = with(comp,comPart(KML,Syndrome))

##clusterCriteria
valsCH <- vector()
valsRT <- vector()
valsSD <- vector()
valsSil <- vector()
for (k in 2:6) {
  # Perform the kmeans algorithm
  cl <- kmeans(cData1[,-1], k,nstart=50)
  # Compute the Calinski_Harabasz index
  valsCH <- c(valsCH,as.numeric(intCriteria(as.matrix(cData1[,-1]),cl$cluster,c("Calinski_Harabasz"))))
  valsRT <- c(valsRT,as.numeric(intCriteria(as.matrix(cData1[,-1]),cl$cluster,c("Ray_Turi"))))
  valsSD <- c(valsSD,as.numeric(intCriteria(as.matrix(cData1[,-1]),cl$cluster,c("SD_Scat"))))
  valsSil <- c(valsSil,as.numeric(intCriteria(as.matrix(cData1[,-1]),cl$cluster,c("Silhouette"))))
  }

idx[counter,1] <- bestCriterion(valsCH,"Calinski_Harabasz")+1
idx[counter,2] <- bestCriterion(valsRT,"Ray_Turi")+1
idx[counter,3] <- bestCriterion(valsSD,"SD_Scat")+1
idx[counter,4] <- bestCriterion(valsSil,"Silhouette")+1
counter = counter+1
}

scen = data.frame(clus = idx,TRand = Index[,2],KMRand = IndexKM[,2],S1Rand = IndexS1[,2],S2Rand =IndexS2[,2],S3Rand =IndexS3[,2])

#saveRDS(scen,paste("scenario_",pScen,sdScen,varScen,".RDA",sep=""))


