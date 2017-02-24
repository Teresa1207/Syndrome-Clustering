### Syndrome Classification    ###
### Model Setup and Simulation ###
# S2: Syndrome 2: Vascular-->HPCV, Amy is idependent. 
# In syndrome 1, Amyloid advances to the 50th percentile prior to HPCV followed by vascular. 
# In syndrome 2, Vascular advances to the 50th percentile prior to HPCV. Amyloid is indepenent. 
setwd("/Users/Teresa/Documents/Dissertation/R_Dissertation/Syndrome")
#########################################
library(mvtnorm)
library(MASS)
library(ggplot2)
library(GGally)
library(reshape2)
library(plyr)
library(geepack)
library(flexclust)
library(kml3d)

`%ni%` <- Negate(`%in%`)
###############################################################################
####### Cluster Design #######
###############################################################################
pScen = 1
sdScen = 1
varScen = 1

#set.seed(12071983)
counter = 1
nRun = 5

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
               b13 = (L1-L3)*(obs-1)+(S1 - S3)*(obs-1)^2/2,
               b23 = (L2-L3)*(obs-1)+(S2 - S3)*(obs-1)^2/2
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
vals <- vector()
for (k in 2:6) {
  # Perform the kmeans algorithm
  cl <- kmeans(cData1[,-1], k,nstart=50)
  # Compute the Calinski_Harabasz index
  vals <- c(vals,as.numeric(intCriteria(as.matrix(cData1[,c(2:4)]),cl$cluster,"Calinski_Harabasz")))
}

idx[counter,1] <- bestCriterion(vals,"Calinski_Harabasz")+1
idx[counter,2] <- bestCriterion(vals,"Ray_Turi")+1
idx[counter,3] <- bestCriterion(vals,"SD_Scat")+1
idx[counter,4] <- bestCriterion(vals,"Silhouette")+1
counter = counter+1
}

scen1 = data.frame(clus = idx,TRand = Index[,2],KMRand = IndexKM[,2],S1Rand = IndexS1[,2],S2Rand =IndexS2[,2],S3Rand =IndexS3[,2])

saveRDS(scen1,"resultsSim/scen1.RDA")



k2 = merge(cData1,k1,by = 'RID')
k2a = merge(k2,subset(k1.s[!duplicated(k1.s$RID),],select = c(RID,Syndrome)),by = 'RID',all.x=TRUE)
b2 = merge(clust,k1,by = 'RID')
b2a = merge(b2,k2a,by = 'RID',all.x=TRUE)

sum1 = ddply(k2a,.(clustB),summarise,
            d12 = mean(b12),
            d13 = mean(b13),
            d23 = mean(b23))

sum2 = ddply(k2a,.(Syndrome),summarise,
             d12 = mean(b12),
             d13 = mean(b13),
             d23 = mean(b23))







m1Sim = merge(subset(s4,RID%in%k1$RID),k1,by = "RID",all.x = TRUE)
m2 = merge(subset(setLong,RID%in%k1$RID),k1,by = 'RID',all.x = TRUE)
m=m2
ggplot(m,aes(x = M1,y = M3,group = RID,color = factor(clustB)))+geom_line()+facet_grid(~clustT)
ggplot(m,aes(x = M1,y = M3,group = RID,color = factor(clustB)))+geom_line()+
  stat_smooth(aes(group = Syndrome),color = 'black')
myColors <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

s = m1Sim
ggplot(s, aes(x = dAge,y = ve,group=c(RID),col = variable))+ggtitle('Clust Group 2')+xlim(0,30)+
  geom_line(data = subset(s,variable=='M1'),aes(x=dAge,y = ve,group=RID,color = "M1"))+
  geom_line(data = subset(s,variable=='M2'),aes(x=dAge,y = ve,group=RID,color = "M2"))+
  geom_line(data = subset(s,variable=='M3'),aes(x=dAge,y = ve,group=RID,color = "M3"))+
  stat_smooth(aes(group = variable),size = 2)+scale_color_manual(values=myColors)

s = subset(m1Sim,clustB == 2)
ggplot(s, aes(x = dAge,y = ve,group=c(RID),col = variable))+ggtitle('Clust Group 2')+xlim(0,30)+
  geom_line(data = subset(s,variable=='M1'),aes(x=dAge,y = ve,group=RID,color = "M1"))+
  geom_line(data = subset(s,variable=='M2'),aes(x=dAge,y = ve,group=RID,color = "M2"))+
  geom_line(data = subset(s,variable=='M3'),aes(x=dAge,y = ve,group=RID,color = "M3"))+
  stat_smooth(aes(group = variable),size = 2)+scale_color_manual(values=myColors)

s = subset(m1Sim,clustB == 1)
ggplot(s, aes(x = dAge,y = ve,group=c(RID),col = variable))+ggtitle('Clust Group 2')+xlim(0,30)+
  geom_line(data = subset(s,variable=='M1'),aes(x=dAge,y = ve,group=RID,color = "M1"))+
  geom_line(data = subset(s,variable=='M2'),aes(x=dAge,y = ve,group=RID,color = "M2"))+
  geom_line(data = subset(s,variable=='M3'),aes(x=dAge,y = ve,group=RID,color = "M3"))+
  stat_smooth(aes(group = variable),size = 2)+scale_color_manual(values=myColors)

s = m1Sim
ggplot(s, aes(x = tAge,y = ve,group=c(RID),col = variable))+ggtitle('Whole Group')+
  geom_line(data = subset(s,variable=='M1'),aes(x=tAge,y = ve,group=RID,color = "M1"))+
  geom_line(data = subset(s,variable=='M2'),aes(x=tAge,y = ve,group=RID,color = "M2"))+
  geom_line(data = subset(s,variable=='M3'),aes(x=tAge,y = ve,group=RID,color = "M3"))+
  #stat_smooth(aes(group = variable),size = 2)
  scale_color_manual(values=myColors)

s = subset(m1Sim,clustB == 2)
ggplot(s, aes(x = tAge,y = ve,group=c(RID),col = variable))+ggtitle('Clust Group 2')+
  geom_line(data = subset(s,variable=='M1'),aes(x=tAge,y = ve,group=RID,color = "M1"))+
  geom_line(data = subset(s,variable=='M2'),aes(x=tAge,y = ve,group=RID,color = "M2"))+
  geom_line(data = subset(s,variable=='M3'),aes(x=tAge,y = ve,group=RID,color = "M3"))+
  stat_smooth(aes(group = variable),size = 2)+scale_color_manual(values=myColors)

s = m1Sim
ggplot(s, aes(x = oTime,y = ve,group=c(RID),col = factor(clustB)))+ggtitle('Clust Group 2')+
  geom_line(data = subset(s,variable=='M1'),aes(x=oTime,y = ve,group=RID,color = "M1"))+
  geom_line(data = subset(s,variable=='M2'),aes(x=oTime,y = ve,group=RID,color = "M2"))+
  geom_line(data = subset(s,variable=='M3'),aes(x=oTime,y = ve,group=RID,color = "M3"))+
  #stat_smooth(aes(group = variable),size = 2)+
  scale_color_manual(values=myColors)


ggplot(m,aes(x = M2,y = M3,group = RID,color = factor(Syndrome)))+geom_line()+
  stat_smooth(aes(group = factor(Syndrome)), color = 'black')

ggplot(m,aes(x = M2,y = M3,group = RID,color = factor(clustB)))+geom_line()+
  stat_smooth(aes(group = factor(clustB)), color = 'black')

