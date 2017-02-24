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
`%ni%` <- Negate(`%in%`)
###############################################################################
####### Cluster Design #######
###############################################################################
counter = 1
nRun = 10
errT1 = matrix(0,nrow = nRun,ncol = 2)
errT2 = matrix(0,nrow = nRun,ncol = 2)
errU = matrix(0,nrow = nRun,ncol = 2)
silsum = matrix(0,nrow = nRun,ncol = 3)
while(counter<(nRun+1)){
source('ModelSetup2.R')
source('DataSetup2.R')

set1 = readRDS('set1a.RDA') #Syndrome 1
set2 = readRDS('set2a.RDA') #Syndrome 2
set3 = readRDS('set3a.RDA') #df.long
set4 = readRDS('set4a.RDA') #df.melt
model <- function(x){
  fit1 = try(lm(qnorm(ve) ~ oTime, data=x))
  data.frame(qL=coef(fit1)[[1]], qS=coef(fit1)[[2]])
}


s4 = set4
clust = ddply(s4,.(RID,variable),function(x){model(x)})

clustL = dcast(clust[,c("RID",'variable',"qL","qS")],RID~variable,value.var = c('qL'))
colnames(clustL) = c('RID','L1',"L2","L3")
clustS = dcast(clust[,c("RID",'variable',"qL","qS")],RID~variable,value.var = c('qS'))
colnames(clustS) = c('RID','S1',"S2","S3")
cCast = merge(clustL,clustS,by = 'RID',all = TRUE)

##Transform Betas##
trans = function(x){scale(x,scale=FALSE)}
cComp = cCast[complete.cases(cCast),]
cData = data.frame(cComp$RID,t(data.frame(apply(cComp[,c('L1','L2','L3')],1,scale))),
                    t(data.frame(apply(cComp[,c('S1','S2','S3')],1,scale))))
colnames(cData) = c("RID",'L1t','L2t','L3t','S1t','S2t','S3t')
obs = nobs
cData1 = ddply(cComp,.(RID),summarise,
               b12 = (L1-L2)*(obs-1)+(S1 - S2)*(obs-1)^2/2,
               b13 = (L1-L3)*(obs-1)+(S1 - S3)*(obs-1)^2/2,
               b23 = (L2-L3)*(obs-1)+(S2 - S3)*(obs-1)^2/2
)

dist1 = dist(cData1[,-1])
k1 = data.frame(RID = cData1$RID,clustB =kmeans(cData1[,-1],2,nstart=20)$cluster, clustT = kmeans(cData[,-1],2,nstart = 20)$cluster,clustU = kmeans(cCast[,-1],2,nstart = 20)$cluster)
k1.s = merge(k1,set3[,c('RID',"Syndrome",'dGrp')],by = 'RID',all.x= TRUE)
#k1.d = merge(k1,d1,by = 'RID')
k3 = data.frame(RID = cData$RID,clustB3 =kmeans(cData1[,-1],3)$cluster, clustT = kmeans(cData[,-1],3)$cluster,clustU = kmeans(cCast[,-1],3)$cluster)
k4 = data.frame(RID = cData$RID,clustB4 =kmeans(cData1[,-1],4)$cluster, clustT = kmeans(cData[,-1],4)$cluster,clustU = kmeans(cCast[,-1],4)$cluster)
k1.full = merge(k1,k3[,c("RID",'clustB3')],by = "RID",all.x  = TRUE)
k1.full  = merge(k1.full,k4[,c("RID","clustB4")],by = 'RID',all.x = TRUE)


cT = xtabs(~clustT+Syndrome,k1.s)/6
cB = xtabs(~clustB+Syndrome,k1.s)/6
cU = xtabs(~clustU+Syndrome,k1.s)/6
silsum[counter,] =round(c(summary(silhouette(k1.full$clustB,dist1))$avg.width,
                    summary(silhouette(k1.full$clustB3,dist1))$avg.width,
                    summary(silhouette(k1.full$clustB4,dist1))$avg.width),4)
errT1[counter,] = round(c(min(cT[,1])/sum(cT[,1]),min(cT[,2])/sum(cT[,2]))*100,2)
errT2[counter,] = round(c(min(cB[,1])/sum(cB[,1]),min(cB[,2])/sum(cB[,2]))*100,2)
errU[counter,] = round(c(min(cU[,1])/sum(cU[,1]),min(cU[,2])/sum(cU[,2]))*100,2)
counter = counter+1
}
err1 = data.frame(errT1,errT2,errU)

k2 = merge(cData1,k1,by = 'RID')
k2a = merge(k2,subset(k1.s[!duplicated(k1.s$RID),],select = c(RID,Syndrome)),by = 'RID',all.x=TRUE)
sum1 = ddply(k2a,.(clustB),summarise,
            d12 = mean(b12),
            d13 = mean(b13),
            d23 = mean(b23))

sum2 = ddply(k2a,.(Syndrome),summarise,
             d12 = mean(b12),
             d13 = mean(b13),
             d23 = mean(b23))







m1 = merge(subset(s4,RID%in%k1$RID),k1,by = "RID",all.x = TRUE)
m2 = merge(subset(set3,RID%in%k1$RID),k1,by = 'RID',all.x = TRUE)
m=m2
ggplot(m,aes(x = M1,y = M3,group = RID,color = factor(clustB)))+geom_line()+facet_grid(~clustT)
ggplot(m,aes(x = M1,y = M3,group = RID,color = factor(clustB)))+geom_line()+
  stat_smooth(aes(group = Syndrome),color = 'black')

s = m1
ggplot(s, aes(x = dAge,y = ve,group=c(RID),col = variable))+ggtitle('Total')+xlim(0,30)+
  geom_line(aes(subset(s,variable=='M1')$dAge,subset(s,variable=='M1')$ve,group=subset(s,variable=='M1')$RID),color='darkred')+
  geom_line(aes(subset(s,variable=='M2')$dAge,subset(s,variable=='M2')$ve,group=subset(s,variable=='M2')$RID),color='darkgreen')+
  geom_line(aes(subset(s,variable=='M3')$dAge,subset(s,variable=='M3')$ve,group=subset(s,variable=='M3')$RID),color='darkblue')+
  stat_smooth(aes(group = variable),size = 2)

s = subset(m1,clustB == 2)
ggplot(s, aes(x = dAge,y = ve,group=c(RID),col = variable))+ggtitle('Clust Group 2')+xlim(0,30)+
  geom_line(aes(subset(s,variable=='M1')$dAge,subset(s,variable=='M1')$ve,group=subset(s,variable=='M1')$RID),color='darkred')+
  geom_line(aes(subset(s,variable=='M2')$dAge,subset(s,variable=='M2')$ve,group=subset(s,variable=='M2')$RID),color='darkgreen')+
  geom_line(aes(subset(s,variable=='M3')$dAge,subset(s,variable=='M3')$ve,group=subset(s,variable=='M3')$RID),color='darkblue')+
  stat_smooth(aes(group = variable),size = 2)

s = subset(m1,clustB == 1)
ggplot(s, aes(x = dAge,y = ve,group=c(RID),col = variable))+ggtitle('Clust Group 1')+xlim(0,30)+
  geom_line(aes(subset(s,variable=='M1')$dAge,subset(s,variable=='M1')$ve,group=subset(s,variable=='M1')$RID),color='darkred')+
  geom_line(aes(subset(s,variable=='M2')$dAge,subset(s,variable=='M2')$ve,group=subset(s,variable=='M2')$RID),color='darkgreen')+
  geom_line(aes(subset(s,variable=='M3')$dAge,subset(s,variable=='M3')$ve,group=subset(s,variable=='M3')$RID),color='darkblue')+
  stat_smooth(aes(group = variable),size = 2)


s = m1
ggplot(s, aes(x = tAge,y = ve,group=c(RID),col = variable))+ggtitle('Total')+
  geom_line(aes(subset(s,variable=='M1')$tAge,subset(s,variable=='M1')$ve,group=subset(s,variable=='M1')$RID),color='darkred')+
  geom_line(aes(subset(s,variable=='M2')$tAge,subset(s,variable=='M2')$ve,group=subset(s,variable=='M2')$RID),color='darkgreen')+
  geom_line(aes(subset(s,variable=='M3')$tAge,subset(s,variable=='M3')$ve,group=subset(s,variable=='M3')$RID),color='darkblue')+
  stat_smooth(aes(group = variable),size = 2)

ggplot(m,aes(x = M2,y = M3,group = RID,color = factor(Syndrome)))+geom_line()+
  stat_smooth(aes(group = factor(Syndrome)), color = 'black')

ggplot(m,aes(x = M2,y = M3,group = RID,color = factor(clustB)))+geom_line()+
  stat_smooth(aes(group = factor(clustB)), color = 'black')
