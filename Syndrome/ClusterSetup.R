### Syndrome Classification    ###
### Model Setup and Simulation ###
v# S2: Syndrome 2: Vascular-->HPCV, Amy is idependent. 
# In syndrome 1, Amyloid advances to the 50th percentile prior to HPCV followed by vascular. 
# In syndrome 2, Vascular advances to the 50th percentile prior to HPCV. Amyloid is indepenent. 

#########################################
library(mvtnorm)
library(MASS)
library(ggplot2)
library(GGally)
library(reshape2)
library(plyr)
`%ni%` <- Negate(`%in%`)
## Data ##
setwd("/Users/Teresa/Documents/Dissertation/R_Dissertation/Syndrome")

source('ModelSetup.R')
source('DataSetup.R')

set1 = readRDS('set1a.RDA') #Syndrome 1
set2 = readRDS('set2a.RDA') #Syndrome 2
set3 = readRDS('set3a.RDA') #df.long
set4 = readRDS('set4a.RDA') #df.melt

n = nrow(set1)

## Plot by Marker
s4 = subset(set4,dTime<30&dTime>0)

d1 = subset(s4,variable=='M1')
d2 = subset(s4,variable=='M2')
d3 = subset(s4,variable=='M3')
s1 = subset(set4,Syndrome==1)
s2 = subset(set4,Syndrome==2)

ggplot(d1, aes(x = tAge,y = ve,group=RID,col = factor(Syndrome)))+geom_line()
ggplot(d2, aes(x = tAge,y = ve,group=RID,col = factor(Syndrome)))+geom_line()
ggplot(d3, aes(x = tAge,y = ve,group=RID,col = factor(Syndrome)))+geom_line()

ggplot(s1, aes(x = dAge,y = ve,group=c(RID),col = variable))+ggtitle('Syndrome 1')+
  geom_line(aes(subset(s1,variable=='M1')$dAge,subset(s1,variable=='M1')$ve,group=subset(s1,variable=='M1')$RID),color='darkred')+
  geom_line(aes(subset(s1,variable=='M2')$dAge,subset(s1,variable=='M2')$ve,group=subset(s1,variable=='M2')$RID),color='darkgreen')+
  geom_line(aes(subset(s1,variable=='M3')$dAge,subset(s1,variable=='M3')$ve,group=subset(s1,variable=='M3')$RID),color='darkblue')+
  stat_smooth(aes(group = variable,color = variable),size=2)


ggplot(s2, aes(x = dAge,y = ve,group=c(RID),col = variable))+ggtitle('Syndrome 2')+
  geom_line(aes(subset(s2,variable=='M1')$dAge,subset(s2,variable=='M1')$ve,group=subset(s2,variable=='M1')$RID),color='darkred')+
  geom_line(aes(subset(s2,variable=='M2')$dAge,subset(s2,variable=='M2')$ve,group=subset(s2,variable=='M2')$RID),color='darkgreen')+
  geom_line(aes(subset(s2,variable=='M3')$dAge,subset(s2,variable=='M3')$ve,group=subset(s2,variable=='M3')$RID),color='darkblue')+
  stat_smooth(aes(group = variable),size=2)
ggplot(s2, aes(x = tAge,y = ve,group=c(RID),col = variable))+ggtitle('Syndrome 2')+
  geom_line(aes(subset(s2,variable=='M1')$tAge,subset(s2,variable=='M1')$ve,group=subset(s2,variable=='M1')$RID),color='darkred')+
  geom_line(aes(subset(s2,variable=='M2')$tAge,subset(s2,variable=='M2')$ve,group=subset(s2,variable=='M2')$RID),color='darkgreen')+
  geom_line(aes(subset(s2,variable=='M3')$tAge,subset(s2,variable=='M3')$ve,group=subset(s2,variable=='M3')$RID),color='darkblue')+
  stat_smooth(aes(group = variable),size=2)
s = set4

ggplot(s, aes(x = dAge,y = ve,group=c(RID),col = variable))+ggtitle('Total')+
  geom_line(aes(subset(s,variable=='M1')$dAge,subset(s,variable=='M1')$ve,group=subset(s,variable=='M1')$RID),color='darkred')+
  geom_line(aes(subset(s,variable=='M2')$dAge,subset(s,variable=='M2')$ve,group=subset(s,variable=='M2')$RID),color='darkgreen')+
  geom_line(aes(subset(s,variable=='M3')$dAge,subset(s,variable=='M3')$ve,group=subset(s,variable=='M3')$RID),color='darkblue')+
  stat_smooth(aes(group = variable),size = 2)

ggplot(s, aes(x = tAge,y = ve,group=c(RID),col = variable))+ggtitle('Total')+
  geom_line(aes(subset(s,variable=='M1')$tAge,subset(s,variable=='M1')$ve,group=subset(s,variable=='M1')$RID),color='darkred')+
  geom_line(aes(subset(s,variable=='M2')$tAge,subset(s,variable=='M2')$ve,group=subset(s,variable=='M2')$RID),color='darkgreen')+
  geom_line(aes(subset(s,variable=='M3')$tAge,subset(s,variable=='M3')$ve,group=subset(s,variable=='M3')$RID),color='darkblue')+
  stat_smooth(aes(group = variable),size = 2)

ggplot(set3, aes(x = M3,y = M2,group=RID,color=factor(Syndrome)))+ggtitle('Total')+geom_line()

###############################################################################
####### Cluster Design #######
###############################################################################
setwd("/Users/Teresa/Documents/Dissertation/R_Dissertation/Syndrome")
counter = 1
nRun = 10
errT1 = matrix(0,nrow = nRun,ncol = 2)
errT2 = matrix(0,nrow = nRun,ncol = 2)
errU = matrix(0,nrow = nRun,ncol = 2)
while(counter<(nRun+1)){
source('ModelSetup2.R')
source('DataSetup2.R')

set1 = readRDS('set1a.RDA') #Syndrome 1
set2 = readRDS('set2a.RDA') #Syndrome 2
set3 = readRDS('set3a.RDA') #df.long
set4 = readRDS('set4a.RDA') #df.melt
model <- function(x){
  fit1 = try(lm(qnorm(ve) ~ oTime+I(tAge-65), data=x))
  data.frame(qL=coef(fit1)[[1]], qS=coef(fit1)[[2]])
}

#s4 = subset(set4,value<.80&value>.25)
s4 = subset(set4,dTime<30&dTime>0)

s4 = set4
clust = ddply(s4,.(RID,variable),function(x){model(x)})

clustL = dcast(clust[,c("RID",'variable',"qL","qS")],RID~variable,value.var = c('qL'))
colnames(clustL) = c('RID','L1',"L2","L3")
clustS = dcast(clust[,c("RID",'variable',"qL","qS")],RID~variable,value.var = c('qS'))
colnames(clustS) = c('RID','S1',"S2","S3")
cCast = merge(clustL,clustS,by = 'RID',all = TRUE)

#mark.cast = dcast(markers[,c("RID",'type',"","qS")],RID~type,value.var = c('qL'))

#ggpairs(M1cast[,-1])

##Transform Betas##
trans = function(x){scale(x,scale=FALSE)}
cComp = cCast[complete.cases(cCast),]
cData = data.frame(cComp$RID,t(data.frame(apply(cComp[,c('L1','L2','L3')],1,scale))),
                    t(data.frame(apply(cComp[,c('S1','S2','S3')],1,scale))))
colnames(cData) = c("RID",'L1t','L2t','L3t','S1t','S2t','S3t')

cData1 = ddply(cComp,.(RID),summarise,
               b12 = (L1-L2)*(nobs-1)+(S1 - S2)*(nobs-1)^2/2,
               b13 = (L1-L3)*(nobs-1)+(S1 - S3)*(nobs-1)^2/2,
               b23 = (L2-L3)*(nobs-1)+(S2 - S3)*(nobs-1)^2/2
)


bFrame = data.frame(b1=-cComp[,2]/cComp[,5],b2=-cComp[,3]/cComp[,6],b3=-cComp[,4]/cComp[,7])
bFrame1 = t(apply(bFrame,1,rank))
d1 = cData
d2 = subset(cCast,RID%in%cData$RID)

k1 = data.frame(RID = d1$RID,clustB =kmeans(cData1[,-1],2)$cluster, clustT = kmeans(d1[,-1],2)$cluster,clustU = kmeans(d2[,-1],2)$cluster)
#k2 = Ω
k1.s = merge(k1,set3[,c('RID',"Syndrome",'dGrp')],by = 'RID',all.x= TRUE)
k1.d = merge(k1,d1,by = 'RID')
#ddply(k1.d,.(clustT),summarise, mL1 = mean(L1t),mL2 = mean(L2t), mL3 = mean(L3t),mS1 = mean(S1t),mS2 = mean(S2t),mS3 = mean(S3t))


cT = xtabs(~clustT+Syndrome,k1.s)/6
cB = xtabs(~clustB+Syndrome,k1.s)/6
cU = xtabs(~clustU+Syndrome,k1.s)/6

errT1[counter,] = round(c(min(cT[,1])/sum(cT[,1]),min(cT[,2])/sum(cT[,2]))*100,2)
errT2[counter,] = round(c(min(cB[,1])/sum(cB[,1]),min(cB[,2])/sum(cB[,2]))*100,2)
errU[counter,] = round(c(min(cU[,1])/sum(cU[,1]),min(cU[,2])/sum(cU[,2]))*100,2)
counter = counter+1
}
err1 = data.frame(errT1,errT2,errU)
k2 = merge(cData1,k1,by = 'RID')



xtabs(~clustU+Syndrome,k1.s)/6
xtabs(~clustU+dGrp,k1.s)/6


m1 = merge(subset(s4,RID%in%k1$RID),k1,by = "RID",all.x = TRUE)
m2 = merge(subset(set3,RID%in%k1$RID),k1,by = 'RID',all.x = TRUE)
m=m2
ggplot(m,aes(x = M1,y = M3,group = RID,color = factor(clustB)))+geom_line()+facet_grid(~clustT)
ggplot(m,aes(x = M1,y = M3,group = RID,color = factor(clustB)))+geom_line()+
  stat_smooth(aes(group = Syndrome),color = 'black')

par(mfrow = c(2,1))
s = subset(m1,clustB == 1)
ggplot(s, aes(x = dAge,y = ve,group=c(RID),col = variable))+ggtitle('Total')+
  geom_line(aes(subset(s,variable=='M1')$dAge,subset(s,variable=='M1')$ve,group=subset(s,variable=='M1')$RID),color='darkred')+
  geom_line(aes(subset(s,variable=='M2')$dAge,subset(s,variable=='M2')$ve,group=subset(s,variable=='M2')$RID),color='darkgreen')+
  geom_line(aes(subset(s,variable=='M3')$dAge,subset(s,variable=='M3')$ve,group=subset(s,variable=='M3')$RID),color='darkblue')+
  stat_smooth(aes(group = variable),size = 2)

s = subset(m1,clustB == 2)
ggplot(s, aes(x = dAge,y = ve,group=c(RID),col = variable))+ggtitle('Total')+
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
