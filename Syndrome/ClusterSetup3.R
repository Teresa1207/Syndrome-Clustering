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

`%ni%` <- Negate(`%in%`)
###############################################################################
####### Cluster Design #######
###############################################################################
#set.seed(12071983)
counter = 1
nRun = 500
errT1 = matrix(0,nrow = nRun,ncol = 2)
errT2 = matrix(0,nrow = nRun,ncol = 2)
errU = matrix(0,nrow = nRun,ncol = 2)
silsum = matrix(0,nrow = nRun,ncol = 3)
while(counter<(nRun+1)){
source('ModelSetup3.R')
source('DataSetup2.R')

set1 = readRDS('set1a.RDA') #Syndrome 1
set2 = readRDS('set2a.RDA') #Syndrome 2
set3 = readRDS('set3a.RDA') #df.long
set4 = readRDS('set4a.RDA') #df.melt
model <- function(x){
  fit1 = try(lm(qnorm(ve) ~ oTime,data=x))
  #fit2 = try(lm(qnorm(value) ~ oTime+oAge, data=x))
  data.frame(qL=coef(fit1)[[1]], qS=coef(fit1)[[2]],qL1=coef(fit1)[[1]], qS1=coef(fit1)[[2]])
}


s4 = set4
clust = ddply(s4,.(RID,variable),function(x){model(x)})

clustL = dcast(clust[,c("RID",'variable',"qL","qS")],RID~variable,value.var = c('qL'))
colnames(clustL) = c('RID','L1',"L2","L3")
clustS = dcast(clust[,c("RID",'variable',"qL","qS")],RID~variable,value.var = c('qS'))
colnames(clustS) = c('RID','S1',"S2","S3")
cCast = merge(clustL,clustS,by = 'RID',all = TRUE)

clustL1 = dcast(clust[,c("RID",'variable',"qL1","qS1")],RID~variable,value.var = c('qL1'))
colnames(clustL1) = c('RID','L1',"L2","L3")
clustS1 = dcast(clust[,c("RID",'variable',"qL1","qS1")],RID~variable,value.var = c('qS1'))
colnames(clustS1) = c('RID','S1',"S2","S3")
cCast1 = merge(clustL1,clustS1,by = 'RID',all = TRUE)



##Transform Betas##
trans = function(x){scale(x,scale=FALSE)}
cComp = cCast[complete.cases(cCast1),]
cData = data.frame(cComp$RID,t(data.frame(apply(cComp[,c('L1','L2','L3')],1,scale))),
                    t(data.frame(apply(cComp[,c('S1','S2','S3')],1,scale))))
colnames(cData) = c("RID",'L1t','L2t','L3t','S1t','S2t','S3t')
obs = 6
cData1 = ddply(cComp,.(RID),summarise,
               b12 = (L1-L2)*(obs-1)+(S1 - S2)*(obs-1)^2/2,
               b13 = (L1-L3)*(obs-1)+(S1 - S3)*(obs-1)^2/2,
               b23 = (L2-L3)*(obs-1)+(S2 - S3)*(obs-1)^2/2
)

k1 = data.frame(RID = cData1$RID,clustB =kmeans(cData1[,-1],2,nstart=50)$cluster)#, clustT = kmeans(cData[,-1],2,nstart = 50)$cluster,clustU = kmeans(cCast1[,-1],2,nstart = 50)$cluster)
k1.s = merge(k1,set3[,c('RID',"Syndrome",'dGrp')],by = 'RID',all.x= TRUE)

cB = xtabs(~clustB+Syndrome,k1.s)/6

#errT1[counter,] = round(c(min(cT[,1])/sum(cT[,1]),min(cT[,2])/sum(cT[,2]))*100,2)
errT2[counter,] = round(c(min(cB[,1])/sum(cB[,1]),min(cB[,2])/sum(cB[,2]))*100,4)
#errU[counter,] = round(c(min(cU[,1])/sum(cU[,1]),min(cU[,2])/sum(cU[,2]))*100,2)
counter = counter+1
}
err1 = data.frame(errT2)

e5 = err1
e5$MR = apply(e5,1,sum)
saveRDS(e5,file = 'e5.RDA')


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
m2 = merge(subset(set3,RID%in%k1$RID),k1,by = 'RID',all.x = TRUE)
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

