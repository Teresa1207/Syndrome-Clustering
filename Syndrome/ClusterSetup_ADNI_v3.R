### Syndrome Classification    ###
### Model Setup and Simulation ###
# S2: Syndrome 2: Vascular-->HPCV, Amy is idependent. 
# In syndrome 1, Amyloid advances to the 50th percentile prior to HPCV followed by vascular. 
# In syndrome 2, Vascular advances to the 50th percentile prior to HPCV. Amyloid is indepenent. 
#source('~/Documents/Dissertation/R_Dissertation/Syndrome/dataClus_ADNI.R')
#########################################
library(mvtnorm)
library(MASS)
library(ggplot2)
library(GGally)
library(reshape2)
library(plyr)
library(geepack)
library(flexclust)
`%ni%` <- Negate(`%in%`)
###############################################################################
####### Cluster Design #######
###############################################################################
#set.seed(12071983)

source('~/Documents/Dissertation/R_Dissertation/Syndrome/dataClus_ADNI.R')
cS = 3
set4 = subset(mark.long,base3dx!='3:AD')
set4 =ddply(set4,.(type),mutate,
  raw.scl = scale(raw)
) 
model <- function(x){
  fit1 = try(lm(qnorm(quant) ~ time,data=x))
  #fit1 = try(lm(raw.scl ~ time, data=x))
  #data.frame(qL=coef(fit1)[[1]]+coef(fit)[[3]]*unique(x$age65), qS=coef(fit1)[[2]])
  data.frame(qL=coef(fit1)[[1]], qS=coef(fit1)[[2]])
  
}

s4 = set4
s4$quant[which(s4$quant==0)]=0.0001
s4$quant[which(s4$quant==1)]=0.9999

#s4 = subset(s4,base3dx=='2:MCI')

clust = ddply(s4,.(RID,type),function(x){model(x)})

clustL = dcast(clust[,c("RID",'type',"qL","qS")],RID~type,value.var = c('qL'))
colnames(clustL) = c('RID','L1',"L2","L3",'L4')
clustS = dcast(clust[,c("RID",'type',"qL","qS")],RID~type,value.var = c('qS'))
colnames(clustS) = c('RID','S1',"S2","S3",'S4')
cCast = merge(clustL,clustS,by = 'RID',all = TRUE)

##Transform Betas##
cComp = cCast[complete.cases(cCast),]
obs = 3.5
cData1 = ddply(cComp,.(RID),summarise,
               b12 = (L1-L2)*(obs-1)+(S1 - S2)*(obs-1)^2/2,
               b13 = (L1-L3)*(obs-1)+(S1 - S3)*(obs-1)^2/2,
               b14 = (L1-L4)*(obs-1)+(S1 - S4)*(obs-1)^2/2,
               b23 = (L2-L3)*(obs-1)+(S2 - S3)*(obs-1)^2/2,
               b24 = (L2-L4)*(obs-1)+(S2 - S4)*(obs-1)^2/2,
               b34 = (L3-L4)*(obs-1)+(S3 - S4)*(obs-1)^2/2
)

#vars = c('b12','b13','b14','b23','b24','b34')
vars = c('b23','b24','b34')
#vars = c('b12','b13','b14')
cDat = subset(cData1,select = c('RID',vars))
#cDat = clusts
gsP.Z <- clusGap(cDat[,-1], FUN = kmeans,nstart = 50, K.max = 8, B = 500)
fviz_gap_stat(gsP.Z)
res.nb <- NbClust(cDat[,-1], distance = "euclidean",min.nc = 2, max.nc = 8, 
                  method = 'kmeans', index ='alllong') 
# All gap statistic values
res.nb$All.index

# Best number of clusters
res.nb$Best.nc
resDis = res.nb
png('~/Documents/Dissertation/Dissertation_2017/figures/nbclust.png', height = 800, width = 1000)
fviz_nbclust(res.nb) + theme_classic() + theme(title = element_text(size=16),
                                               axis.text.x = element_text(size = 14),
                                               axis.text.y = element_text(size = 14))
dev.off()


k1 = data.frame(RID = cDat$RID,clustB =kmeans(cDat[,-1],cS,nstart=50)$cluster)#, clustT = kmeans(cData[,-1],3,nstart = 50)$cluster,clustU = kmeans(cCast1[,-1],2,nstart = 50)$cluster)
nbPart = data.frame(RID =cDat$RID,nbPart = res.nb$Best.partition, clustB = k1$clustB)

comPart(nbPart$nbPart,nbPart$clustB)

k2 = merge(cDat,k1,by = 'RID')

k2a = merge(k2,dc.u,by= 'RID')

graceC3 = merge(set4,k1)
graceC4 = merge(set4,k1)

##clusterCriteria
library(clusterCrit)
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

idx = c(bestCriterion(valsCH,"Calinski_Harabasz")+1,bestCriterion(valsRT,"Ray_Turi")+1,bestCriterion(valsSD,"SD_Scat")+1,bestCriterion(valsSil,"Silhouette")+1)



