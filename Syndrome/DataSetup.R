### Syndrome Classification    ###
### Model Setup and Simulation ###
setwd("/Users/Teresa/Documents/Dissertation/R_Dissertation/Syndrome")

###

############# Description ###############
# Our goal is to define two syndromes that lead to dementia (AD).
# A syndrome is defined by the progression pattern of key AD Biomarkers
# More specifically, when are these markers hitting certain percentiles
# Which marker is advancing first?
# For this analysis we will have 2 syndromes. 
# S1: Syndrome 1: Amy-->HPCV-->Vascular
# S2: Syndrome 2: Vascular-->HPCV, Amy is idependent. 
# In syndrome 1, Amyloid advances to the 50th percentile prior to HPCV followed by vascular. 
# In syndrome 2, Vascular advances to the 50th percentile prior to HPCV. Amyloid is indepenent. 

#########################################
library(mvtnorm)
library(MASS)
library(ggplot2)
library(GGally)
library(reshape2)
## Data ##

set1 = readRDS('set1.RDA') #Syndrome 1
set2 = readRDS('set2.RDA') #Syndrome 2
n = nrow(set1)
nobs = 6 #number (yrs) of observations

df = data.frame(RID = 1:n,dGrp = sample(1:4,prob=c(.50,.50,00,0),size = n,replace=TRUE))
#Groups = 1:Early (50), 2:Medium (65), 3: Late (80) 4:None(100)
df$dType = sapply(1:n,function(x){if(df$dGrp[x]==1){
  return(50)}else{
    if(df$dGrp[x]==2){
      return(65)}else{
        if(df$dGrp[x]==3){
          return(80)}else{
            return(100)}
          }
        }
      })
df$dStart = df$dType+rnorm(n, mean = 0, sd = 3)
df$oAge = runif(n,55,85)
df$dTime = df$oAge - df$dStart
df$Syndrome = sample(1:2,n,replace=TRUE) #proportion of syndrome1,2
df.long = merge(data.frame(RID  = unlist(lapply(1:n, function(x){rep(x,nobs)})),oTime = rep(0:(nobs-1),n)),
                          df,by = 'RID', all.x = TRUE)
df.long$dAge = with(df.long, dTime+oTime)
df.long$tAge = with(df.long, oAge+oTime)

fun1 = function(t,B0,B1){pnorm(B0+B1*t+rnorm(1,0,.1))}
Ifun1 = function(t,B0,B1){(qnorm(t) - B0)/B1}

df.long$M1 = unlist(lapply(1:n, function(x){
  dfsub = subset(df.long,RID==x)
  s = unique(dfsub$Syndrome)
  if(s==1){setsub=set1[x,]}else{setsub=set2[x,]}
M1 = fun1(t = dfsub$dAge,B0 = setsub$B0,B1 = setsub$B1)
  }))
df.long$M2 = unlist(lapply(1:n, function(x){
  dfsub = subset(df.long,RID==x)
  s = unique(dfsub$Syndrome)
  if(s==1){setsub=set1[x,]}else{setsub=set2[x,]}
  M2 = fun1(t = dfsub$dAge,B0 = setsub$C0,B1 = setsub$C1)
}))

df.long$M3 = unlist(lapply(1:n, function(x){
  dfsub = subset(df.long,RID==x)
  s = unique(dfsub$Syndrome)
  if(s==1){setsub=set1[x,]}else{setsub=set2[x,]}
  M3 = fun1(t = dfsub$dAge,B0 = setsub$D0,B1 = setsub$D1)
}))

  
  
df.melt = melt(data = df.long, id.vars = c(1:10),measure.vars = c(11:13))
df.melt$error = rnorm(n = nrow(df.melt),mean = 0, sd = 0.01)
df.melt$ve = with(df.melt, value+error)
df.melt$ve[df.melt$ve<0]=0.0001
df.melt$ve[df.melt$ve>1]=.9999
### Plot ###
samp1 = sample(x = 1:n,n,replace=FALSE)
df.melt  = df.melt[order(df.melt$RID,df.melt$variable),]
dataP = subset(df.melt,RID%in%samp1)
ggplot(dataP, aes(x = dAge,y = ve,group=RID,col=factor(Syndrome)))+geom_line()+facet_grid(~variable)+xlim(-10,30)
ggplot(dataP, aes(x = oTime,y = ve,group=RID,col=factor(Syndrome)))+geom_line()+facet_grid(Syndrome~variable)+xlim(0,5)
ggplot(dataP, aes(x = tAge,y = ve,group=RID,col=factor(Syndrome)))+geom_line()+facet_grid(Syndrome~variable)+xlim(60,95)
ggplot(dataP, aes(x = tAge,y = ve,group=RID,col=factor(Syndrome)))+geom_line()+facet_grid(~variable)+xlim(60,95)
ggplot(dataP, aes(x = tAge,y = ve,group=RID,col=factor(Syndrome)))+geom_line()+facet_grid(~variable)+xlim(60,95)

with(df.long,plot(oAge,dStart))

set3 = saveRDS(df.long,'set3.RDA')
set4 = saveRDS(df.melt,'set4.RDA')


