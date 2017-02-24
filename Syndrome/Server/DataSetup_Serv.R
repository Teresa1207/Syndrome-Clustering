### Syndrome Classification    ###
### Model Setup and Simulation ###
#setwd("/Users/Teresa/Documents/Dissertation/R_Dissertation/Syndrome")

###

############# Description ###############
# Our goal is to define two syndromes that lead to dementia (AD).
# A syndrome is defined by the progression pattern of key AD Biomarkers
# More specifically, when are these markers hitting certain percentiles
# Which marker is advancing first?
# For this analysis we will have 2 syndromes. 
# S1: Syndrome 1: Amy-->HPCV-->Vascular
# S2: Syndrome 2: Vascular-->HPCV---->AMY
# In syndrome 1, Amyloid advances to the 50th percentile prior to HPCV followed by vascular. 
# In syndrome 2, Vascular advances to the 50th percentile prior to HPCV, followed by AMY
# In syndrome 3, HPCV advances to the 50th percentile prior to AMY, followed by Vascular

#########################################
library(mvtnorm)
library(MASS)
library(reshape2)
## Data ##


sdSim = c(.01,.05,.1,.5)
pMat = as.matrix(rbind(rep(1/3,3),c(.5,.25,.25),c(.5,.40,.1)),nrow = 3)

set1 = readRDS('set1a.RDA') #Syndrome 1
set2 = readRDS('set2a.RDA') #Syndrome 2
set3 = readRDS('set3a.RDA') #Syndrome 3

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
df$Syndrome = sample(1:3,prob = pMat[pScen,],n,replace=TRUE) #proportion of syndrome 1,2
df.long = merge(data.frame(RID  = unlist(lapply(1:n, function(x){rep(x,nobs)})),oTime = rep(0:(nobs-1),n)),
                          df,by = 'RID', all.x = TRUE)
df.long$dAge = with(df.long, dTime+oTime)
df.long$tAge = with(df.long, oAge+oTime)

fun1 = function(t,B0,B1){pnorm(B0+B1*t+rnorm(n=1,mean = 0,sd = sdSim[sdScen]))}

df.long$M1 = unlist(lapply(1:n, function(x){
  dfsub = subset(df.long,RID==x)
  s = unique(dfsub$Syndrome)
  if(s==1){
    setsub=set1[x,]
    del = setsub$delta}else{
      if(s==2){
      setsub=set2[x,]
      del = setsub$delta-setsub$lag1-setsub$lag2}else{
      setsub=set3[x,]
      del = setsub$delta-setsub$lag1}}
M1 = fun1(t = (dfsub$dAge+del),B0 = setsub$B0,B1 = setsub$B1)
  }))
df.long$M2 = unlist(lapply(1:n, function(x){
  dfsub = subset(df.long,RID==x)
  s = unique(dfsub$Syndrome)
  if(s==1){
    setsub=set1[x,]
    del = setsub$delta-setsub$lag1}else{
      if(s==2){
      setsub=set2[x,]
      del = setsub$delta-setsub$lag1}else{
        setsub=set3[x,]
        del = setsub$delta}}
  M2 = fun1(t = (dfsub$dAge+del),B0 = setsub$C0,B1 = setsub$C1)
}))

df.long$M3 = unlist(lapply(1:n, function(x){
  dfsub = subset(df.long,RID==x)
  s = unique(dfsub$Syndrome)
  if(s==1){
    setsub=set1[x,]
    del = setsub$delta-setsub$lag1-setsub$lag2}else{
      if(s==2){
      setsub=set2[x,]
    del = setsub$delta}else{
      setsub=set3[x,]
      del = setsub$delta-setsub$lag1-setsub$lag2}}
  M3 = fun1(t = (dfsub$dAge+del),B0 = setsub$D0,B1 = setsub$D1)
}))

  
  
df.melt = melt(data = df.long, id.vars = c(1:10),measure.vars = c(11:13))
df.melt$error = rnorm(n = nrow(df.melt),mean = 0.00,0.00)#, sd = sdSim[sdScen])
df.melt$ve = with(df.melt, value+error)
df.melt$ve[df.melt$ve==0]=0.0001
df.melt$ve[df.melt$ve==1]=.9999

saveRDS(df.long,'setL.RDA')
saveRDS(df.melt,'setM.RDA')



