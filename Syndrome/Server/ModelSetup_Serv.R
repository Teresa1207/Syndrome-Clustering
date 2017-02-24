### Syndrome Classification    ###
### Model Setup and Simulation ###
#setwd("/Users/Teresa/Documents/Dissertation/R_Dissertation/Syndrome/Server")
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
#library(GGally)
## Data ##
myColors <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
t = seq(-10,30,by = 1)
n = 500

sim = matrix(cbind(c(1,2,3,1,1,1,1,3),c(1,1,1,2,3,4,5,3)),ncol = 2)
varPred = c(0.01,0.02,0.03,.04,.05)
sig1 = varPred[sim[varScen,1]]
sig2 = varPred[sim[varScen,2]]                                                    
Corr = matrix(.25,nrow = 6,ncol = 6)
diag(Corr) =1
CovR = diag(rep(c(sig1,sig2),3))%*%Corr%*%diag(rep(c(sig1,sig2),3))
sL = .1 #start level
on = .5 #onset

#S1  AMY: B0 = qnorm(.2), B1 = (qnorm(.5)-qnorm(.2))/10
#S1 HPCV: C0,C1
#S1 Vasc: D0,D1
# Random Start levels, mean B0

iT = 10 #time for first marker
B0 = qnorm(sL)
B1 = (qnorm(on)-qnorm(sL))/iT
A1 = function(t,B0,B1){pnorm(B0+B1*t)}
I1 = function(t,B0,B1){(qnorm(t) - B0)/B1}
I2 = function(t,B0,B1,d,l){(qnorm(t) - B0)/B1-d+l}
A1.t = function(t){A1(t,B0=B0,B1=B1)}

## Generate random levels and slopes for each individual for Marker 1,2,3
b = mvrnorm(n = n,mu = c(0,0,0,0,0,0),Sigma = CovR)
###
syn1 = data.frame(B0 = B0+b[,1], B1 = B1+b[,2])
syn1[,'t1'] = mapply(function(x,y){I1(on,B0=x,B1 = y)},x = syn1$B0,y = syn1$B1)
syn1$delta = syn1$t1-iT

# iT2 = lag of second marker
iT2 = 5
data2 = data.frame(C0 = B0+b[,3], C1 = B1+b[,4])
#data2$lag1 = rgamma(n,shape = 10,rate = 2)
data2$lag1 = 5
syn1[,'t2'] = mapply(function(x,y,z,l){I2(on,B0=x,B1 = y,d = z,l = l )},x = data2$C0,y = data2$C1,z = syn1$delta,l = data2$lag1)
syn1[,'t2a'] = mapply(function(x,y,z,l){I2(on,B0=x,B1 = y,d = z,l = l )},x = data2$C0,y = data2$C1,z = 0,l = data2$lag1)

A2.t = function(t){A1(t-5,B0,B1)}

## Generate random levels and slopes for each individual for Marker 3
iT3 = 5 #lag for marker 3
data3 = data.frame(D0 = B0+b[,5], D1 = B1+b[,6])
#data3$lag2 = rgamma(n,shape = 6,rate = 2)
data3$lag2 = 5
syn1[,'t3'] = mapply(function(x,y,z,l){I2(on,B0=x,B1 = y,d = z,l = l )},x = data3$D0,y = data3$D1,z = syn1$delta,l = (data3$lag2+data2$lag1))
A3.t = function(t){A1(t-10,B0,B1)}

set1 = data.frame(syn1,data2,data3)
set1 = set1[,order(colnames(set1))]

####################################################
####################################################
## Data ##
#S2  AMY: B0 = qnorm(.2), B1 = (qnorm(.5)-qnorm(.2))/10
#S2 HPCV: C0,C1
#S2 Vasc: D0 = qnorm(.2), D1 = (qnorm(.5)-qnorm(.2))/10
# Random Start levels for S3 (Marker 3 initiates)
iT=10
D0 = qnorm(sL)
D1 = (qnorm(on)-qnorm(sL))/iT
B1.t = function(t){A1(t,B0=D0,B1=D1)}
##S2: Generate random levels and slopes for each individual for Marker 3
c = mvrnorm(n = n,mu = c(0,0,0,0,0,0),Sigma = CovR)
syn2 = data.frame(D0 = D0+c[,1], D1 = D1+c[,2])
syn2[,'t3'] = mapply(function(x,y){I1(on,B0=x,B1 = y)},x = syn2$D0,y = syn2$D1)
syn2$delta = syn2$t3-iT

## Generate random levels and slopes for each individual for Marker 2
d2 = data.frame(C0 = B0+c[,3], C1 = B1+c[,4])
#d2$lag1 = rgamma(n,shape = 10,rate = 2)
d2$lag1 = 5
syn2[,'t2'] = mapply(function(x,y,z,l){I2(on,B0=x,B1 = y,d = z,l = l )},x = d2$C0,y = d2$C1,z = syn2$delta,l = d2$lag1)
B2.t = function(t){A1(t-5,D0,D1)}


##S2: Generate random levels and slopes for each individual for Marker 3
## AMY is third, not independent in this case. 
iT3 = 5 #lag for marker 3
d3 = data.frame(B0 = D0+c[,5], B1 = D1+c[,6])
#d3$lag2 = runif(n =n, min = -5,max = 7)
d3$lag2 = 5
syn2[,'t1'] = mapply(function(x,y,z,l){I2(on,B0=x,B1 = y,d = z,l = l )},x = d3$B0,y = d3$B1,z = syn2$delta,l = d2$lag1+d3$lag2)
B3.t = function(t){A1(t-10,D0,D1)}

set2 = data.frame(syn2,d2,d3)
set2 = set2[,order(colnames(set2))]

####################################################
####################################################
## Data ##
#S3
# Random Start levels for S3 (Marker 2 initiates)
iT=10
C0 = qnorm(sL)
C1 = (qnorm(on)-qnorm(sL))/iT
B1.t = function(t){A1(t,B0=C0,B1=C1)}
##S2: Generate random levels and slopes for each individual for Marker 3
d = mvrnorm(n = n,mu = c(0,0,0,0,0,0),Sigma = CovR)
syn3 = data.frame(C0 = C0+d[,1], C1 = C1+d[,2])
syn3[,'t2'] = mapply(function(x,y){I1(on,B0=x,B1 = y)},x = syn3$C0,y = syn3$C1)
syn3$delta = syn3$t2-iT

## Generate random levels and slopes for each individual for Marker 1
d1 = data.frame(B0 = C0+c[,3], B1 = C1+c[,4])
d1$lag1 = 5
syn3[,'t1'] = mapply(function(x,y,z,l){I2(on,B0=x,B1 = y,d = z,l = l )},x = d1$B0,y = d1$B1,z = syn3$delta,l = d1$lag1)
B2.t = function(t){A1(t-5,C0,C1)}


##S2: Generate random levels and slopes for each individual for Marker 3
## AMY is third, not independent in this case. 
iT3 = 5 #lag for marker 3
d3 = data.frame(D0 = C0+c[,5], D1 = C1+c[,6])
#d3$lag2 = runif(n =n, min = -5,max = 7)
d3$lag2 = 5
syn3[,'t3'] = mapply(function(x,y,z,l){I2(on,B0=x,B1 = y,d = z,l = l )},x = d3$D0,y = d3$D1,z = syn3$delta,l = d1$lag1+d3$lag2)
B3.t = function(t){A1(t-10,D0,D1)}

set3 = data.frame(syn3,d1,d3)
set3 = set3[,order(colnames(set3))]



## Save out Set 1, Set 2 data sets ##
saveRDS(set1,file = 'set1a.RDA')
saveRDS(set2,file = 'set2a.RDA')
saveRDS(set3,file = 'set3a.RDA')

