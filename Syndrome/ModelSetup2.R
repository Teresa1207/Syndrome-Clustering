### Syndrome Classification    ###
### Model Setup and Simulation ###
setwd("/Users/Teresa/Documents/Dissertation/R_Dissertation/Syndrome")
library(cluster)
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
## Data ##
t = seq(-10,30,by = 1)
n = 1000
#S1  AMY: B0 = qnorm(.2), B1 = (qnorm(.5)-qnorm(.2))/10
#S1 HPCV: C0,C1
#S1 Vasc: D0,D1
# Random Start levels, mean B0
sig = .01
sig1 = 0.01
Corr = matrix(.5,nrow = 6,ncol = 6)
diag(Corr) =1
CovR = diag(rep(c(sig,sig1),3))%*%Corr%*%diag(rep(c(sig,sig1),3))
sL = .1 #start level
on = .5 #onset
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
data2$lag1 = rgamma(n,shape = 10,rate = 2)
#syn1[,'etp2'] = mapply(function(x,y){A1(iT+iT2,B0=x,B1 = y)},x = data2$C0,y = data2$C1)
syn1[,'t2'] = mapply(function(x,y,z,l){I2(on,B0=x,B1 = y,d = z,l = l )},x = data2$C0,y = data2$C1,z = syn1$delta,l = data2$lag1)

A2.t = function(t){A1(t-5,B0,B1)}

## Generate random levels and slopes for each individual for Marker 3
iT3 = 3 #lag for marker 3
data3 = data.frame(D0 = B0+b[,5], D1 = B1+b[,6])
data3$lag2 = rgamma(n,shape = 6,rate = 2)
#syn1[,'etp3'] = mapply(function(x,y){A1(iT+iT2+iT3,B0=x,B1 = y)},x = data3$D0,y = data3$D1)
syn1[,'t3'] = mapply(function(x,y,z,l){I2(on,B0=x,B1 = y,d = z,l = l )},x = data3$D0,y = data3$D1,z = syn1$delta,l = (data3$lag2+data2$lag1))
A3.t = function(t){A1(t-8,B0,B1)}

set1 = data.frame(syn1,data2,data3)
set1 = set1[,order(colnames(set1))]


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
d2$lag1 = rgamma(n,shape = 10,rate = 2)
syn2[,'t2'] = mapply(function(x,y,z,l){I2(on,B0=x,B1 = y,d = z,l = l )},x = d2$C0,y = d2$C1,z = syn2$delta,l = d2$lag1)
B2.t = function(t){A1(t-5,D0,D1)}


##S3: Generate random levels and slopes for each individual for Marker 3
## Amy is totally independent, happens anywhere between 5 and 20
iT3 = 3 #lag for marker 3
d3 = data.frame(B0 = D0+c[,5], B1 = D1+c[,6])
d3$lag2 = runif(n =n, min = -5,max = 7)
syn2[,'t1'] = mapply(function(x,y,z,l){I2(on,B0=x,B1 = y,d = z,l = l )},x = d3$B0,y = d3$B1,z = rep(0,n),l = d3$lag2)
B3.t = function(t){A1(t-1,D0,D1)}

set2 = data.frame(syn2,d2,d3)
set2 = set2[,order(colnames(set2))]


## Plot the two syndromes
dontrun2= function(d){
with(syn1,hist(t2-t1,main = 'Syndrome 1, Marker 1,2, Onset Differential'))
with(syn1,hist(t3-t2,main = 'Syndrome 1, Marker 2,3, Onset Differential'))
with(syn1,hist(t3-t1,main = 'Syndrome 1, Marker 1,3, Onset Differential'))
with(syn2,hist(t2-t3,main = 'Syndrome 2, Marker 3,2, Onset Differential'))
with(syn2,hist(t1-t3,main = 'Syndrome 2, Marker 3,1, Onset Differential'))
with(syn2,hist(t2-t1,main = 'Syndrome 2, Marker 2,1, Onset Differential'))
## Save out Set 1, Set 2 data sets ##
}
saveRDS(set1,file = 'set1a.RDA')
saveRDS(set2,file = 'set2a.RDA')

#set1 = readRDS('set1.RDA')
#set2 = readRDS('set2.RDA')

### CURVES
dontrun = function(n){
  curve(A1.t,-10,30,n=n,col='red')
  counter = 1
  while(counter<=nrow(syn1)){
    A1.p = function(t){A1((t+syn1$delta[counter]),B0=syn1[counter,1],B1=syn1[counter,2])}
    A2.p = function(t){A1((t+syn1$delta[counter]-data2$lag1[counter]),B0=data2[counter,1],B1=data2[counter,2])}
    A3.p = function(t){A1((t+syn1$delta[counter]-data2$lag1[counter]-data3$lag2[counter]),B0=data3[counter,1],B1=data3[counter,2])}
    curve(A1.p,-10,30,n=n,add = TRUE)
    curve(A2.p,-10,30,n = n, add = TRUE, col = 'darkgreen')
    curve(A3.p,-10,30,n = n, add = TRUE, col = 'darkblue')
    counter = counter+1
  }
  curve(A1.t,-10,30,n=n,col='red',add = TRUE,lwd = 2)
  curve(A2.t,-10,30,n=n,col='red',add = TRUE,lwd = 2)
  curve(A3.t,-10,30,n=n,col='red',add = TRUE,lwd = 2)
  
  curve(B1.t,-10,30,n=n,col='red')
  counter = 1
  while(counter<=nrow(syn2)){
    B1.p = function(t){A1(t+syn2$delta[counter],B0=syn2[counter,1],B1=syn2[counter,2])}
    B2.p = function(t){A1((t+syn2$delta[counter]-d2$lag1[counter]),B0=d2[counter,1],B1=d2[counter,2])}
    B3.p = function(t){A1((t+d3$lag2[counter]),B0=d3[counter,1],B1=d3[counter,2])}
    
    curve(B1.p,-10,30,n=n,add = TRUE,col = 'darkblue')
    curve(B2.p,-10,30,n=n,add = TRUE,col = 'darkgreen')
    curve(B3.p,-10,30,n=n,add = TRUE,col = 'black')
    counter = counter+1
  }
  curve(B1.t,-10,30,n=n,col='red',add = TRUE,lwd=3)
  curve(B2.t,-10,30,n=n,col='red',add = TRUE,lwd=3)
  curve(B3.t,-10,30,n=n,col='red',add = TRUE,lwd=3)
  
 }
