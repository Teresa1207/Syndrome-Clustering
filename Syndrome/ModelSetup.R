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
## Data ##
t = seq(-10,30,by = 1)
n = 1000
#S1  AMY: B0 = qnorm(.2), B1 = (qnorm(.5)-qnorm(.2))/10
#S1 HPCV: C0,C1
#S1 Vasc: D0,D1
# Random Start levels, mean B0
sL = .1 #start level
on = .5 #onset
iT = 10 #time for first marker
B0 = qnorm(sL)
B1 = (qnorm(on)-qnorm(sL))/iT
A1 = function(t,B0,B1){pnorm(B0+B1*t)}
I1 = function(t,B0,B1){(qnorm(t) - B0)/B1}
A1.t = function(t){A1(t,B0=B0,B1=B1)}
curve(A1.t,-10,30,n = 100)

## Generate random levels and slopes for each individual for Marker 1
b = mvrnorm(n = n,mu = c(0,0),Sigma = matrix(c(.01^2,0,0,.01^2),nrow=2,ncol=2))
data1 = data.frame(B0 = B0+b[,1], B1 = B1+b[,2])
data1[,'etp1'] = mapply(function(x,y){A1(iT,B0=x,B1 = y)},x = data1$B0,y = data1$B1)
data1[,'t1'] = mapply(function(x,y){I1(on,B0=x,B1 = y)},x = data1$B0,y = data1$B1)

## Generate random levels and slopes for each individual for Marker 2
# iT2 = lag of second marker
iT2 = 5
c = data.frame(c0 = B0+rnorm(n,0,.01^2),c1 = rgamma(n,shape = 10,rate = 2), e2 = rnorm(n,0,1^2))

data2 = with(c, data.frame(C0 = c0, C1 = -c0/(data1$t1+c1+e2)))
data1[,'etp2'] = mapply(function(x,y){A1(iT+iT2,B0=x,B1 = y)},x = data2$C0,y = data2$C1)
data1[,'t2'] = mapply(function(x,y){I1(on,B0=x,B1 = y)},x = data2$C0,y = data2$C1)
with(data1,hist(t2-t1))
with(data1,range(t2-t1))
A2.t = function(t){A1(t,B0,(qnorm(on)-qnorm(sL))/(iT+iT2))}

## Generate random levels and slopes for each individual for Marker 3
iT3 = 3 #lag for marker 3
d = data.frame(d0 = B0+rnorm(n,0,.01^2),d1 = rgamma(n,shape = 6,rate = 2), e3 = rnorm(n,0,1^2))

data3 = with(d, data.frame(D0 = d0, D1 = -d0/(data1$t2+d1+e3)))
data1[,'etp3'] = mapply(function(x,y){A1(iT+iT2+iT3,B0=x,B1 = y)},x = data3$D0,y = data3$D1)
data1[,'t3'] = mapply(function(x,y){I1(on,B0=x,B1 = y)},x = data3$D0,y = data3$D1)
with(data1,hist(t3-t2))
with(data1,range(t3-t2))
with(data1,hist(t3-t1))
with(data1,range(t3-t1))
with(data1,mean(t3 - t1))
A3.t = function(t){A1(t,B0,(qnorm(on)-qnorm(sL))/sum(c(iT,iT2,iT3)))}

set1 = data.frame(data1,data2,data3)
set1 = set1[,order(colnames(set1))]



## Data ##
#S2  AMY: B0 = qnorm(.2), B1 = (qnorm(.5)-qnorm(.2))/10
#S2 HPCV: C0,C1
#S2 Vasc: D0 = qnorm(.2), D1 = (qnorm(.5)-qnorm(.2))/10
# Random Start levels for S3 (Marker 3 initiates)
iT=10
D0 = qnorm(sL)
D1 = (qnorm(on)-qnorm(sL))/iT
A1 = function(t,B0,B1){pnorm(B0+B1*t)}
I1 = function(t,B0,B1){(qnorm(t) - B0)/B1}
B1.t = function(t){A1(t,B0=D0,B1=D1)}
curve(B1.t,-10,30,n = 100)
##S2: Generate random levels and slopes for each individual for Marker 3
b3 = mvrnorm(n = n,mu = c(0,0),Sigma = matrix(c(.01^2,0,0,.01^2),nrow=2,ncol=2))
data12 = data.frame(D0 = D0+b3[,1], D1 = D1+b3[,2])
data12[,'etp3'] = mapply(function(x,y){A1(iT,B0=x,B1 = y)},x = data12$D0,y = data12$D1)
data12[,'t3'] = mapply(function(x,y){I1(on,B0=x,B1 = y)},x = data12$D0,y = data12$D1)

## Generate random levels and slopes for each individual for Marker 2
c2 = data.frame(c0 = D0+rnorm(n,0,.01^2),c1 = rgamma(n,shape = 10,rate = 2), e2 = rnorm(n,0,1.5^2))
data22 = with(c2, data.frame(C0 = c0, C1 = -c0/(data12$t3+c1+e2)))

data12[,'etp2'] = mapply(function(x,y){A1(sum(iT,iT2),B0=x,B1 = y)},x = data22$C0,y = data22$C1)
data12[,'t2'] = mapply(function(x,y){I1(on,B0=x,B1 = y)},x = data22$C0,y = data22$C1)
with(data12,hist(t2-t3))
with(data12,range(t2-t3))
B2.t = function(t){A1(t,D0,(qnorm(on)-qnorm(sL))/sum(iT,iT2))}

##S3: Generate random levels and slopes for each individual for Marker 3
## Amy is totally independent, happens anywhere between 5 and 20
data23 = data.frame(B0 = D0+rnorm(n,0,.01^2), tA = runif(n = n, min = 5,max = 15))
data23$B1 = mapply(function(x,y){-x/y},x = data23$B0,y = data23$tA)
data23 = data23[,order(colnames(data23))]
data12[,'t1'] = mapply(function(x,y){I1(on,B0=x,B1 = y)},x = data23$B0,y = data23$B1)
B3.t = function(t){A1(t,B0,(qnorm(on)-qnorm(sL))/10)} #14 average time to onset for M1,S2

set2 = data.frame(data12,data22,data23)
set2 = set2[,order(colnames(set2))]


## Plot the two syndromes

## Save out Set 1, Set 2 data sets ##

saveRDS(set1,file = 'set1.RDA')
saveRDS(set2,file = 'set2.RDA')

#set1 = readRDS('set1.RDA')
#set2 = readRDS('set2.RDA')

### CURVES
dontrun = function(n){curve(A1.t,-10,30,n=n,col='red')
counter = 1
while(counter<=nrow(data1)){
  A1.p = function(t){A1(t,B0=data1[counter,1],B1=data1[counter,2])}
  curve(A1.p,-10,30,n=n,add = TRUE)
  counter = counter+1
}
curve(A1.t,-10,30,n=n,col='red',add = TRUE)


curve(A2.t,-10,30,n=n,col='red')
counter = 1
while(counter<=nrow(data1)){
  A1.p = function(t){A1(t,B0=data1[counter,1],B1=data1[counter,2])}
  A2.p = function(t){A1(t,B0=data2[counter,1],B1=data2[counter,2])}
  curve(A1.p,-10,30,n=n,add = TRUE)
  curve(A2.p,-10,30,n = n, add = TRUE, col = 'darkgreen')
  counter = counter+1
}
curve(A1.t,-10,30,n=n,col='red',add = TRUE,lwd = 2)
curve(A2.t,-10,30,n=n,col='red',add = TRUE,lwd = 2)

curve(A1.t,-10,30,n=n)
counter = 1
while(counter<=nrow(data1)){
  A1.p = function(t){A1(t,B0=data1[counter,1],B1=data1[counter,2])}
  A2.p = function(t){A1(t,B0=data2[counter,1],B1=data2[counter,2])}
  A3.p = function(t){A1(t,B0=data3[counter,1],B1=data3[counter,2])}
  curve(A1.p,-10,30,n=n,add = TRUE,col = 'black')
  curve(A2.p,-10,30,n = n, add = TRUE, col = 'darkgreen')
  curve(A3.p,-10,30,n = n, add = TRUE, col = 'darkblue')
  counter = counter+1
}
curve(A1.t,-10,30,n=n,col='red',add = TRUE,lwd = 3)
curve(A2.t,-10,30,n=n,col='red',add = TRUE,lwd = 3)
curve(A3.t,-10,30,n=n,col='red',add = TRUE,lwd = 3)


curve(B1.t,-10,30,n=n,col='red')
counter = 1
while(counter<=nrow(data1)){
  B1.p = function(t){A1(t,B0=data12[counter,1],B1=data12[counter,2])}
  curve(B1.p,-10,30,n=n,add = TRUE,col = 'darkblue')
  counter = counter+1
}
curve(B1.t,-10,30,n=n,col='red',add = TRUE,lwd=3)

curve(B2.t,-10,30,n=n,col='red')
counter = 1
while(counter<=nrow(data12)){
  B1.p = function(t){A1(t,B0=data12[counter,1],B1=data12[counter,2])}
  B2.p = function(t){A1(t,B0=data22[counter,1],B1=data22[counter,2])}
  curve(B1.p,-10,30,n=n,add = TRUE,col = 'darkblue')
  curve(B2.p,-10,30,n = n, add = TRUE, col = 'darkgreen')
  counter = counter+1
}
curve(B1.t,-10,30,n=n,col='red',add = TRUE,lwd = 2)
curve(B2.t,-10,30,n=n,col='red',add = TRUE,lwd = 2)

frm = -20
to = 50
par(mfrow = c(1,2))
curve(A1.t,frm,to,n=n)
counter = 1
while(counter<=nrow(data1)){
  A1.p = function(t){A1(t,B0=data1[counter,1],B1=data1[counter,2])}
  A2.p = function(t){A1(t,B0=data2[counter,1],B1=data2[counter,2])}
  A3.p = function(t){A1(t,B0=data3[counter,1],B1=data3[counter,2])}
  curve(A1.p,frm,to,n=n,add = TRUE,col = 'black')
  curve(A2.p,frm,to,n = n, add = TRUE, col = 'darkgreen')
  curve(A3.p,frm,to,n = n, add = TRUE, col = 'darkblue')
  counter = counter+1
}
curve(A1.t,frm,to,n=n,col='red',add = TRUE,lwd = 3)
curve(A2.t,frm,to,n=n,col='red',add = TRUE,lwd = 3,lty=5)
curve(A3.t,frm,to,n=n,col='red',add = TRUE,lwd = 3,lty=4)
#
curve(B3.t,frm,to,n=n,col='red')
counter = 1
while(counter<=nrow(data12)){
  B1.p = function(t){A1(t,B0=data12[counter,1],B1=data12[counter,2])}
  B2.p = function(t){A1(t,B0=data22[counter,1],B1=data22[counter,2])}
  B3.p = function(t){A1(t,B0=data23[counter,1],B1=data23[counter,2])}
  
  curve(B1.p,frm,to,n=n,add = TRUE,col = 'darkblue')
  curve(B2.p,frm,to,n = n, add = TRUE, col = 'darkgreen')
  curve(B3.p,frm,to,n = n, add = TRUE, col = 'black')
  counter = counter+1
}
curve(B1.t,frm,to,n=n,col='red',add = TRUE,lwd = 3,lty = 4)
curve(B2.t,frm,to,n=n,col='red',add = TRUE,lwd = 3,lty = 5)
curve(B3.t,frm,to,n=n,col='red',add = TRUE,lwd = 3)

curve(B3.t,-10,30,n=n,col='red')
counter = 1
while(counter<=nrow(data12)){
  B1.p = function(t){A1(t,B0=data12[counter,1],B1=data12[counter,2])}
  B2.p = function(t){A1(t,B0=data22[counter,1],B1=data22[counter,2])}
  B3.p = function(t){A1(t,B0=data23[counter,1],B1=data23[counter,2])}
  
  curve(B1.p,-10,30,n=n,add = TRUE,col = 'darkblue')
  curve(B2.p,-10,30,n = n, add = TRUE, col = 'darkgreen')
  curve(B3.p,-10,30,n = n, add = TRUE, col = 'black')
  
  counter = counter+1
}
curve(B1.t,-10,30,n=n,col='red',add = TRUE,lwd = 3)
curve(B2.t,-10,30,n=n,col='red',add = TRUE,lwd = 3)
curve(B3.t,-10,30,n=n,col='red',add = TRUE,lwd = 3)
}
