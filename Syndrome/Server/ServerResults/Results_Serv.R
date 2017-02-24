## Server Results
setwd("~/Documents/Dissertation/R_Dissertation/Syndrome/Server/ServerResults")

#scenarios = c(111,121,131,141,
#              211,221,231,241,
#              311,321,331,341,
#              112,113,114,115) STILL NEED THESE
#poisson: 121,141
#gumbel: ,
#gcox: ,321,331,341

scenDat = data.frame(Serv = c(111,121,131,141,
                              211,221,231,241,
                              311,321,331,341,
                              112,113,114,115,
                              116,117,118,238,338),Scen = c(1:21))

scenarios = scenDat$Serv
resServ = lapply(scenarios,function(x){
  readRDS(paste('scenario_',x,'.RDA',sep = ''))})

names(resServ) = unlist(lapply(scenarios,function(x){
    paste('res',scenDat$Scen[which(scenDat$Serv==x)],sep='')}))

## Data Summary
clustNoT = function(x){
  lapply(2:6,function(n){apply(x[,1:4],2,function(l){length(which(l==n))})})
}

clustNo= function(x,n){
  apply(x[,1:4],2,function(l){length(which(l==n))/nrow(x)})
}


RImean = function(x){
  apply(x[,5:9],2,mean)
}


ResTable = data.frame(t(sapply(1:(length(scenarios)),function(x){
  clustNo(data.frame(resServ[x]),3)})),t(sapply(1:(length(scenarios)),function(x){
    RImean(data.frame(resServ[x]))}))
  )
rownames(ResTable) =  unlist(lapply(scenarios,function(x){
  paste('Scenario ',scenDat$Scen[which(scenDat$Serv==x)],sep='')}))

ResTable[,1]=ResTable[,1]*100
ResTable = round(ResTable,digits = c(rep(2,4),rep(4,5)))
