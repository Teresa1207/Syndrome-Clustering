################################################################################
######    Syndrome Classifiation  		
######	  By: Teresa Filshtein
######	  Date: April 22, 2016			                      
######	  Modified: 
######    EDITED BY: 
################################################################################
options(stringsAsFactors = FALSE)
library(date)
library(plyr)
library(xtable)
library(corrplot)
library(Hmisc)
library(reshape2)
library(GGally)
setwd('~/Documents/Davis2015-2016/ADNIDanielle/Data/DataApril')
dApril = read.csv("alldata_go2_20160419.csv",header=T)
nlids=with(dApril, unique(RID[base3dx=='1:Normal']))
mciids=with(dApril, unique(RID[base3dx=='2:MCI']))
adids=with(dApril, unique(RID[base3dx=='3:AD']))
`%ni%` <- Negate(`%in%`) 
################################################################################
setwd("~/Documents/Dissertation/R_Dissertation/Syndrome")
dc = dApril[,c('RID','viscode2','base3dx',"dx3cat","age65","fs51date",'qcvar',"lag_fs51",'hoc','hoc0',"hpcv","etrv","etrt","ventricles","wbrain",
                 "wmhdate","lag_wmh","WHITMATHYP","fmridate","lag_fmri","PDMNRV","asldate","lag_asl","precuneus_asl","postcing_asl",
                 "dtidate","lag_dti","cing_hipp_dti","bsi3dt","lag_bsi","kmndbcbbsi3","vbsi3","hbsi3", "tbmdate","lag_tbm",
                 "TBMSYNSCOR","fdgdate","lag_fdg","rall","av45date","lag_av45","suvr","lag_adas","TOTAL13","lag_mmse","MMSCORE","lag_ravlt","RAVLTtotal","lag_cdr","cdrsum","lag_ecogpt","ecog_score",
                 "lag_ecogsp","specog_score")]
mrivar1 = c('RID','viscode2','base3dx','dx3cat',"fs51date",'qcvar',"lag_fs51",'hoc',"hpcv","etrv","etrt","ventricles","wbrain")
mrivar2 = c('RID','viscode2','base3dx',"wmhdate",'qcvar',"lag_wmh","WHITMATHYP")
mrivar3 = c('RID','viscode2','base3dx',"bsi3dt","lag_bsi","kmndbcbbsi3","vbsi3","hbsi3")
mrivar4 = c('RID','viscode2','base3dx',"tbmdate","lag_tbm",
            "TBMSYNSCOR")
fmrivar = c('RID','viscode2','base3dx',"lag_fmri","PDMNRV")
dtivar = c('RID','viscode2','base3dx',"dtidate","lag_dti","cing_hipp_dti")
aslvar = c('RID','viscode2','base3dx',"asldate","lag_asl","precuneus_asl","postcing_asl")
fdgvar =c('RID','viscode2','base3dx','dx3cat',"lag_fdg","rall")
av45var = c('RID','viscode2','base3dx',"av45date","lag_av45","suvr")
cogvar = c('RID','viscode2','base3dx','dx3cat',"lag_adas","lag_mmse","lag_ravlt","lag_cdr","lag_ecogpt","lag_ecogsp",
           "TOTAL13","MMSCORE","RAVLTtotal","cdrsum","ecog_score","specog_score")

types = list(mrivar1,mrivar2,mrivar3,mrivar4,fmrivar,dtivar,aslvar,fdgvar,av45var,cogvar)
databreak = sapply(types,function(x){
  dc[,x]
})
dxRates = c(.5,.3,.2)

#### HOC ####
hoc = data.frame(databreak[1])
hoc = hoc[!is.na(hoc$hoc),]
hocDx = xtabs(~hoc$dx3cat)
hoc$dx3cat[which(hoc$dx3cat=='')] = hoc$base3dx[which(hoc$dx3cat=='')] 
hocDx = xtabs(~hoc$dx3cat)
hoc$weight = ifelse(hoc$dx3cat=='1:Normal',dxRates[1]/hocDx[1],ifelse(hoc$dx3cat=='2:MCI',dxRates[2]/hocDx[2],dxRates[3]/hocDx[3]))
## ECDF ##
hoc.w = with(hoc,data.frame(wtd.Ecdf(hoc,weights=weight)))
hoc.nw = with(hoc,data.frame(Ecdf(hoc,what = '1-F')))
hoc.w$S = 1-hoc.w$ecdf
colnames(hoc.w) = c('hoc','hocF','hocQ')
plot(hoc.w$hoc,hoc.w$hocQ,type = 'l')

### ADAS ###
adas = data.frame(databreak[10])
adas = adas[!is.na(adas$TOTAL13),]
adasDx = xtabs(~adas$dx3cat)
adas$dx3cat[which(adas$dx3cat=='')] = adas$base3dx[which(adas$dx3cat=='')] 
adasDx = xtabs(~adas$dx3cat)
adas$weight = with(adas, ifelse(dx3cat=='1:Normal',dxRates[1]/adasDx[1],
                                ifelse(dx3cat=='2:MCI',dxRates[2]/adasDx[2],dxRates[3]/adasDx[3])))
## ECDF ##
adas.w = with(adas,data.frame(wtd.Ecdf(TOTAL13,weights=weight)))
#adas.nw = with(hoc,data.frame(Ecdf(hoc,what = '1-F')))
colnames(adas.w) = c('adas','adasQ')
plot(adas.w$adas,adas.w$adasQ,type = 'l')

### FDG ###
fdg = data.frame(databreak[8])
fdg = fdg[!is.na(fdg$lag_fdg),]
fdgDx = xtabs(~fdg$dx3cat)
fdg$dx3cat[which(fdg$dx3cat=='')] = fdg$base3dx[which(fdg$dx3cat=='')] 
fdgDx = xtabs(~fdg$dx3cat)
fdg$weight = with(fdg, ifelse(dx3cat=='1:Normal',dxRates[1]/fdgDx[1],
                                ifelse(dx3cat=='2:MCI',dxRates[2]/fdgDx[2],dxRates[3]/fdgDx[3])))

## ECDF ##
fdg.w = with(fdg,data.frame(wtd.Ecdf(rall,weights=weight)))
fdg.nw = with(fdg,data.frame(Ecdf(rall,what = '1-F')))
fdg.w$S = 1-fdg.w$ecdf
colnames(fdg.w) = c('fdg','fdgF','fdgQ')
plot(fdg.w$fdg,fdg.w$fdgQ,type = 'l')

## ECDF ##
adas.w = with(adas,data.frame(wtd.Ecdf(TOTAL13,weights=weight)))
#adas.nw = with(hoc,data.frame(Ecdf(hoc,what = '1-F')))
colnames(adas.w) = c('adas','adasQ')
plot(adas.w$adas,adas.w$adasQ,type = 'l')

####
hocM = merge(hoc, hoc.w[-1,c('hoc','hocQ')], by.x = 'hoc',by.y = 'hoc',all.x = TRUE)
adasM = merge(adas, adas.w[-1,], by.x = 'TOTAL13',by.y = 'adas',all.x = TRUE)
fdgM = merge(fdg, fdg.w[-1,c('fdg','fdgQ')], by.x = 'rall',by.y = 'fdg',all.x = TRUE)

hocM = with(hocM,hocM[order(RID,lag_fs51),c('RID','hoc','hocQ','lag_fs51')])
adasM = with(adasM,adasM[order(RID,lag_adas),c('RID','TOTAL13','adasQ','lag_adas')])
fdgM = with(fdgM,fdgM[order(RID,lag_fdg),c('RID','rall','fdgQ','lag_fdg')])

hocM$type = 'hoc'
adasM$type = 'adas'
fdgM$type = 'fdg'
colnames(fdgM)=colnames(hocM)=colnames(adasM) = c('RID','raw','quant','time','type')




### Create frame of variables ###
markers = rbind(hocM,adasM,fdgM)
markers = with(markers,markers[order(RID,type,time),])
markers$RID = factor(markers$RID)
## Levels and Slopes ##
model <- function(x){
  fit1 = try(lm(raw ~ time, data=x))
  fit2 = try(lm(quant ~ time, data=x))
  data.frame(rL=coef(fit1)[[1]], rS=coef(fit1)[[2]],
             qL=coef(fit2)[[1]], qS=coef(fit2)[[2]])
  }
model1 <- function(x){
  fit1 = try(lm(raw ~ time, data=x))
  fit2 = try(lm(quant ~ time, data=x))
  data.frame(rL=coef(fit1)[[1]], rS=coef(fit1)[[2]],qL=coef(fit2)[[1]], qS=coef(fit2)[[2]])
}
adasU = subset(markers,type=='adas'&!is.na(time))
hocU = subset(markers,type=='hoc'&!is.na(time))
fdgU = subset(markers,type=='fdg'&!is.na(time))
cA = subset(count(adasU,vars = "RID"),freq>1)
cH = subset(count(hocU,vars = "RID"),freq>1)
cF = subset(count(fdgU,vars = "RID"),freq>1)
cross = intersect(cA$RID,cH$RID)
cross = intersect(cross,cF$RID)
mark.long = subset(markers,RID%in%cross)

M1 = ddply(mark.long,.(RID,type),function(x){model(x)})

M1L = dcast(M1[,c("RID",'type',"qL","qS")],RID~type,value.var = c('qL'))
colnames(M1L) = c('RID','aL','fL','hL')
M1S = dcast(M1[,c("RID",'type',"qL","qS")],RID~type,value.var = c('qS'))
colnames(M1S) = c('RID','aS','fS','hS')
M1cast = merge(M1L,M1S,by = 'RID',all = TRUE)

mark.cast = dcast(markers[,c("RID",'type',"","qS")],RID~type,value.var = c('qL'))

#ggpairs(M1cast[,-1])

##Transform Betas##
trans = function(x){scale(x,scale=FALSE)}
M1cent = data.frame(M1cast$RID,t(data.frame(apply(M1cast[,c('aL','fL','hL')],1,trans))),
                    t(data.frame(apply(M1cast[,c('aS','fS','hS')],1,trans))))
colnames(M1cent) = c("RID",'aLt','fLt','hLt','aSt','fSt','hSt')


### plots ###
library(ggplot2)

ggplot(subset(mark.long,RID%in%cross[1:100]),aes(x = time,y = quant,group = RID,color=RID))+geom_line()+facet_grid(~type)

d1 = subset(M1cent)
k1 = data.frame(RID = d1$RID,clust = kmeans(d1[-1],3)$cluster)
m1 = merge(mark.long,k1,by = 'RID',all.y = TRUE)
m1 = merge(m1,unique(dc[,c("RID","age65","base3dx")]),by = 'RID',all.x = TRUE)

m1$age = m1$time+m1$age65
ggplot(m1,aes(x = age,y = quant,group = RID,color = factor(base3dx)))+geom_line()+
  stat_smooth(aes(group = 1),color = 'black')+facet_grid(clust~type) 
ggplot(m1,aes(x = time,y = quant,group = RID,color = factor(base3dx)))+geom_line()+
  stat_smooth(aes(group = 1),color = 'black')+facet_grid(clust~type) 

ggplot(m1,aes(x = age,y = quant,group = RID,color = factor(base3dx)))+geom_line()+
  stat_smooth(aes(group = 1),color = 'black')+facet_grid(base3dx~type) 

## Kendall's rank
kD = merge(d1,k1,by = 'RID')
levs = kD[,c('aLt','fLt','hLt')]
rownames(levs) = kD[,"RID"]

slps = kD[,c('aSt','fSt','hSt')]
rownames(kD) = kD[,"RID"]

levMat = lapply(1:max(kD$clust),function(x){cor(levs[kD$clust==x,],method = "kendall",use= 'pairwise')})


slpMat = cor(slps,method = "kendall",use= 'pairwise')
totMat = cbind(levMat,slpMat)
lM = t(apply(levs,1,rank))
sM = t(apply(slps,1,rank))
disL = matrix(0,nrow = nrow(kD),ncol = nrow(kD))
disS = matrix(0,nrow = nrow(kD),ncol = nrow(kD))
#A1 = outer(rep(1,nrow(kD)),1:(nrow(kD)))
cFun = function(df,x,y){cor(df[x,],df[y,],method='kendall')}

for(i in 1:nrow(kD)){
  for(j in 1:nrow(kD)){
    disL[i,j] = cFun(lM,i,j)
    disS[i,j] = cFun(sM,i,j)
  }
}

disT = disL+disS

disM1 = 2 - disT
disM1d = as.dist(disM1)

c1 = hclust(disM1d)
plot(c1)
ct<- cutree(c1,k=3)
ctd = data.frame(RID = kD$RID,ct = ct)
rect.hclust(c1, k = 3)
m1d = merge(m1,ctd,by = 'RID',all.x = TRUE)
ggplot(m1d,aes(x = age,y = quant,group = RID,color = factor(base3dx)))+geom_line()+
  stat_smooth(aes(group = 1),color = 'black')+facet_grid(clust~type) 
ggplot(m1d,aes(x = age,y = quant,group = RID,color = factor(base3dx)))+geom_line()+
  stat_smooth(aes(group = 1),color = 'black')+facet_grid(ct~type) 

