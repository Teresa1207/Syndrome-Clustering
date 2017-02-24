################################################################################
######    Syndrome Classifiation  		
######	  By: Teresa Filshtein
######	  Date: April 22, 2016			                      
######	  Modified: August 22, 2016
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
av45var = c('RID','viscode2','base3dx','dx3cat',"av45date","lag_av45","suvr")
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

### AV45 ###
av45 = data.frame(databreak[9])
av45 = av45[!is.na(av45$lag_av45),]
av45Dx = xtabs(~av45$dx3cat)
av45$dx3cat[which(av45$dx3cat=='')] = av45$base3dx[which(av45$dx3cat=='')] 
av45Dx = xtabs(~av45$dx3cat)
av45$weight = with(av45, ifelse(dx3cat=='1:Normal',dxRates[1]/av45Dx[1],
                              ifelse(dx3cat=='2:MCI',dxRates[2]/av45Dx[2],dxRates[3]/av45Dx[3])))



## ECDF ##
hoc.w = with(hoc,data.frame(wtd.Ecdf(hoc,weights=weight)))
hoc.nw = with(hoc,data.frame(Ecdf(hoc,what = '1-F')))
hoc.w$S = 1-hoc.w$ecdf
colnames(hoc.w) = c('hoc','hocF','hocQ')
plot(hoc.w$hoc,hoc.w$hocQ,type = 'l')

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


## ECDF ##
av45.w = with(av45,data.frame(wtd.Ecdf(suvr,weights=weight)))
colnames(av45.w) = c('av45','av45Q')
plot(av45.w$av45,av45.w$av45Q,type = 'l')

####
hocM = merge(hoc, hoc.w[-1,c('hoc','hocQ')], by.x = 'hoc',by.y = 'hoc',all.x = TRUE)
adasM = merge(adas, adas.w[-1,], by.x = 'TOTAL13',by.y = 'adas',all.x = TRUE)
fdgM = merge(fdg, fdg.w[-1,c('fdg','fdgQ')], by.x = 'rall',by.y = 'fdg',all.x = TRUE)
av45M = merge(av45, av45.w[-1,], by.x = 'suvr',by.y = 'av45',all.x = TRUE)

hocM = with(hocM,hocM[order(RID,lag_fs51),c('RID','hoc','hocQ','lag_fs51')])
adasM = with(adasM,adasM[order(RID,lag_adas),c('RID','TOTAL13','adasQ','lag_adas')])
fdgM = with(fdgM,fdgM[order(RID,lag_fdg),c('RID','rall','fdgQ','lag_fdg')])
av45M = with(av45M,av45M[order(RID,lag_av45),c('RID','suvr','av45Q','lag_av45')])

hocM$type = 'hoc'
adasM$type = 'adas'
fdgM$type = 'fdg'
av45M$type = 'av45'

colnames(fdgM)=colnames(hocM)=colnames(adasM) =colnames(av45M)= c('RID','raw','quant','time','type')




### Create frame of variables ###
markers = rbind(hocM,adasM,fdgM,av45M)
markers = with(markers,markers[order(RID,type,time),])
markers$RID = factor(markers$RID)
## Levels and Slopes ##
adasU = subset(markers,type=='adas'&!is.na(time))
hocU = subset(markers,type=='hoc'&!is.na(time))
fdgU = subset(markers,type=='fdg'&!is.na(time))
av45U = subset(markers,type=='av45'&!is.na(time))
cA = subset(count(adasU,vars = "RID"),freq>1)
cH = subset(count(hocU,vars = "RID"),freq>1)
cF = subset(count(fdgU,vars = "RID"),freq>1)
cAv = subset(count(av45U,vars = "RID"),freq>1)

cross = intersect(cA$RID,cH$RID)
cross = intersect(cross,cF$RID)
cross = intersect(cross,cAv$RID)
mark.long = subset(markers,RID%in%cross)

dc.u = dApril[!duplicated(dc$RID),c('RID','base3dx','age65','apoe4','female')]
mark.long = merge(mark.long,dc.u,by = 'RID',all.x = TRUE)
