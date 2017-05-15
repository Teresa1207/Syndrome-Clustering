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
library(ggplot2)
library(Rmisc)
library(gridExtra)
library(grid)
library(geepack)
#setwd('~/Documents/Davis2015-2016/ADNIDanielle/Data/DataApril')
setwd("~/Documents/Dissertation/R_Dissertation/Syndrome")
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
#dxRates = c(.7,.2,.1)

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
av45.r = subset(av45, RID!='4386') #duplicate and huge outlier, prob not real 
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

cA = subset(as.data.frame(table(adasU$RID)),Freq>1)
cH = subset(as.data.frame(table(hocU$RID)),Freq>1)
cF = subset(as.data.frame(table(fdgU$RID)),Freq>1)
cAv = subset(as.data.frame(table(av45U$RID)),Freq>1)

cross = intersect(cA$Var1,cH$Var1)
cross = intersect(cross,cF$Var1)
cross = intersect(cross,cAv$Var1)
mark.long = subset(markers,RID%in%cross)

dc.u = dApril[!duplicated(dc$RID),c('RID','base3dx','age65','apoe4','female')]
mark.long = merge(mark.long,dc.u,by = 'RID',all.x = TRUE)


### Interesting plots
themeP =theme_bw()+
        theme(legend.title=element_blank(), legend.text = element_text(size = 30,face = 'bold'),
        legend.position = 'bottom',
        axis.text.x = element_text(size = 24,face = 'bold'),
        axis.text.y = element_text(size = 24,face = 'bold'),
        axis.title = element_text(size = 24,face = 'bold'))

  

wP = .03
hP = .03
# The palette with black:
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
colP = cbbPalette[c(2,4,6)]
fdgP = merge(fdgM, fdg[,c('RID',"base3dx","dx3cat")],by = 'RID')
hocP = merge(hocM, hoc[,c('RID',"base3dx","dx3cat")],by = 'RID')
av45P = merge(av45M, av45.r[,c('RID',"base3dx","dx3cat")],by = 'RID')
av45P = subset(av45P, RID!='4386') #duplicate and huge outlier, prob not real 

cdf1 = ggplot(fdgP,aes(x = raw, y = quant))+geom_jitter(aes(color = dx3cat), size = 2, width = wP, height = hP)+
  xlab('FDG')+ylab('Prop > x')+
  scale_colour_manual(values=colP)+
  geom_line(size = 2) + coord_cartesian(ylim = c(0,1))+
  themeP +
  guides(color = guide_legend(override.aes = list(size=4)))

cdf2 = ggplot(hocP,aes(x = raw, y = quant,group = dx3cat))+geom_jitter(aes(color = dx3cat), size = 2, width = wP, height = hP)+
  xlab('HOC')+ylab('Prop > x')+
    scale_colour_manual(values=colP)+
    geom_line(size = 2) + coord_cartesian(ylim = c(0,1))+
  themeP +
  guides(color = guide_legend(override.aes = list(size=4)))
cdf3 = ggplot(av45P,aes(x = raw, y = quant))+geom_jitter(aes(color = base3dx), size = 2, width = wP, height = hP)+
  xlab('AV45')+ylab('Prop > x')+
  scale_colour_manual(values=colP)+
  geom_line(size = 2) + coord_cartesian(ylim = c(0,1))+
themeP +
  guides(color = guide_legend(override.aes = list(size=4)))

mylegend6<-g_legend(cdf1)



png('~/Documents/Dissertation/Dissertation_2017/figures/cdfPlots.png',height = 1000, width = 1500)

gCDF =  grid.arrange(arrangeGrob(cdf3 + theme(legend.position="none"),
                                 cdf1 + theme(legend.position="none"),#+ylab(""),
                                 cdf2 + theme(legend.position="none"),#+ylab(""),
                                 nrow=3),
                     mylegend6,nrow = 2,heights=c(5,.5))
dev.off()
#### Tables
library(dplyr)
tabav = tbl_df(av45P[!duplicated(av45P$RID),])
tabav %>% group_by(base3dx) %>%
  summarise(
    avsum  = sum(quant)
    #av45 = length(RID)
  )


tabfd = tbl_df(fdgP[!duplicated(fdgP$RID),])%>% group_by(base3dx) %>%
  summarise(
    FDG = length(RID)
  )

tabho = tbl_df(hocP[!duplicated(hocP$RID),])%>% group_by(base3dx) %>%
  summarise(
    HOC = length(RID)
  )

tabCDF  = merge(tabav,tabfd, by = 'base3dx')
tabCDF = merge(tabCDF,tabho, by = 'base3dx')

print(xtable(tabCDF),include.rownames = FALSE)

##### STOP ####

g1=readRDS('~/Documents/Dissertation/R_Dissertation/Syndrome/g2T1.RDA')
g2=readRDS('~/Documents/Dissertation/R_Dissertation/Syndrome/g2T2.RDA')

h1=readRDS('~/Documents/Dissertation/R_Dissertation/Syndrome/g3T1.RDA')
h2=readRDS('~/Documents/Dissertation/R_Dissertation/Syndrome/g3T2.RDA')
h3=readRDS('~/Documents/Dissertation/R_Dissertation/Syndrome/g3T3.RDA')

clustersADNI2 = data.frame(RID = c(as.character(g1$fits$adas$subset$id),
                                     as.character(g2$fits$adas$subset$id)),
                             group2 = c(as.character(g1$fits$adas$subset$group),
                                       as.character(g2$fits$adas$subset$group)))
clustersADNI2 = clustersADNI2[!duplicated(clustersADNI2$RID),]
clustersADNI3 = data.frame(RID = c(as.character(h1$fits$adas$subset$id),
                                   as.character(h2$fits$adas$subset$id),
                                   as.character(h3$fits$adas$subset$id)),
                           group3 = c(as.character(h1$fits$adas$subset$group),
                                     as.character(h2$fits$adas$subset$group),
                                     as.character(h3$fits$adas$subset$group)))
clustersADNI3 = clustersADNI3[!duplicated(clustersADNI3$RID),]

markClus = merge(mark.long,clustersADNI2,by = 'RID')
markClus = merge(markClus,clustersADNI3,by = 'RID')

clusAll = subset(markClus,type!='adas')
clusMCI = subset(markClus,base3dx=='2:MCI'&type!='adas')
clusCN = subset(markClus,base3dx=='1:Normal'&type!='adas')

dfClus = clusMCI[order(clusMCI$RID,clusMCI$type,clusMCI$time),]
#dfClus = clusAll[order(clusAll$RID,clusAll$type,clusAll$time),]
#dfClus = clusCN[order(clusCN$RID,clusCN$type,clusCN$time),]
library(nlme)
library(geepack)
geeG = geeglm(quant~I(time+age65)*type,family = binomial(),id = RID,data = markClus)
geeG2a = function(x){geeglm(quant~I(time+age65)*type,family = binomial(link = 'probit'),id = RID,corstr="ar1",data = subset(dfClus,group2==x))}
geeG3a = function(x){geeglm(quant~I(time+age65)*type,family = binomial(link = 'probit'),id = RID,corstr="ar1",data = subset(dfClus,group3==x))}
geeG2 = function(x){geeglm(quant~I(time+age65)*type,family = binomial(link = 'probit'),id = RID,data = subset(dfClus,group2==x))}
geeG3 = function(x){geeglm(quant~I(time+age65)*type,family = binomial(link = 'probit'),id = RID,data = subset(dfClus,group3==x))}




fun1 = function(f,time,type){
  if(type=='av45'){
  eta = coef(f)[1]+coef(f)[2]*time}else{
    if(type=='fdg'){
      eta = coef(f)[1]+coef(f)[3]+(coef(f)[2]+coef(f)[5])*time}else{
        if(type=='hoc'){
          eta = coef(f)[1]+coef(f)[4]+(coef(f)[2]+coef(f)[6])*time}else{
    }
  }}
  pnorm(eta)
}

# fun1 = function(f,time,type){
#   if(type=='av45'){
#     eta = fixed.effects(f)[1]+fixed.effects(f)[2]*time}else{
#       if(type=='fdg'){
#         eta = fixed.effects(f)[1]+fixed.effects(f)[3]+(fixed.effects(f)[2]+fixed.effects(f)[5])*time}else{
#           if(type=='hoc'){
#             eta = fixed.effects(f)[1]+fixed.effects(f)[4]+(fixed.effects(f)[2]+fixed.effects(f)[6])*time}else{
#
#             }
#         }}
#   pnorm(eta)
# }

#geePlot = function(x,t){fun1(geeG3a(1),time = t,type = x)}
#av45P1 = function(t){geePlot1('av45',t)}
#hocP1 = function(t){geePlot1('hoc',t)}
#fdgP1 = function(t){geePlot1('fdg',t)}


geePlot1 = function(x,t){fun1(geeG3a(1),time = t,type = x)}
av45P1 = function(t){geePlot1('av45',t)}
hocP1 = function(t){geePlot1('hoc',t)}
fdgP1 = function(t){geePlot1('fdg',t)}

geePlot2 = function(x,t){fun1(geeG3a(2),time = t,type = x)}
av45P2 = function(t){geePlot2('av45',t)}
hocP2 = function(t){geePlot2('hoc',t)}
fdgP2 = function(t){geePlot2('fdg',t)}

geePlot3 = function(x,t){fun1(geeG3a(3),time = t,type = x)}
av45P3 = function(t){geePlot3('av45',t)}
hocP3 = function(t){geePlot3('hoc',t)}
fdgP3 = function(t){geePlot3('fdg',t)}

r = ggplot(data.frame(x = c(-30,30)), aes(x))

r1 = r + theme_classic()+ coord_cartesian(ylim = c(0,1)) +theme(legend.position = 'bottom')+
  stat_function(aes(linetype = 'AV45'),size = 1,fun = av45P1) +
  stat_function(aes(linetype = 'HOC'),size = 1,fun = hocP1) +
  stat_function(aes(linetype = 'FDG'),size = 1,fun = fdgP1)+
  theme(legend.position='bottom',plot.title = element_text(size = 20),legend.title=element_blank(),
        legend.text = element_text(size = 20),legend.key.size = unit(1.45, "cm"))+
  theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
  theme(axis.text.x = element_text(size = 15,face = 'bold',color = 'black'))+
  theme(axis.text.y = element_text(size = 15,face = 'bold',color = 'black'))+
  ylab("Biomarker Severity")+xlab("Time")+scale_linetype_manual(values=c('solid','dotted','longdash','solid'))

r2 = r + theme_classic()+ coord_cartesian(ylim = c(0,1))+
  stat_function(aes(linetype = 'AV45'),size = 1,fun = av45P2) +
  stat_function(aes(linetype = 'HOC'),size = 1,fun = hocP2) +
  stat_function(aes(linetype = 'FDG'),size = 1,fun = fdgP2)+
  theme(legend.position='bottom',plot.title = element_text(size = 20),legend.title=element_blank(),
        legend.text = element_text(size = 20),legend.key.size = unit(1.45, "cm"))+
  theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
  theme(axis.text.x = element_text(size = 15,face = 'bold',color = 'black'))+
  theme(axis.text.y = element_text(size = 15,face = 'bold',color = 'black'))+
  ylab("Biomarker Severity")+xlab("Time")+scale_linetype_manual(values=c('solid','dotted','longdash','solid'))


r3 = r + theme_classic()+ coord_cartesian(ylim = c(0,1))+
  stat_function(aes(linetype = 'AV45'),size = 1,fun = av45P3) +
  stat_function(aes(linetype = 'HOC'),size = 1,fun = hocP3) +
  stat_function(aes(linetype = 'FDG'),size = 1,fun = fdgP3)+
  theme(legend.position='bottom',plot.title = element_text(size = 20),legend.title=element_blank(),
        legend.text = element_text(size = 20),legend.key.size = unit(1.45, "cm"))+
  theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
  theme(axis.text.x = element_text(size = 15,face = 'bold',color = 'black'))+
  theme(axis.text.y = element_text(size = 15,face = 'bold',color = 'black'))+
  ylab("Biomarker Severity")+xlab("Time")+scale_linetype_manual(values=c('solid','dotted','longdash','solid'))


mylegendR = g_legend(r1)
png('~/Documents/Dissertation/Dissertation_2017/figures/gee3.png',width = 1600, height = 800)

hTTS3 = grid.arrange(arrangeGrob(r1 + theme(legend.position="none")+ggtitle('Group 1a'),
                                 r3 + theme(legend.position="none")+ylab("")+ggtitle('Group 1b'),
                                 r2 + theme(legend.position="none")+ylab("")+ggtitle('Group 2'),
                                 nrow=1),
                     mylegendR,nrow = 2,heights=c(10, 1))


dev.off()

geePlot1 = function(x,t){fun1(geeG2a(1),time = t,type = x)}
av45P1 = function(t){geePlot1('av45',t)}
hocP1 = function(t){geePlot1('hoc',t)}
fdgP1 = function(t){geePlot1('fdg',t)}

geePlot2 = function(x,t){fun1(geeG2a(2),time = t,type = x)}
av45P2 = function(t){geePlot2('av45',t)}
hocP2 = function(t){geePlot2('hoc',t)}
fdgP2 = function(t){geePlot2('fdg',t)}

r = ggplot(data.frame(x = c(-30,30)), aes(x))

r1 = r + theme_classic() + coord_cartesian(ylim = c(0,1)) +theme(legend.position = 'bottom')+
  stat_function(aes(linetype = 'AV45'),size = 1,fun = av45P1) +
  stat_function(aes(linetype = 'HOC'),size = 1,fun = hocP1) +
  stat_function(aes(linetype = 'FDG'),size = 1,fun = fdgP1)+
  theme(legend.position='bottom',plot.title = element_text(size = 20),legend.title=element_blank(),
        legend.text = element_text(size = 20),legend.key.size = unit(1.45, "cm"))+
  theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
  theme(axis.text.x = element_text(size = 15,face = 'bold',color = 'black'))+
  theme(axis.text.y = element_text(size = 15,face = 'bold',color = 'black'))+
  ylab("Biomarker Severity")+xlab("Time")+scale_linetype_manual(values=c('solid','dotted','longdash','solid'))

r2 = r + theme_classic() + coord_cartesian(ylim = c(0,1))+
  stat_function(aes(linetype = 'AV45'),size = 1,fun = av45P2) +
  stat_function(aes(linetype = 'HOC'),size = 1,fun = hocP2) +
  stat_function(aes(linetype = 'FDG'),size = 1,fun = fdgP2) +
  theme(legend.position='bottom',plot.title = element_text(size = 20),legend.title=element_blank(),
        legend.text = element_text(size = 20),legend.key.size = unit(1.45, "cm"))+
  theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
  theme(axis.text.x = element_text(size = 15,face = 'bold',color = 'black'))+
  theme(axis.text.y = element_text(size = 15,face = 'bold',color = 'black'))+
  ylab("Biomarker Severity")+xlab("Time")+scale_linetype_manual(values=c('solid','dotted','longdash','solid'))



mylegendR = g_legend(r1)
png('~/Documents/Dissertation/Dissertation_2017/figures/gee2.png',width = 1600, height = 800)
hTTS2 = grid.arrange(arrangeGrob(r2 + theme(legend.position="none")+ggtitle('Group 1'),
                                 r1 + theme(legend.position="none")+ylab("")+ggtitle('Group 2'),
                                 nrow=1),
                     mylegendR,nrow = 2,heights=c(10, 1))

dev.off()