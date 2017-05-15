###################################
###### Graphics Grace #############
###################################

library(ggplot2)
library(Rmisc)
library(gridExtra)
library(grid)
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
myColors <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Show ADNI clusters
#pull from ClusterSetup_ADNI.R
pF = function(f){
dd1 = f$fits$av45$subset
dd2 = f$fits$fdg$subset
dd3 = f$fits$hoc$subset
d1 = ggplot() 
d1+ theme_classic()+
  geom_line(data = dd1,aes(argvals + gamma0, y, color = group,group = id)) + 
  geom_point(data = dd1,aes(argvals + gamma0, y,color = group), size = 1, alpha = 1)+ 
  geom_line(data = dd2,aes(argvals + gamma0, y, color = group,group = id)) + 
  geom_point(data = dd2,aes(argvals + gamma0, y,color = group), size = 1, alpha = 1)+
  geom_line(data = dd3,aes(argvals + gamma0, y, color = group,group = id)) + 
  geom_point(data = dd3,aes(argvals + gamma0, y,color = group), size = 1, alpha = 1)+
  geom_line(data = dd1,aes(argvals +gamma0, ghat, group = NULL,linetype = 'AV45'), size = 1.5)+
  geom_line(data = dd2,aes(argvals +gamma0, ghat, group = NULL,linetype = 'FDG'), size = 1.5)+
  geom_line(data = dd3,aes(argvals +gamma0, ghat, group = NULL,linetype = 'HOC'), size = 1.5)+
  theme(legend.position='bottom',plot.title = element_text(size = 20),legend.title=element_blank(),
  legend.text = element_text(size = 20),legend.key.size = unit(1.45, "cm"))+
  theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
  theme(axis.text.x = element_text(size = 15,face = 'bold',color = 'black'))+
  theme(axis.text.y = element_text(size = 15,face = 'bold',color = 'black'))+
  ylab("Biomarker Severity")+xlab("Adjusted Time")+scale_linetype_manual(values=c('solid','dotted','longdash'))+facet_grid(~group)
}


pF1 = function(f){
  dd1 = f$fits$av45$subset
  dd2 = f$fits$fdg$subset
  dd3 = f$fits$hoc$subset
  d1 = ggplot() 
  d1+ theme_classic()+
  geom_line(data = dd1,aes(argvals +gamma0, ghat, group = NULL,linetype = 'AV45'), size = 1.5)+
  geom_line(data = dd2,aes(argvals +gamma0, ghat, group = NULL, linetype = 'FDG'),size = 1.5)+
  geom_line(data = dd3,aes(argvals +gamma0, ghat, group = NULL, linetype = 'HOC'),size = 1.5)+
    theme(legend.position='bottom',plot.title = element_text(size = 20),legend.title=element_blank(),
          legend.text = element_text(size = 20),legend.key.size = unit(1.45, "cm"))+
    theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
    theme(axis.text.x = element_text(size = 15,face = 'bold',color = 'black'))+
    theme(axis.text.y = element_text(size = 15,face = 'bold',color = 'black'))+
    ylab("Biomarker Severity")+xlab("Adjusted Time")+scale_linetype_manual(values=c('solid','dotted','longdash'))

}

pFcomp = function(f){
  dd1 = f$fits$av45$subset
  dd2 = f$fits$fdg$subset
  dd3 = f$fits$hoc$subset
  d1 = ggplot() 
  d1+ theme_classic()+
    geom_line(data = dd1,aes(argvals + gamma0, ghat, group = NULL,linetype = 'AV45'), size = 1.5)
    geom_line(data = dd2,aes(argvals +gamma0, ghat, group = NULL, linetype = 'FDG'),size = 1.5)+
    geom_line(data = dd3,aes(argvals +gamma0, ghat, group = NULL, linetype = 'HOC'),size = 1.5)+
      geom_point(data = dd1,aes(argvals + gamma0, ghat, color = group), size = 1.5)+
    theme(legend.position='bottom',plot.title = element_text(size = 20),legend.title=element_blank(),
          legend.text = element_text(size = 20),legend.key.size = unit(1.45, "cm"))+
    theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
    theme(axis.text.x = element_text(size = 15,face = 'bold',color = 'black'))+
    theme(axis.text.y = element_text(size = 15,face = 'bold',color = 'black'))+
    ylab("Biomarker Severity")+xlab("Adjusted Time")+scale_linetype_manual(values=c('solid','dotted','longdash'))
  
}


pFcomp(gcomb)


pF1a = function(f){
  dd1 = f$fits$av45$subset
  dd2 = f$fits$fdg$subset
  dd3 = f$fits$hoc$subset
  dd4 = f$fits$adas$subset
  d1 = ggplot() 
  d1+ theme_classic()+
    geom_line(data = dd1,aes(argvals +gamma0, ghat, group = NULL,linetype = 'AV45'), size = 1.5)+
    geom_line(data = dd2,aes(argvals +gamma0, ghat, group = NULL, linetype = 'FDG'),size = 1.5)+
    geom_line(data = dd3,aes(argvals +gamma0, ghat, group = NULL, linetype = 'HOC'),size = 1.5)+
    geom_line(data = dd4,aes(argvals +gamma0, ghat, group = NULL, linetype = 'ADAS'),color = 'red',size = 2)+
    theme(legend.position='bottom',plot.title = element_text(size = 20),legend.title=element_blank(),
          legend.text = element_text(size = 20),legend.key.size = unit(1.45, "cm"))+
    theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
    theme(axis.text.x = element_text(size = 15,face = 'bold',color = 'black'))+
    theme(axis.text.y = element_text(size = 15,face = 'bold',color = 'black'))+
    ylab("Biomarker Severity")+xlab("Adjusted Time")+scale_linetype_manual(values=c('solid','dotted','longdash','solid'))
  
}



#png('graphic6a.png',width = 800)
#g6A<- grid.arrange(arrangeGrob(g6a + theme(legend.position="none"),
#                               g6b + theme(legend.position="none")+ylab(""),
#                                g6c + theme(legend.position="none")+ylab(""),
#                                nrow=1),
#                    mylegend6,nrow = 2,heights=c(10, 1))
#dev.off()



#png('graphic6b.png',width = 800)
#g6B<- grid.arrange(arrangeGrob(g6e + theme(legend.position="none"),
#                               g6f + theme(legend.position="none")+ylab(""),
#                               g6g + theme(legend.position="none")+ylab(""),
#                               nrow=1),
#                   mylegend6,nrow = 2,heights=c(10, 1))
#dev.off()
gcomb=readRDS('~/Documents/Dissertation/R_Dissertation/Syndrome/gTSplit.RDA')

g1=readRDS('~/Documents/Dissertation/R_Dissertation/Syndrome/g2T1.RDA')
g2=readRDS('~/Documents/Dissertation/R_Dissertation/Syndrome/g2T2.RDA')

h1=readRDS('~/Documents/Dissertation/R_Dissertation/Syndrome/g3T1.RDA')
h2=readRDS('~/Documents/Dissertation/R_Dissertation/Syndrome/g3T2.RDA')
h3=readRDS('~/Documents/Dissertation/R_Dissertation/Syndrome/g3T3.RDA')

gCombplot = pF1(gcomb)

gT1 = pF(g1)
gT2 = pF(g2)  

gTT1 = pF1(g1)
gTT2 = pF1(g2)  
mylegend6<-g_legend(gT1)
gTS2 =  grid.arrange(arrangeGrob(gT2 + theme(legend.position="none"),
                               gT1 + theme(legend.position="none")+ylab(""),
                               nrow=1),
                   mylegend6,nrow = 2,heights=c(10, 1))

png("~/Documents/Dissertation/Dissertation_2017/figures/syn2b.png",width = 1600, height = 800)
gTTS2 = grid.arrange(arrangeGrob(gTT2 + theme(legend.position="none")+ggtitle('Group 1'),
                               gTT1 + theme(legend.position="none")+ylab("")+ggtitle('Group 2'),
                               nrow=1),
                   mylegend6,nrow = 2,heights=c(10, 1))
dev.off()



hT1 = pF(h1)
hT2 = pF(h2)  
hT3 = pF(h3)  
hTT1 = pF1(h1)
hTT2 = pF1(h2)  
hTT3 = pF1(h3) 

hTa1 = pF1a(h1)
hTa2 = pF1a(h2)  
hTa3 = pF1a(h3) 
mylegendA<-g_legend(hTa1)
hTA =  grid.arrange(arrangeGrob(hTa1 + theme(legend.position="none"),
                                 hTa3 + theme(legend.position="none")+ylab(""),
                                 hTa2 + theme(legend.position="none")+ylab(""),
                                 nrow=1),
                     mylegendA,nrow = 2,heights=c(10, 1))


mylegend6<-g_legend(hT1)
png('~/Documents/Dissertation/Dissertation_2017/figures/syn3bscatt.png',width = 1600, height = 800)
hTS3 =  grid.arrange(arrangeGrob(hT1 + theme(legend.position="none"),
                                 hT3 + theme(legend.position="none")+ylab(""),
                                 hT2 + theme(legend.position="none")+ylab(""),
                                 nrow=1),
                     mylegend6,nrow = 2,heights=c(10, 1))
dev.off()
png("~/Documents/Dissertation/Dissertation_2017/figures/syn3b.png",width = 1600, height = 800)

hTTS3 = grid.arrange(arrangeGrob(hTT1 + theme(legend.position="none")+ggtitle('Group 1a'),
                                 hTT3 + theme(legend.position="none")+ylab("")+ggtitle('Group 1b'),
                                 hTT2 + theme(legend.position="none")+ylab("")+ggtitle('Group 2'),
                                 nrow=1),
                     mylegend6,nrow = 2,heights=c(10, 1))

dev.off()
mylegendcomb<-g_legend(gCombplot)
png("~/Documents/Dissertation/Dissertation_2017/figures/total.png",width = 1600, height = 800)
hTTcomb= grid.arrange(arrangeGrob(gCombplot + theme(legend.position="none")+ggtitle('Total Group'),
                                 nrow=1),
                     mylegendcomb,nrow = 2,heights=c(10, 1))

dev.off()
### Look closer at ADNI groups. 
adniGroups2 = data.frame(id = c(as.character(g1$fits$adas$subset$id),
                                as.character(g2$fits$adas$subset$id)),
                         group =c(as.character(g1$fits$adas$subset$group),
                                  as.character(g2$fits$adas$subset$group)))
adniGroups2 = adniGroups2[!duplicated(adniGroups2$id),]
adniGroups3 = data.frame(id = c(as.character(h1$fits$adas$subset$id),
                                as.character(h2$fits$adas$subset$id),
                                as.character(h3$fits$adas$subset$id)),
                         group =c(as.character(h1$fits$adas$subset$group),
                                  as.character(h2$fits$adas$subset$group),
                                  as.character(h3$fits$adas$subset$group)))
adniGroups3 = adniGroups3[!duplicated(adniGroups3$id),]


source('~/Documents/Dissertation/R_Dissertation/Syndrome/ClusterSetup_ADNI_v3.R')
library(dplyr)

### look at final diagnosis
dFinal = tbl_df(dApril) %>% group_by(RID) %>%
  summarise(
    FinalDiag = max(dx3cat, na.rm = TRUE)
  )

adniG2 = merge(graceC3[!duplicated(graceC3$RID),],adniGroups2, by.x = "RID",by.y = 'id')
adniG2 = merge(adniG2,dApril[!duplicated(dApril$RID),c('RID','PTEDUCAT',"total130","abeta0","tau0")],by = "RID")
adniG2_dt = tbl_df(adniG2)
adniG2_dt$groupName = ifelse(adniG2_dt$group == 1,'Group 2','Group 1')
adniG2_dt = merge(adniG2_dt,dFinal, by = 'RID')

adniG3 = merge(graceC3[!duplicated(graceC3$RID),],adniGroups3, by.x = "RID",by.y = 'id')
adniG3 = merge(adniG3,dApril[!duplicated(dApril$RID),c('RID','PTEDUCAT',"total130","abeta0","tau0")],by = "RID")
adniG3_dt = tbl_df(adniG3)
adniG3_dt$groupName = NA
adniG3_dt$groupName = ifelse(adniG3_dt$group == 1,'Group 1a',ifelse(adniG3_dt$group == 2, 'Group 2','Group 1b'))
adniG3_dt = merge(adniG3_dt,dFinal, by = 'RID')


tab3 = adniG3_dt[!duplicated(adniG3_dt$RID),] %>% group_by(groupName) %>%
  summarise(
    Count = length(RID),
    Age  = paste(round(mean(age65+65),2)," (",round(sd(age65+65),2),")",sep = ""),
    Education  = paste(round(mean(PTEDUCAT,na.rm = TRUE),2)," (",round(sd(PTEDUCAT,na.rm = TRUE),2),")",sep = ""),
    ADAS  = paste(round(mean(total130,na.rm = TRUE),2)," (",round(sd(total130,na.rm = TRUE),2),")",sep = ""),
    '%Male' = round(100*length(RID[female == '0:Male'])/length(RID),2),
    '%APOE4' =  round(100*length(RID[apoe4 == 1])/length(RID),2),
    '%MCI' = round(100*length(RID[FinalDiag != '1:Normal'])/length(RID),2),
    abeta  = paste(round(mean(abeta0,na.rm = TRUE),2)," (",round(sd(abeta0,na.rm = TRUE),2),")",sep = ""),
    tau = paste(round(mean(tau0,na.rm = TRUE),2)," (",round(sd(tau0,na.rm = TRUE),2),")",sep = "")
  )

tab2 = adniG2_dt[!duplicated(adniG2_dt$RID),]%>% group_by(groupName) %>%
  summarise(
    Count = length(RID),
    Age  = paste(round(mean(age65+65),2)," (",round(sd(age65+65),2),")",sep = ""),
    Education  = paste(round(mean(PTEDUCAT,na.rm = TRUE),2)," (",round(sd(PTEDUCAT,na.rm = TRUE),2),")",sep = ""),
    ADAS  = paste(round(mean(total130,na.rm = TRUE),2)," (",round(sd(total130,na.rm = TRUE),2),")",sep = ""),
    '%Male' = round(100*length(RID[female == '0:Male'])/length(RID),2),
    '%APOE4' =  round(100*length(RID[apoe4 == 1])/length(RID),2),
    '%MCI' = round(100*length(RID[FinalDiag != '1:Normal'])/length(RID),2),
    abeta  = paste(round(mean(abeta0,na.rm = TRUE),2)," (",round(sd(abeta0,na.rm = TRUE),2),")",sep = ""),
    tau = paste(round(mean(tau0,na.rm = TRUE),2)," (",round(sd(tau0,na.rm = TRUE),2),")",sep = "")
  )
    

adniTotal = adniG3_dt[!duplicated(adniG3_dt$RID),]
adniTotal$Diagnosis = NA
adniTotal$Diagnosis = ifelse(adniTotal$FinalDiag=='1:Normal','1:Normal','2:MCI')

tab1 = adniTotal %>% group_by(Diagnosis) %>%
  summarise(
    Count = length(RID),
    Age  = paste(round(mean(age65+65),2)," (",round(sd(age65+65),2),")",sep = ""),
    Education  = paste(round(mean(PTEDUCAT,na.rm = TRUE),2)," (",round(sd(PTEDUCAT,na.rm = TRUE),2),")",sep = ""),
    ADAS  = paste(round(mean(total130,na.rm = TRUE),2)," (",round(sd(total130,na.rm = TRUE),2),")",sep = ""),
    '%Male' = round(100*length(RID[female == '0:Male'])/length(RID),2),
    '%APOE4' =  round(100*length(RID[apoe4 == 1])/length(RID),2),
    #'%MCI' = length(RID[base3dx == '2:MCI'])/length(RID),
    abeta  = paste(round(mean(abeta0,na.rm = TRUE),2)," (",round(sd(abeta0,na.rm = TRUE),2),")",sep = ""),
    tau = paste(round(mean(tau0,na.rm = TRUE),2)," (",round(sd(tau0,na.rm = TRUE),2),")",sep = "")
  )

library(xtable)
xtable(t(tab1),caption = "Desriptive Summary of ADNI sample")
xtable(t(tab2),caption = "Two Groups: Desriptive Summary of the subgroups found")
xtable(t(tab3),caption = "Three Groups: Desriptive Summary of the subgroups found")


### Test for differences in groups
tab3a = adniG3_dt[!duplicated(adniG3_dt$RID),] %>% group_by(groupName) %>%
  summarise(
    Count = length(RID),
    Age  = paste(round(mean(age65+65),2)," (",round(sd(age65+65),2),")",sep = ""),
    Education  = paste(round(mean(PTEDUCAT,na.rm = TRUE),2)," (",round(sd(PTEDUCAT,na.rm = TRUE),2),")",sep = ""),
    ADAS  = paste(round(mean(total130,na.rm = TRUE),2)," (",round(sd(total130,na.rm = TRUE),2),")",sep = ""),
    '#Male' = length(RID[female == '0:Male']),
    '#APOE4' =  length(RID[apoe4 == 1]),
    '#MCI' = length(RID[FinalDiag != '1:Normal']),
    abeta  = paste(round(mean(abeta0,na.rm = TRUE),2)," (",round(sd(abeta0,na.rm = TRUE),2),")",sep = ""),
    tau = paste(round(mean(tau0,na.rm = TRUE),2)," (",round(sd(tau0,na.rm = TRUE),2),")",sep = "")
  )

tab2a = adniG2_dt[!duplicated(adniG2_dt$RID),]%>% group_by(groupName) %>%
  summarise(
    Count = length(RID),
    Age  = paste(round(mean(age65+65),2)," (",round(sd(age65+65),2),")",sep = ""),
    Education  = paste(round(mean(PTEDUCAT,na.rm = TRUE),2)," (",round(sd(PTEDUCAT,na.rm = TRUE),2),")",sep = ""),
    ADAS  = paste(round(mean(total130,na.rm = TRUE),2)," (",round(sd(total130,na.rm = TRUE),2),")",sep = ""),
    '#Male' = length(RID[female == '0:Male']),
    '#APOE4' =  length(RID[apoe4 == 1]),
    '#MCI' = length(RID[FinalDiag != '1:Normal']),
    abeta  = paste(round(mean(abeta0,na.rm = TRUE),2)," (",round(sd(abeta0,na.rm = TRUE),2),")",sep = ""),
    tau = paste(round(mean(tau0,na.rm = TRUE),2)," (",round(sd(tau0,na.rm = TRUE),2),")",sep = "")
  )

tab3a = adniG3_dt[!duplicated(adniG3_dt$RID),]%>% group_by(groupName) %>%
  summarise(
    Count = length(RID),
    Age  = paste(round(mean(age65+65),2)," (",round(sd(age65+65),2),")",sep = ""),
    Education  = paste(round(mean(PTEDUCAT,na.rm = TRUE),2)," (",round(sd(PTEDUCAT,na.rm = TRUE),2),")",sep = ""),
    ADAS  = paste(round(mean(total130,na.rm = TRUE),2)," (",round(sd(total130,na.rm = TRUE),2),")",sep = ""),
    '#Male' = length(RID[female == '0:Male']),
    '#APOE4' =  length(RID[apoe4 == 1]),
    '#MCI' = length(RID[FinalDiag != '1:Normal']),
    abeta  = paste(round(mean(abeta0,na.rm = TRUE),2)," (",round(sd(abeta0,na.rm = TRUE),2),")",sep = ""),
    tau = paste(round(mean(tau0,na.rm = TRUE),2)," (",round(sd(tau0,na.rm = TRUE),2),")",sep = "")
  )



aov.out1 = aov(age65~groupName,adniG3_dt[!duplicated(adniG3_dt$RID),])
a1 = TukeyHSD(aov.out1)

aov.out2 = aov(PTEDUCAT~groupName,adniG3_dt[!duplicated(adniG3_dt$RID),])
a2 = TukeyHSD(aov.out2)

aov.out3 = aov(total130~groupName,adniG3_dt[!duplicated(adniG3_dt$RID),])
a3 = TukeyHSD(aov.out3)
a4 = pairwise.prop.test(tab3a$`#Male`,tab3a$Count)
a5 = pairwise.prop.test(tab3a$`#APOE4`,tab3a$Count)
a6 = pairwise.prop.test(tab3a$`#MCI`,tab3a$Count)

aov.out7 = aov(abeta0~groupName,adniG3_dt[!duplicated(adniG3_dt$RID),])
a7 = TukeyHSD(aov.out7)

aov.out8 = aov(tau0~groupName,adniG3_dt[!duplicated(adniG3_dt$RID),])
a8 = TukeyHSD(aov.out8)

aovP = function(x){
  df = data.frame(x$groupName)
  df$p.adj
}

propP = function(x){
  c(x$p.value[,1],x$p.value[2,2])
}
lev3 = c('Age','Education','ADAS','%Male','%APOE4','%MCI','abeta','tau')


p_groups3= data.frame(
  factor = unlist(lapply(lev3, function(x){rep(x,3)})),
  comparison = rep(row.names(a1$groupName),8),
  pVal = c(aovP(a1),aovP(a2),aovP(a3),propP(a4),propP(a5),propP(a6),aovP(a7),aovP(a8)),
  pVal.round = round(c(aovP(a1),aovP(a2),aovP(a3),propP(a4),propP(a5),propP(a6),aovP(a7),aovP(a8)),4))
p_groups3$p.adj = round(p.adjust(p_groups3$pVal,method = 'holm'),4)
p_groups3$sig = ifelse(p_groups3$p.adj<0.05,1,0)

dunc = lapply(list(aov.out1,aov.out2,aov.out3,aov.out7,aov.out8),
       function(x){duncan.test(x,"groupName",alpha = .05)})







t1 = t.test(age65~groupName,data = adniG2_dt[!duplicated(adniG2_dt$RID),])
t2 = t.test(PTEDUCAT~groupName,data = adniG2_dt[!duplicated(adniG2_dt$RID),])
t3 = t.test(total130~groupName,data = adniG2_dt[!duplicated(adniG2_dt$RID),]) #yes, but slight p = 0.04
t4 = prop.test(tab2a$`#Male`,tab2a$Count)
t5 = prop.test(tab2a$`#APOE4`,tab2a$Count)
t6 = prop.test(tab2a$`#MCI`,tab2a$Count)
t7 = t.test(abeta0~groupName,data = adniG2_dt[!duplicated(adniG2_dt$RID),]) 
t8 = t.test(tau0~groupName,data = adniG2_dt[!duplicated(adniG2_dt$RID),])

p_groups2 = data.frame(
  factor = colnames(tab2)[-c(1:2)],
  pVal = round(unlist(lapply(list(t1,t2,t3,t4,t5,t6,t7,t8),function(x){
  x$p.value})),4), 
  adjPval = round(p.adjust(unlist(lapply(list(t1,t2,t3,t4,t5,t6,t7,t8),function(x){x$p.value})),method = 'bonferroni'),4)
)

tab2c = data.frame(t(tab2),'Adj. p-val' = c("","",p_groups2$adjPval))
colnames(tab2c) = c(tab2c[1,1:2],'Adj. p-val')
tab2c = tab2c[-1,]
xtable(tab2c)

require(graphics)





c2g1sim =ddply(data.frame(rbind(g1$fits$av45$subset,g1$fits$fdg$subset,g1$fits$hoc$subset)),'outcome',function(x){sum(x$smooth.resids^2)})
c2g2sim =ddply(data.frame(rbind(g2$fits$av45$subset,g2$fits$fdg$subset,g2$fits$hoc$subset)),'outcome',function(x){sum(x$smooth.resids^2)})

c2h1sim =ddply(data.frame(rbind(h1$fits$av45$subset,h1$fits$fdg$subset,h1$fits$hoc$subset)),'outcome',function(x){sum(x$smooth.resids^2)})
c2h2sim =ddply(data.frame(rbind(h2$fits$av45$subset,h2$fits$fdg$subset,h2$fits$hoc$subset)),'outcome',function(x){sum(x$smooth.resids^2)})
c2h3sim =ddply(data.frame(rbind(h3$fits$av45$subset,h3$fits$fdg$subset,h3$fits$hoc$subset)),'outcome',function(x){sum(x$smooth.resids^2)})
