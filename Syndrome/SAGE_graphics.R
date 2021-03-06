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
  geom_line(data = dd1,aes(argvals + gamma0, y, linetype = 'AV45',group = id)) + 
  geom_point(data = dd1,aes(argvals + gamma0, y,linetype = 'AV45'), size = 1, alpha = 1)+ 
  geom_line(data = dd2,aes(argvals + gamma0, y, linetype = 'FDG',group = id)) + 
  geom_point(data = dd2,aes(argvals + gamma0, y,linetype = 'FDG'), size = 1, alpha = 1)+
  geom_line(data = dd3,aes(argvals + gamma0, y, linetype = 'HOC',group = id)) + 
  geom_point(data = dd3,aes(argvals + gamma0, y,linetype = 'HOC'), size = 1, alpha = 1)+
  theme(legend.position='bottom',plot.title = element_text(size = 20),legend.title=element_blank(),
  legend.text = element_text(size = 20),legend.key.size = unit(1.45, "cm"))+
  theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
  theme(axis.text.x = element_text(size = 15,face = 'bold',color = 'black'))+
  theme(axis.text.y = element_text(size = 15,face = 'bold',color = 'black'))+
  ylab("Biomarker Severity")+xlab("Adjusted Time")+scale_linetype_manual(values=c('solid','dotted','longdash'))
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





g1=readRDS('~/Documents/Dissertation/R_Dissertation/Syndrome/g2T1.RDA')
g2=readRDS('~/Documents/Dissertation/R_Dissertation/Syndrome/g2T2.RDA')

h1=readRDS('~/Documents/Dissertation/R_Dissertation/Syndrome/g3T1.RDA')
h2=readRDS('~/Documents/Dissertation/R_Dissertation/Syndrome/g3T2.RDA')
h3=readRDS('~/Documents/Dissertation/R_Dissertation/Syndrome/g3T3.RDA')

gT1 = pF(g1)
gT2 = pF(g2)  

gTT1 = pF1(g1)
gTT2 = pF1(g2)  
mylegend6<-g_legend(gT1)
gTS2 =  grid.arrange(arrangeGrob(gT2 + theme(legend.position="none"),
                               gT1 + theme(legend.position="none")+ylab(""),
                               nrow=1),
                   mylegend6,nrow = 2,heights=c(10, 1))

png("~/Documents/Dissertation/SyndromeClassification/SAGE/figures/syn2b.png",width = 800, height = 800)
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
mylegend6<-g_legend(hT1)
png('~/Documents/Dissertation/SyndromeClassification/AnnalsOfStat/figures/syn3bscatt.png',width = 800, height = 800)
hTS3 =  grid.arrange(arrangeGrob(hT1 + theme(legend.position="none"),
                                 hT3 + theme(legend.position="none")+ylab(""),
                                 hT2 + theme(legend.position="none")+ylab(""),
                                 nrow=1),
                     mylegend6,nrow = 2,heights=c(10, 1))
dev.off()
png("~/Documents/Dissertation/SyndromeClassification/SAGE/figures/syn3b.png",width = 800, height = 800)

hTTS3 = grid.arrange(arrangeGrob(hTT1 + theme(legend.position="none")+ggtitle('Group 1a'),
                                 hTT3 + theme(legend.position="none")+ylab("")+ggtitle('Group 1b'),
                                 hTT2 + theme(legend.position="none")+ylab("")+ggtitle('Group 2'),
                                 nrow=1),
                     mylegend6,nrow = 2,heights=c(10, 1))

dev.off()



c2g1sim =ddply(data.frame(rbind(g1$fits$av45$subset,g1$fits$fdg$subset,g1$fits$hoc$subset)),'outcome',function(x){sum(x$smooth.resids^2)})
c2g2sim =ddply(data.frame(rbind(g2$fits$av45$subset,g2$fits$fdg$subset,g2$fits$hoc$subset)),'outcome',function(x){sum(x$smooth.resids^2)})

c2h1sim =ddply(data.frame(rbind(h1$fits$av45$subset,h1$fits$fdg$subset,h1$fits$hoc$subset)),'outcome',function(x){sum(x$smooth.resids^2)})
c2h2sim =ddply(data.frame(rbind(h2$fits$av45$subset,h2$fits$fdg$subset,h2$fits$hoc$subset)),'outcome',function(x){sum(x$smooth.resids^2)})
c2h3sim =ddply(data.frame(rbind(h3$fits$av45$subset,h3$fits$fdg$subset,h3$fits$hoc$subset)),'outcome',function(x){sum(x$smooth.resids^2)})
