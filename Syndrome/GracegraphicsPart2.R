###################################
###### Graphics Grace #############
###################################

library(ggplot2)
library(Rmisc)
library(gridExtra)
library(grid)
source('~/Documents/Dissertation/R_Dissertation/Syndrome/ClusterSetup_ADNI_v3.R')

g1=readRDS('~/Documents/Dissertation/R_Dissertation/Syndrome/g2T1.RDA')
g2=readRDS('~/Documents/Dissertation/R_Dissertation/Syndrome/g2T2.RDA')

h1=readRDS('~/Documents/Dissertation/R_Dissertation/Syndrome/g3T1.RDA')
h2=readRDS('~/Documents/Dissertation/R_Dissertation/Syndrome/g3T2.RDA')
h3=readRDS('~/Documents/Dissertation/R_Dissertation/Syndrome/g3T3.RDA')
dxFunc = function(gx){
  dx = data.frame(RID = unique(as.character(gx$fits$adas$subset$id)))
  dx = merge(dx,mark.long[!duplicated(mark.long$RID),c('RID',"base3dx")],by = 'RID')
}

g1$dx = dxFunc(g1)
g2$dx = dxFunc(g2)

h1$dx = dxFunc(h1)
h2$dx = dxFunc(h2)
h3$dx = dxFunc(h3)

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
colP = cbbPalette[c(2,4,6)]

pFtype = function(f,type){
  dd1 = f[['fits']][[type]][['subset']]
  dd1 = merge(dd1,f$dx,by.x = 'id',by.y = 'RID',all.x  = TRUE)
  d1 = ggplot() 
  d1+ theme_classic()+
  geom_line(data = dd1,aes(argvals +gamma0, ghat, group = NULL), size = 2)+
    #geom_line(data = dd1,aes(argvals + gamma0, y,group = id,color = base3dx)) + 
    geom_point(data = dd1,aes(argvals + gamma0, y,color = base3dx), size = 3, alpha = 1)+
    theme(legend.position='bottom',plot.title = element_text(size = 20),legend.title=element_blank(),
          legend.text = element_text(size = 20),legend.key.size = unit(1.45, "cm"))+
    theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
    theme(axis.text.x = element_text(size = 15,face = 'bold',color = 'black'))+
    theme(axis.text.y = element_text(size = 15,face = 'bold',color = 'black'))+
    ylab("Biomarker Severity")+xlab("Adjusted Time")+
    scale_colour_manual(values=colP)
}


hT1a = pFtype(h1,'av45')
hT1b = pFtype(h1,'fdg')
hT1c = pFtype(h1,'hoc') 


hT2a = pFtype(h2,'av45')
hT2b = pFtype(h2,'fdg')
hT2c = pFtype(h2,'hoc') 


hT3a = pFtype(h3,'av45')
hT3b = pFtype(h3,'fdg')
hT3c = pFtype(h3,'hoc') 


mylegendA<-g_legend(hT1a)

hTA1 =  grid.arrange(arrangeGrob(hT1a + theme(legend.position="none")+ggtitle('AV45'),
                                 hT1b + theme(legend.position="none")+ylab("")+ggtitle('FDG'),
                                 hT1c + theme(legend.position="none")+ylab("")+ggtitle('HOC'),
                                 nrow=1),
                     mylegendA,nrow = 2,heights=c(10, 1),widths = c(10))

hTA2 =  grid.arrange(arrangeGrob(hT2a + theme(legend.position="none")+ggtitle('AV45'),
                                 hT2b + theme(legend.position="none")+ylab("")+ggtitle('FDG'),
                                 hT2c + theme(legend.position="none")+ylab("")+ggtitle('HOC'),
                                 nrow=1),
                     mylegendA,nrow = 2,heights=c(10, 1),widths = c(10))
hTA3 =  grid.arrange(arrangeGrob(hT3a + theme(legend.position="none")+ggtitle('AV45'),
                                 hT3b + theme(legend.position="none")+ylab("")+ggtitle('FDG'),
                                 hT3c + theme(legend.position="none")+ylab("")+ggtitle('HOC'),
                                 nrow=1),
                     mylegendA,nrow = 2,heights=c(10, 1),widths = c(10))
