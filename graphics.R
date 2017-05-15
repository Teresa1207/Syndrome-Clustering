################################################################################
######    Syndrome Classifiation  		
######	  By: Teresa Filshtein
######	  Date: March 14, 2016			                      
######	  Modified: 
######    EDITED BY: 
################################################################################
#setwd("~/Documents/Dissertation/R_Dissertation/Syndrome")
################################################################################

#### Syndrome Classification
#### Define 2 Syndromes; 2 onsets (Early/Late) - 4 groups 
#### S1: X < Y < Z
#### S2: Y < X < Z
#### X initiates Y initiates Z, Y initiates X initiates Z.
#### Assume time lags of 5 years between process inititation
#### Early Onset: 60, Late Onset 80

#### Syndrome 1
b0 = -1.28
b1 = 1.28/5


s60 = function(t){pnorm(b0+b1*(t-60)*I(t>=60))}
s65 = function(t){pnorm(b0+b1*(t-65)*I(t>=65))}
s66 = function(t){pnorm(b0+b1*(t-66)*I(t>=66))}
s70 = function(t){pnorm(b0+b1*(t-70)*I(t>=70))}
s71 = function(t){pnorm(b0+b1*(t-71)*I(t>=71))}

s75 = function(t){pnorm(b0+b1*(t-75)*I(t>=75))}
s77 = function(t){pnorm(b0+b1*(t-77)*I(t>=77))}

s80 = function(t){pnorm(b0+b1*(t-80)*I(t>=80))}
s82 = function(t){pnorm(b0+b1*(t-82)*I(t>=82))}

setwd("~/Documents/Dissertation/SyndromeClassification/AnnalsOfStat/figures")

df1 = data.frame(Age = seq(60,90,by = .5))
size1 = 1.2
graphics1 = theme_classic() + 
  theme(axis.text.x = element_text(size = 18,face = 'bold'),
        axis.text.y = element_text(size = 18,face = 'bold'),
        axis.title.y = element_text(size = 22, face = 'bold'),
        axis.title.x = element_text(size = 16, face = 'bold')
  )
g1 = ggplot(df1, aes(x = Age)) + 
  stat_function(fun = s75,aes(linetype = 'FDG'),size = size1) + xlim(c(73,95)) + graphics1 + ylab('')
g2 = ggplot(df1, aes(x = Age)) + 
  stat_function(fun = s75,aes(linetype = 'FDG'),size = size1) + xlim(c(73,95)) + graphics1 + ylab('Biomarker Severity')
g3 = ggplot(df1, aes(x = Age)) + 
  stat_function(fun = s82,aes(linetype = 'FDG'),size = size1) + xlim(c(73,100)) + graphics1 + ylab('')+coord_cartesian(xlim = c(73,95))

g11 = g1 +
 stat_function(fun = s70, aes(linetype = 'Tau'), size = size1) + 
  xlim(c(70,95)) +
  scale_linetype_manual(values = c('solid','dotdash'), guide_legend(title = '')) +
  theme(legend.position = 'bottom')+
  theme(legend.key.size = unit(1.5, "cm"),
          legend.text = element_text(size =16))


g22 = g2  + stat_function(fun = s80, aes(linetype = 'Tau'),size = size1)+ xlim(c(70,95)) +
  scale_linetype_manual(values = c('solid','dotdash'), guide_legend(title = ''))+theme(legend.position = 'bottom')


g33 = g3  + stat_function(fun = s77, aes(linetype = 'Tau'),size = size1)+ xlim(c(70,95)) +
  scale_linetype_manual(values = c('solid','dotdash'), guide_legend(title = ''))+theme(legend.position = 'bottom')


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

leg = g_legend(g11)


png("~/Documents/Dissertation/SyndromeClassification/SAGE/figures/plotSage.png",height = 800, width = 500)
grid.arrange(arrangeGrob(g1 + theme(legend.position="none",plot.title = element_text(size = 20))+ggtitle('FDG'),
                         g11 + theme(legend.position="none",plot.title = element_text(size = 20))+ggtitle('FDG, Tau'),
                         g2 + theme(legend.position = 'none'),
                         g22 + theme(legend.position="none")+ylab(''),
                         g3 + theme(legend.position="none"),
                         g33 + theme(legend.position = 'none'),
                         nrow=3),leg,nrow = 2,heights=c(10,1))
dev.off()


library(ggplot2)
library(Rmisc)
library(gridExtra)
library(grid)
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

#### Syndrome 1
b0 = -1
b1 = .2


f1 = function(t){pnorm(b0+b1*(t-60))}
f2 = function(t){pnorm(b0+b1*(t-65))}
f3 = function(t){pnorm(b0+b1*(t-68))}
f4 = function(t){pnorm(b0+b1*(t-80))}
f5 = function(t){pnorm(b0+b1*(t-85))}
f6 = function(t){pnorm(b0+b1*(t-88))}

l1 = function(t){b0+b1*(t-60)}
l2 = function(t){b0+b1*(t-65)}
###################################
###### Graphic 1a and 1b ##########
###################################
### Classic Amyloid Syndrome
### 2 people with the same Syndrome but different times 
### (early/late onset)par(mfrow = c(2,2))
myColors <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
df = data.frame(AgeE = seq(55,80,by = .1),AgeL =seq(75,100,by = .1))
base = ggplot()
g1a = base+
  stat_function(data = df,fun = f1,aes(x = AgeE,color = 'CSF ABeta'),size =2.5)+
  stat_function(data = df,fun = f2,aes(x = AgeE,color = 'CSF Tau'),size =2.5)+
  stat_function(data = df,fun = f3,aes(x = AgeE,color = 'MRI+FDGPET'),size =2.5)+
  theme(legend.position='bottom',plot.title = element_text(size = 20),legend.title=element_blank(),
        legend.text = element_text(size = 20),legend.key.size = unit(1.45, "cm"))+
  theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
  theme(axis.text.x = element_text(size = 15,face = 'bold',color = 'black'))+
  theme(axis.text.y = element_text(size = 15,face = 'bold',color = 'black'))+
  ylab("Biomarker Severity")+xlab("Age")+ylim(0.001,.999)+scale_color_manual(values=myColors)+
  ggtitle("Amyloid Syndrome")

mylegend1a<-g_legend(g1a)


png('graphic1a.png',width = 800)
g1aP<- grid.arrange(arrangeGrob(g1a + theme(legend.position="none"),
                                nrow=1),
                    mylegend1a,nrow = 2,heights=c(10, 1))
dev.off()

g1b = base+
  stat_function(data = df,fun = f1,aes(x = AgeE,color = 'CSF ABeta'),size =2.5)+
  stat_function(data = df,fun = f2,aes(x = AgeE,color = 'CSF Tau'),size =2.5)+
  stat_function(data = df,fun = f3,aes(x = AgeE,color = 'MRI+FDGPET'),size =2.5)+
  stat_function(data = df,fun = f4,aes(x = AgeL,color = 'CSF ABeta'),size =2.5)+
  stat_function(data = df,fun = f5,aes(x = AgeL,color = 'CSF Tau'),size =2.5)+
  stat_function(data = df,fun = f6,aes(x = AgeL,color = 'MRI+FDGPET'),size =2.5)+
  theme(legend.position='bottom',plot.title = element_text(size = 20),legend.title=element_blank(),
        legend.text = element_text(size = 20),legend.key.size = unit(1.45, "cm"))+
  theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
  theme(axis.text.x = element_text(size = 15,face = 'bold',color = 'black'))+
  theme(axis.text.y = element_text(size = 15,face = 'bold',color = 'black'))+
  ylab("Biomarker Severity")+xlab("Age")+ylim(0.001,.999)+scale_color_manual(values=myColors)+
  ggtitle("Amyloid Syndrome: Early/Late Onset")

mylegend1b<-g_legend(g1b)


png('graphic1b.png',width = 800)
g1bP<- grid.arrange(arrangeGrob(g1b + theme(legend.position="none"),
                                nrow=1),
                    mylegend1b,nrow = 2,heights=c(10, 1))
dev.off()


###################################
###### Graphic 2 and 3 ############
###################################
### 2 people with similar Amyloid, but other symptoms are out of order
### (early/late onset)par(mfrow = c(2,2))

df = data.frame(AgeE = seq(55,80,by = .1),AgeL =seq(75,100,by = .1))
base = ggplot()
g2a = base+
  stat_function(data = df,fun = f2,aes(x = AgeE,color = 'CSF ABeta'),size =2.5)+
  theme(legend.position='bottom',plot.title = element_text(size = 20),legend.title=element_blank(),
        legend.text = element_text(size = 20),legend.key.size = unit(1.45, "cm"))+
  theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
  theme(axis.text.x = element_text(size = 15,face = 'bold',color = 'black'))+
  theme(axis.text.y = element_text(size = 15,face = 'bold',color = 'black'))+
  ylab("Biomarker Severity")+xlab("Age")+ylim(0.001,.999)+scale_color_manual(values=myColors)
g2b = base+
  stat_function(data = df,fun = f2,aes(x = AgeE,color = 'CSF ABeta'),size =2.5)+
  theme(legend.position='bottom',plot.title = element_text(size = 20),legend.title=element_blank(),
        legend.text = element_text(size = 20),legend.key.size = unit(1.45, "cm"))+
  theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
  theme(axis.text.x = element_text(size = 15,face = 'bold',color = 'black'))+
  theme(axis.text.y = element_text(size = 15,face = 'bold',color = 'black'))+
  ylab("Biomarker Severity")+xlab("Age")+ylim(0.001,.999)+scale_color_manual(values=myColors)
mylegend1<-g_legend(g2a)

png('graphic2.png',width = 800)
g2ab<- grid.arrange(arrangeGrob(g2a + theme(legend.position="none"),
                                g2b + theme(legend.position="none"),
                                nrow=2),
                    mylegend1,nrow = 2,heights=c(10, 1))

dev.off()
g2c = base+
  stat_function(data = df,fun = f2,aes(x = AgeE,color = 'CSF ABeta'),size =2.5)+
  stat_function(data = df,fun = f1,aes(x = AgeE,color = 'MRI+FDGPET'),size =2.5)+
  stat_function(data = df,fun = f3,aes(x = AgeE,color = 'CSF Tau'),size =2.5)+
  theme(legend.position='bottom',plot.title = element_text(size = 20),legend.title=element_blank(),
        legend.text = element_text(size = 20),legend.key.size = unit(1.45, "cm"))+
  theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
  theme(axis.text.x = element_text(size = 15,face = 'bold',color = 'black'))+
  theme(axis.text.y = element_text(size = 15,face = 'bold',color = 'black'))+scale_color_manual(values=myColors)+
  ylab("Biomarker Severity")+xlab("Age")+ylim(0.001,.999)
g2d = base+
  stat_function(data = df,fun = f2,aes(x = AgeE,color = 'CSF ABeta'),size =2.5)+
  stat_function(data = df,fun = f3,aes(x = AgeE,color = 'MRI+FDGPET'),size =2.5)+
  stat_function(data = df,fun = f1,aes(x = AgeE,color = 'CSF Tau'),size =2.5)+
  theme(legend.position='bottom',plot.title = element_text(size = 20),legend.title=element_blank(),
        legend.text = element_text(size = 20),legend.key.size = unit(1.45, "cm"))+
  theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
  theme(axis.text.x = element_text(size = 15,face = 'bold',color = 'black'))+
  theme(axis.text.y = element_text(size = 15,face = 'bold',color = 'black'))+scale_color_manual(values=myColors)+
  ylab("Biomarker Severity")+xlab("Age")+ylim(0.001,.999)
mylegend2<-g_legend(g2c)

png('graphic3.png',width = 800)
g2cd<- grid.arrange(arrangeGrob(g2c + theme(legend.position="none"),
                                g2d + theme(legend.position="none"),
                                nrow=2),
                    mylegend2,nrow = 2,heights=c(10, 1))
dev.off()


###################################
###### Graphic 4 ##################
###################################
# Show the shaded region time invariant part
# build plot object for rendering 

## 4a, 2 curves
df = data.frame(AgeE = seq(55,80,by = .1),AgeL =seq(75,100,by = .1))
base = ggplot()
g4a = base+theme_classic()+
  stat_function(data = df,fun = f1,aes(x = AgeE,linetype = 'Marker 1'),size =1)+
  stat_function(data = df,fun = f2,aes(x = AgeE,linetype = 'Marker 2'),size =1)+
  theme(legend.position='bottom',plot.title = element_text(size = 10),legend.title=element_blank(),
        legend.text = element_text(size = 10),legend.key.size = unit(1, "cm"))+
  theme(axis.title.x = element_text(size = 10),axis.title.y = element_text(size = 10))+
  theme(axis.text.x = element_text(size = 10,face = 'bold',color = 'black'))+
  theme(axis.text.y = element_text(size = 10,face = 'bold',color = 'black'))+
  ylab("Biomarker Severity")+xlab("Age")+ylim(0.001,.999)+scale_color_manual(values=myColors)+
  theme(legend.position = 'bottom')+
  theme(legend.key.size = unit(1.5, "cm"),
        legend.text = element_text(size =16))+
  theme(axis.text.x= element_text(size = 20), 
        axis.text.y= element_text(size = 20), axis.title.x = element_text(size = 24), 
        axis.title.y = element_text(size = 24)) 
l4a = base+theme_classic()+
  stat_function(data = df,fun = l1,aes(x = AgeE,linetype = 'Marker 1'),size =1)+
  stat_function(data = df,fun = l2,aes(x = AgeE,linetype = 'Marker 2'),size =1)+
  theme(legend.position='bottom',plot.title = element_text(size = 10),legend.title=element_blank(),
        legend.text = element_text(size = 10),legend.key.size = unit(1.45, "cm"))+
  theme(axis.title.x = element_text(size = 10),axis.title.y = element_text(size = 10))+
  theme(axis.text.x = element_text(size = 10,face = 'bold',color = 'black'))+
  theme(axis.text.y = element_text(size = 10,face = 'bold',color = 'black'))+
  ylab("Linear Representation")+xlab("Age")+scale_color_manual(values=myColors)+
  geom_vline(xintercept = c(65,70),linetype = 'longdash',size = 1)
mylegendg4ls = g_legend(g4a)
ll4 <- ggplot_build(l4a)
lf2 <- data.frame(x = ll4$data[[1]]$x[which(ll4$data[[1]]$x>=65&ll4$data[[1]]$x<=70)],
                  ymin = ll4$data[[1]]$y[which(ll4$data[[1]]$x>=65&ll4$data[[1]]$x<=70)],
                  ymax = ll4$data[[2]]$y[which(ll4$data[[1]]$x>=65&ll4$data[[1]]$x<=70)]) 

# use the loess data to add the 'ribbon' to plot 
l4S = l4a + geom_ribbon(data = lf2, aes(x = x, ymin = ymin, ymax = ymax),
                        fill = "grey1", alpha = 0.4)+
  theme(axis.text.x= element_text(size = 20), 
        axis.text.y= element_text(size = 20), axis.title.x = element_text(size = 24), 
        axis.title.y = element_text(size = 24)) 

png('/Users/Teresa/Documents/Dissertation/SyndromeClassification/SAGE/figures/g4ls.png',width = 800)
g4ls<-grid.arrange(arrangeGrob(g4a + theme(legend.position="none"),
                               l4S + theme(legend.position="none"),
                               nrow=1),
                   mylegendg4ls,nrow = 2,heights=c(10, 1))
dev.off()
