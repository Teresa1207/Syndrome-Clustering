### Graphics for AAIC 2016 ###
setwd("/Users/Teresa/Documents/Davis2015-2016/AAIC/figures")
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

l4S

#4b cut off so not infinite
time = qnorm(seq(.10,.9,by=.05))
nt = length(time)
g4b = g4a+geom_vline(xintercept = c(59,77),linetype = 'dotted',size = 1)
g4c = g4a+geom_vline(xintercept = c(59,77),linetype = 'dotted',size = 1)+
  geom_line(arrow = arrow(ends = "both"),linetype = 'dashed',
            aes(x = seq(from = 60.5+(qnorm(.95)-b0)/b1,to = 60+(qnorm(.95)-b0)/b1+3.5,by = .1), y = .95))+
  geom_line(arrow = arrow(ends = "both"),linetype = 'dashed',
            aes(x = seq(from = 60.5+(qnorm(.85)-b0)/b1,to = 60+(qnorm(.85)-b0)/b1+4.5,by = .1), y = .85))+
  geom_line(arrow = arrow(ends = "both"),linetype = 'dashed',
            aes(x = seq(from = 60.5+(qnorm(.75)-b0)/b1,to = 60+(qnorm(.75)-b0)/b1+4.5,by = .1), y = .75))+
  geom_line(arrow = arrow(ends = "both"),linetype = 'dashed',
            aes(x = seq(from = 60.5+(qnorm(.65)-b0)/b1,to = 60+(qnorm(.65)-b0)/b1+4.5,by = .1), y = .65))+
  geom_line(arrow = arrow(ends = "both"),linetype = 'dashed',
            aes(x = seq(from = 60.5+(qnorm(.55)-b0)/b1,to = 60+(qnorm(.55)-b0)/b1+4.5,by = .1), y = .55))+
  geom_line(arrow = arrow(ends = "both"),linetype = 'dashed',
            aes(x = seq(from = 60.5+(qnorm(.45)-b0)/b1,to = 60+(qnorm(.45)-b0)/b1+4.5,by = .1), y = .45))+
  geom_line(arrow = arrow(ends = "both"),linetype = 'dashed',
            aes(x = seq(from = 60.5+(qnorm(.35)-b0)/b1,to = 60+(qnorm(.35)-b0)/b1+4.5,by = .1), y = .35))+
  geom_line(arrow = arrow(ends = "both"),linetype = 'dashed',
            aes(x = seq(from = 60.5+(qnorm(.25)-b0)/b1,to = 60+(qnorm(.25)-b0)/b1+4.5,by = .1), y = .25))+
  geom_line(arrow = arrow(ends = "both"),linetype = 'dashed',
            aes(x = seq(from = 60.5+(qnorm(.15)-b0)/b1,to = 60+(qnorm(.15)-b0)/b1+4.5,by = .1), y = .15))+
  geom_line(arrow = arrow(ends = "both"),linetype = 'dashed',
            aes(x = seq(from = 62.2+(qnorm(.05)-b0)/b1,to = 60+(qnorm(.05)-b0)/b1+4.5,by = .1), y = .05))

# extract data for the loess lines from the 'data' slot

gg4 <- ggplot_build(g4a+xlim(59,77))
df2 <- data.frame(x = gg4$data[[1]]$x,
                  ymin = gg4$data[[1]]$y,
                  ymax = gg4$data[[2]]$y) 

# use the loess data to add the 'ribbon' to plot 
g4d = g4a + geom_vline(xintercept = c(59,77),linetype = 'dotted',size = 1) +
  geom_ribbon(data = df2, aes(x = x, ymin = ymin, ymax = ymax),
              fill = "grey1", alpha = 0.4)

#g4e transformed data, panel data
l1 = function(t){b0+b1*(t-60)}
l2 = function(t){b0+b1*(t-65)}
age =c(65,66,67)
df1 = data.frame(x =age,y1 = qnorm(f1(age+rnorm(3,0,0.5))),y2 = qnorm(f2(age+rnorm(3,0,0.5))))
                 
g4e = 
  ggplot(df1)+
  geom_point(aes(x,y1))+
  geom_point(aes(x = x,y = y2))+
  xlim(62,72)+geom_smooth(data = df1,aes(x = x,y = y1),method = "lm", se = FALSE,color = myColors[1])+
  geom_smooth(data = df1,aes(x = x,y = y2),method = "lm", se = FALSE,color = myColors[2])+ylab("Transformed Y")+xlab("Age")+
  theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
  theme(axis.text.x = element_text(size = 15,face = 'bold',color = 'black'))+
  theme(axis.text.y = element_text(size = 15,face = 'bold',color = 'black'))
  

gg4e = ggplot_build(g4e)
df3 <- data.frame(x = gg4e$data[[3]]$x,
                  ymin = gg4e$data[[3]]$y,
                  ymax = gg4e$data[[4]]$y) 

g4f = g4e+ geom_ribbon(data = df3, aes(x = x, ymin = ymin, ymax = ymax),
                       fill = "grey1", alpha = 0.4)+ylim(-2,2)+ylab("Transformed Y")+xlab("Age")+
  theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
  theme(axis.text.x = element_text(size = 15,face = 'bold',color = 'black'))+
  theme(axis.text.y = element_text(size = 15,face = 'bold',color = 'black'))
  

png('graphic4a.png',width=800)
g4a
dev.off()

png('graphic4b.png',width=800)
g4b
dev.off()

png('graphic4c.png',width=800)
g4c
dev.off()

png('graphic4d.png',width=800)
g4d
dev.off()

png('graphic4e.png',width=800)
g4e+ylim(-2,2)
dev.off()

png('graphic4f.png',width=800)
g4f+ylim(-2,2)
dev.off()



###################################
###### Graphic 5 ##################
###################################
#Show ADNI clusters
#pull from ClusterSetup_ADNI.R
pF = function(f){
dd1 = f$fits$av45$subset
dd2 = f$fits$fdg$subset
dd3 = f$fits$hoc$subset
d1 = ggplot() 
d1+ 
  geom_line(data = dd1,aes(argvals + gamma0, y, colour = 'AV45',group = id)) + 
  geom_point(data = dd1,aes(argvals + gamma0, y,colour = 'AV45'), size = 1, alpha = 1)+ 
  geom_line(data = dd2,aes(argvals + gamma0, y, colour = 'FDG',group = id)) + 
  geom_point(data = dd2,aes(argvals + gamma0, y,colour = 'FDG'), size = 1, alpha = 1)+
  geom_line(data = dd3,aes(argvals + gamma0, y, colour = 'HOC',group = id)) + 
  geom_point(data = dd3,aes(argvals + gamma0, y,colour = 'HOC'), size = 1, alpha = 1)+
  theme(legend.position='bottom',plot.title = element_text(size = 20),legend.title=element_blank(),
  legend.text = element_text(size = 20),legend.key.size = unit(1.45, "cm"))+
  theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
  theme(axis.text.x = element_text(size = 15,face = 'bold',color = 'black'))+
  theme(axis.text.y = element_text(size = 15,face = 'bold',color = 'black'))+
  ylab("Biomarker Severity")+xlab("Adjusted Time")+scale_color_manual(values=myColors)
}


pF1 = function(f){
  dd1 = f$fits$av45$subset
  dd2 = f$fits$fdg$subset
  dd3 = f$fits$hoc$subset
  d1 = ggplot() 
  d1+ 
  geom_line(data = dd1,aes(argvals +gamma0, ghat, group = NULL,color = 'AV45'), size = 1.5)+
  geom_line(data = dd2,aes(argvals +gamma0, ghat, group = NULL, color = 'FDG'),size = 1.5)+
  geom_line(data = dd3,aes(argvals +gamma0, ghat, group = NULL, color = 'HOC'),size = 1.5)+
    theme(legend.position='bottom',plot.title = element_text(size = 20),legend.title=element_blank(),
          legend.text = element_text(size = 20),legend.key.size = unit(1.45, "cm"))+
    theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
    theme(axis.text.x = element_text(size = 15,face = 'bold',color = 'black'))+
    theme(axis.text.y = element_text(size = 15,face = 'bold',color = 'black'))+
    ylab("Biomarker Severity")+xlab("Adjusted Time")+scale_color_manual(values=myColors)
}


g6a = pF(gT1)
g6b = pF(gT2)  
g6c = pF(gT)

mylegend6<-g_legend(g6a)
png('graphic6a.png',width = 800)
g6A<- grid.arrange(arrangeGrob(g6a + theme(legend.position="none"),
                                g6b + theme(legend.position="none")+ylab(""),
                                g6c + theme(legend.position="none")+ylab(""),
                                nrow=1),
                    mylegend6,nrow = 2,heights=c(10, 1))
dev.off()


g6e = pF1(gT1)
g6f = pF1(gT2)  
g6g = pF1(gT3)

png('graphic6b.png',width = 800)
g6B<- grid.arrange(arrangeGrob(g6e + theme(legend.position="none"),
                               g6f + theme(legend.position="none")+ylab(""),
                               g6g + theme(legend.position="none")+ylab(""),
                               nrow=1),
                   mylegend6,nrow = 2,heights=c(10, 1))
dev.off()
