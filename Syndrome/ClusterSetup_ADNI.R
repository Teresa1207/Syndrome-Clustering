### Syndrome Classification    ###
### Model Setup and Simulation ###
# S2: Syndrome 2: Vascular-->HPCV, Amy is idependent. 
# In syndrome 1, Amyloid advances to the 50th percentile prior to HPCV followed by vascular. 
# In syndrome 2, Vascular advances to the 50th percentile prior to HPCV. Amyloid is indepenent. 
setwd("/Users/Teresa/Documents/Dissertation/R_Dissertation/Syndrome")
#########################################
library(mvtnorm)
library(MASS)
library(ggplot2)
library(GGally)
library(reshape2)
library(plyr)
library(geepack)

`%ni%` <- Negate(`%in%`)
###############################################################################
####### Cluster Design #######
###############################################################################
set.seed(12071983)
counter = 1
nRun = 10
errT1 = matrix(0,nrow = nRun,ncol = 2)
errT2 = matrix(0,nrow = nRun,ncol = 2)
errU = matrix(0,nrow = nRun,ncol = 2)
silsum = matrix(0,nrow = nRun,ncol = 3)

source('dataClus_ADNI.R')
set4 = mark.long #df.melt
model <- function(x){
  fit1 = try(lm(qnorm(quant) ~ time,data=x))
  #fit2 = try(lm(qnorm(value) ~ oTime+oAge, data=x))
  data.frame(qL=coef(fit1)[[1]], qS=coef(fit1)[[2]],qL1=coef(fit1)[[1]], qS1=coef(fit1)[[2]])
}

s4 = set4
s4$quant[which(s4$quant==0)]=0.0001
s4$quant[which(s4$quant==1)]=0.9999

#s4 = subset(s4,base3dx=='2:MCI')

clust = ddply(s4,.(RID,type),function(x){model(x)})

clustL = dcast(clust[,c("RID",'type',"qL","qS")],RID~type,value.var = c('qL'))
colnames(clustL) = c('RID','L1',"L2","L3",'L4')
clustS = dcast(clust[,c("RID",'type',"qL","qS")],RID~type,value.var = c('qS'))
colnames(clustS) = c('RID','S1',"S2","S3",'S4')
cCast = merge(clustL,clustS,by = 'RID',all = TRUE)

clustL1 = dcast(clust[,c("RID",'type',"qL1","qS1")],RID~type,value.var = c('qL1'))
colnames(clustL1) = colnames(clustL)
clustS1 = dcast(clust[,c("RID",'type',"qL1","qS1")],RID~type,value.var = c('qS1'))
colnames(clustS1) = colnames(clustS)
cCast1 = merge(clustL1,clustS1,by = 'RID',all = TRUE)



##Transform Betas##
trans = function(x){scale(x,scale=FALSE)}
cComp = cCast[complete.cases(cCast1),]
cData = data.frame(cComp$RID,t(data.frame(apply(cComp[,colnames(clustL)[-1]],1,scale))),
                    t(data.frame(apply(cComp[,colnames(clustS)[-1]],1,scale))))
colnames(cData) = c("RID",unlist(lapply(c(colnames(clustL)[-1],colnames(clustS)[-1]),function(x){paste(x,'t',sep='')})))
obs =3.5
cData1 = ddply(cComp,.(RID),summarise,
               b12 = (L1-L2)*(obs-1)+(S1 - S2)*(obs-1)^2/2,
               b13 = (L1-L3)*(obs-1)+(S1 - S3)*(obs-1)^2/2,
               b14 = (L1-L4)*(obs-1)+(S1 - S4)*(obs-1)^2/2,
               b23 = (L2-L3)*(obs-1)+(S2 - S3)*(obs-1)^2/2,
               b24 = (L2-L4)*(obs-1)+(S2 - S4)*(obs-1)^2/2,
               b34 = (L3-L4)*(obs-1)+(S3 - S4)*(obs-1)^2/2
)


### Calculate Average Lag for each person for each variable

##ADAS,AV45,FDG,HOC
#set4$QN = qnorm(set4$quant)
library(lme4)
#set4$QN[is.infinite(set4$QN)]=NA
#setlme = set4
#lmeClust = sapply(unique(setlme$RID),function(x){
#  ranef(lmer(QN~time+(time|type),data = subset(setlme,RID==x)))})

#rMat= function(t){matrix(c(t,t,t,0,0,0,-t,0,0,t,t,0,0,-t,0,-t,0,t,0,0,-t,0,-t,-t),nrow = 6)}
#rMatC = function(t){cbind(rMat(t),rMat(t^2/2))}
#clusts = data.frame(RID = unique(setlme$RID),t(sapply(lmeClust,function(x){
 # rMatC(2)%*%unlist(x)
#})))




#dist1 = dist(cData1[,-1])
vars = c('b12','b13','b14','b23','b24','b34')
#vars = c('b23','b24','b34')
#vars = c('b12','b13','b14')
cDat = subset(cData1,select = c('RID',vars))
#cDat = clusts
k1 = data.frame(RID = cDat$RID,clustB =kmeans(cDat[,-1],3,nstart=50)$cluster)#, clustT = kmeans(cData[,-1],3,nstart = 50)$cluster,clustU = kmeans(cCast1[,-1],2,nstart = 50)$cluster)
pam1 = pam(x = cDat[,-1],diss = FALSE,k = 2)
k2 = merge(cDat,k1,by = 'RID')
k2a = merge(k2,dc.u,by= 'RID')
k2a$pam = pam1$clustering
distAdni = dist(cDat[,-1])
hcAdni = hclust(distAdni,'single')
cutAdni <- cutree(hcAdni, 3)
k2a$hc = cutAdni
#k2a = merge(k2,cCast1,by = 'RID')
### Scatter plot
myGrey = c('black','grey35','grey55','grey5')

p1 = ggplot(k2a,aes(x = b12,y = b13))+
  geom_point(aes(shape = factor(hc),color  = factor(hc)),size = 4)+
  theme(legend.position="none")+geom_vline(xintercept = 0)+geom_hline(yintercept = 0)+ 
  xlab(expression('Amyloid - FDG,  Increasing Lag' %->% ''))+ylab(expression('FDG - HOC, Increasing Lag' %->% ''))+#scale_color_manual(values=myColors)+
  theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
  theme(axis.text.x = element_text(size = 15,face = 'bold',color = 'black'))+
  theme(axis.text.y = element_text(size = 15,face = 'bold',color = 'black'))
#setwd("/Users/Teresa/Documents/Davis2015-2016/AAIC/figures")
#png('graphic5.png')
p1
#dev.off()


setwd("/Users/Teresa/Documents/Davis2015-2016/AAIC")

ggplot(k2a,aes(x = b23,y = b34,color = factor(clustB)))+geom_point()

png('clustersAdni.png',width=800)
p 
dev.off()




sum1 = ddply(k2a,.(clustB),summarise,
            d12 = mean(b12),
            d13 = mean(b13),
            d14 = mean(b14),
            d23 = mean(b23),
            d24 = mean(b24),
            d34 = mean(b34))
sum2 = ddply(k2a,.(clustB),summarise,
             l1 = mean(L1),
             l2 = mean(L2),
             l3 = mean(L3),
             l4 = mean(L4),
             s1 = mean(S1),
             s2 = mean(S2),
             s3 = mean(S3),
            s4 = mean(S4))


e.lag1 = with(k2a,data.frame(cID = clustB,lg1=mapply(function(x,y){-x/y},x = L1,y = S1)))
e.lag2 = with(k2a,data.frame(cID = clustB,lg2=mapply(function(x,y){-x/y},x = L2,y = S2)))
e.lag3 = with(k2a,data.frame(cID = clustB,lg3=mapply(function(x,y){-x/y},x = L3,y = S3)))
e.lag4 = with(k2a,data.frame(cID = clustB,lg4=mapply(function(x,y){-x/y},x = L4,y = S4)))

lg12 = data.frame(cID = e.lag1$cID,lg12 = e.lag2$lg2 - e.lag1$lg1)
lg13 = data.frame(cID = e.lag1$cID,lg13 =e.lag3$lg3 - e.lag1$lg1)
lg14 = data.frame(cID = e.lag1$cID,lg14 =e.lag4$lg4 - e.lag1$lg1)
lg23 = data.frame(cID = e.lag1$cID,lg23 =e.lag3$lg3 - e.lag2$lg2)
lg24 = data.frame(cID = e.lag1$cID,lg24 =e.lag4$lg4 - e.lag2$lg2)
lg34 = data.frame(cID = e.lag1$cID,lg34 =e.lag4$lg4 - e.lag3$lg3)

lags = data.frame(lg23,lg24=lg24$lg24,lg34 =lg34$lg34)
lags = lags[complete.cases(lags),]
lags = subset(lags,abs(lg23)<20&abs(lg24)<20)
lags = subset(lags,abs(lg34)<20)

library(plotly)
plot_ly(data = lags, x = lg24,y = lg34, mode = "markers",
        color = cID)

ddply(lg12[!is.infinite(lg12$lg12),],'cID',summarise,
      lag12 = mean(lg12))
ddply(lg13[!is.infinite(lg13$lg13),],'cID',summarise,
      lag13 = mean(lg13))
ddply(lg14[!is.infinite(lg14$lg14),],'cID',summarise,
      lag14 = mean(lg14))
ddply(lags,'cID',summarise,
      lag23 = mean(lg23),
      lag24 = mean(lg24),
      lag34 = mean(lg34))
myColors <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

m1 = merge(subset(s4,RID%in%k1$RID),k1,by = "RID",all.x = TRUE)
saveRDS(m1,file = 'm1.RDA')

m1$time1 = round_any(m1$time,.25)

#m2 = merge(subset(set3,RID%in%k1$RID),k1,by = 'RID',all.x = TRUE)
m=m2
ggplot(m,aes(x = M1,y = M3,group = RID,color = factor(clustB)))+geom_line()+facet_grid(~clustT)
ggplot(m,aes(x = M1,y = M3,group = RID,color = factor(clustB)))+geom_line()+
  stat_smooth(aes(group = Syndrome),color = 'black')

s = m1
ggplot(s, aes(x = time,y = quant,group=c(RID),color = type))+ggtitle('Total')+geom_line()+facet_grid(~type)
  geom_line(aes(subset(s,variable=='M1')$dAge,subset(s,variable=='M1')$ve,group=subset(s,variable=='M1')$RID),color='darkred')
  geom_line(aes(subset(s,variable=='M2')$dAge,subset(s,variable=='M2')$ve,group=subset(s,variable=='M2')$RID),color='darkgreen')+
  geom_line(aes(subset(s,variable=='M3')$dAge,subset(s,variable=='M3')$ve,group=subset(s,variable=='M3')$RID),color='darkblue')+
  stat_smooth(aes(group = variable),size = 2)

s = subset(m1,clustB==1)
ggplot()+ggtitle('Clust Groups')+xlim(0,7)+
  geom_line(data = subset(s,type=='adas'),aes(time,quant,group=RID,color = 'ADAS'))+
  geom_line(data = subset(s,type=='av45'),aes(time,quant,group=RID,color = 'AV45'))+
  geom_line(data = subset(s,type=='fdg'),aes(time,quant,group=RID,color = 'FDG'))+
  geom_line(data = subset(s,type=='hoc'),aes(time,quant,group=RID,color = 'HOC'))+
  stat_smooth(data = subset(s,type=='adas'),aes(group = NULL),size = 2)+scale_color_manual(values=myColors)
s = subset(m1,clustB==2)
ggplot()+ggtitle('Clust Groups')+xlim(0,7)+
  geom_line(data = subset(s,type=='adas'),aes(time,quant,group=RID,color = 'ADAS'))+
  geom_line(data = subset(s,type=='av45'),aes(time,quant,group=RID,color = 'AV45'))+
  geom_line(data = subset(s,type=='fdg'),aes(time,quant,group=RID,color = 'FDG'))+
  geom_line(data = subset(s,type=='hoc'),aes(time,quant,group=RID,color = 'HOC'))+
  stat_smooth(data = subset(s,type=='adas'),aes(group = NULL),size = 2)+scale_color_manual(values=myColors)


s = m1
ggplot(s, aes(x = tAge,y = ve,group=c(RID),col = variable))+ggtitle('Total')+
  geom_line(aes(subset(s,variable=='M1')$tAge,subset(s,variable=='M1')$ve,group=subset(s,variable=='M1')$RID),color='darkred')+
  geom_line(aes(subset(s,variable=='M2')$tAge,subset(s,variable=='M2')$ve,group=subset(s,variable=='M2')$RID),color='darkgreen')+
  geom_line(aes(subset(s,variable=='M3')$tAge,subset(s,variable=='M3')$ve,group=subset(s,variable=='M3')$RID),color='darkblue')+
  stat_smooth(aes(group = variable),size = 2)

ggplot(m,aes(x = M2,y = M3,group = RID,color = factor(Syndrome)))+geom_line()+
  stat_smooth(aes(group = factor(Syndrome)), color = 'black')

ggplot(m,aes(x = M2,y = M3,group = RID,color = factor(clustB)))+geom_line()+
  stat_smooth(aes(group = factor(clustB)), color = 'black')

library(plotly)

f <- list(
  family = "Courier New, monospace",
  size = 20,
  color = "black"
)
x <- list(
  title = "Amyloid/FDG",
  titlefont = f
)
y <- list(
  title = "FDG/HOC",
  titlefont = f
)
z <- list(
  title = "Amyloid/HOC",
  titlefont = f
)
p1 = plot_ly(data = k2a, x = b23,symbol = factor(clustB),color = NULL,colors = cols,name = 'ADNI Clusters',mode = "markers")%>%
  layout(xaxis = list(title = "Amyloid/FDG")) %>% layout(showlegend = FALSE)
p1
p2 = plot_ly(data = k2a, x = b24,name = 'ADNI Clusters',mode = "markers",color =factor(clustB))%>%
  layout(xaxis = list(title = "Amyloid/HOC")) %>% layout(showlegend = FALSE)
p2

p3 = plot_ly(data = k2a, x = b34,name = 'ADNI Clusters',mode = "markers",color =factor(clustB))%>%
  layout(xaxis = list(title = 'FDG/HOC')) %>% layout(showlegend = FALSE)
p3
p = plot_ly(data = k2a, x = b23,y = b34,name = 'ADNI Clusters',mode = "markers",color =factor(clustB))%>%
  layout(xaxis = x, yaxis = y) %>% layout(showlegend = FALSE)
p
p4 = plot_ly(data = k2a, x = b23,y = b24,name = 'ADNI Clusters',mode = "markers",color =factor(clustB))%>%
  layout(xaxis = list(title = "FDG lag"), yaxis =list(title = "HOC lag")) %>% layout(showlegend = FALSE)
p4
#source('plotly_infor.R')
setwd("/Users/Teresa/Documents/Davis2015-2016/AAIC")

plotly_IMAGE(p1, format = "png", out_file = "graphic5a.png",width = 1000,height = 800)
plotly_IMAGE(p2, format = "png", out_file = "graphic5b.png",width = 1000,height = 800)
plotly_IMAGE(p3, format = "png", out_file = "graphic5c.png",width = 1000,height = 800)
plotly_IMAGE(p, format = "png", out_file = "graphic5.png",width = 1000,height = 800)
plotly_IMAGE(p4, format = "png", out_file = "graphic6.png",width = 1000,height = 800)
