
setwd("~/Documents/Dissertation/R_Dissertation/Syndrome")

a1 = fdgM[,c('RID','quant','time')]
a2 = adasM[,c('RID','quant','time')]
a3 = hocM[,c('RID','quant','time')]
colnames(a1) = c('RID','fdg','time')
colnames(a2) = c('RID','adas','time')
colnames(a3) = c('RID','hoc','time')

a1$t2 = round_any(a1$time,.25)
a2$t2 = round_any(a2$time,.25)
a3$t2 = round_any(a3$time,.25)

a12=merge(a1[,-3],a2[,-3],by = c('RID','t2'),all=TRUE)
a123 = merge(a12,a3[,-3],by = c('RID','t2'),all=TRUE)
a23 = merge(a2[,-3],a3[,-3],by = c('RID','t2'),all=TRUE)

a.com = a23[complete.cases(a23),]
a.com = a.com[order(a.com$RID,a.com$t2),]
b=3
cData2 = ddply(M1cast,.(RID),summarise,
               b12 = (aL-hL)*(b)+(aS - hS)*(b)^2/2,
               b13 = (aL-hL)*(b)+(aS - hS)*(b)^2/2,
               b23 = (hL-fL)*(b)+(hS - fS)*(b)^2/2
)



k1 = data.frame(RID = d1$RID,clustR = kmeans(cData2[,-c(1,4)],3)$cluster,clustU = kmeans(d1[,-1],3)$cluster,clustT = kmeans(d2[,-1],3)$cluster)
r = kmeans(cData2[,-c(1,4)],3)
a.m = merge(a.com,k1,by = 'RID',all.x = TRUE)
a.m = merge(a.m,unique(dc[,c("RID","age65","base3dx")]),by = 'RID',all.x = TRUE)

ggplot(a.m,aes(adas,hoc,group = RID,color = factor(clustR)))+geom_line()

ggplot(subset(a.m,!is.na(clustU)),aes(x = adas,y = hoc,group = RID,color = factor(clustU)))+geom_line()+
  stat_smooth(aes(group = factor(clustU)), color = 'black')

ggplot(subset(a.m,!is.na(clustT)),aes(x = adas,y = hoc,group = RID,color = factor(clustT)))+geom_line()+
  stat_smooth(aes(group = factor(clustT)), color = 'black')

ggplot(subset(a.m,!is.na(clustR)),aes(x = adas,y = hoc,group = RID,color = factor(clustR)))+geom_line()+
  stat_smooth(aes(group = factor(clustR)), color = 'black')


