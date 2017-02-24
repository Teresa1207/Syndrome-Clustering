## Not run: 
library(grace)
library(fda)
library(mvtnorm)
setwd("/Users/Teresa/Documents/Dissertation/R_Dissertation/Syndrome")

set.seed(10)
x <- seq(-6, 6, by = 0.1)
f1 <- function(t){ 1 / (1 + exp(-t)) }
f2 <- function(t){ f1(t+5)}
f3 <- function(t){ f1(t -5)}

pdd <- rbind(
  data.frame(x = x, y = f1(x), curve = "f1(x)"),
  data.frame(x = x, y = f2(x), curve = "f2(x)"),
  data.frame(x = x, y = f3(x), curve = "f3(x)"))  

qplot(x, y, data = pdd, colour = curve, linetype = curve, geom = "line")

t <- seq(-1, 1, by = 0.5)
n <- 100
id <- 1:n
gamma0 <- runif(n, min=-5, max=5)
dd <- cbind(id, gamma0)
dd <- merge(dd, data.frame(t))
sig <- 0.01
rho <- 0.005

dd <- data.frame(rbind(
  cbind(dd, Outcome = 1, rmvnorm(n, c(0,0), matrix(c(sig, rho, rho, sig), nrow = 2))),
  cbind(dd, Outcome = 2, rmvnorm(n, c(0,0), matrix(c(sig, rho, rho, sig), nrow = 2))),
  cbind(dd, Outcome = 3, rmvnorm(n, c(0,0), matrix(c(sig, rho, rho, sig), nrow = 2)))))
colnames(dd) <- c("id", "gamma0", "t", "Outcome", "alpha0", "alpha1")
dd$e <- rnorm(nrow(dd), sd = 0.1)


dd$Y0 <- with(dd, 
              ifelse(Outcome == 1, f1(t + gamma0),
                     ifelse(Outcome == 2, f2(t + gamma0),
                            f3(t + gamma0))))
dd$Y <- with(dd, Y0 + alpha0 + alpha1 * t + e)
dd <- dd[with(dd, order(id, Outcome, t)), ]
dd$Outcome <- paste("Outcome", dd$Outcome)

dd1 = dd
dd1$group = 'g1'
dd2 = dd1
dd2$group = 'g2'
dd2$Outcome2 = NA
dd2$Outcome2[which(dd2$Outcome=='Outcome 1')]='Outcome 3'
dd2$Outcome2[which(dd2$Outcome=='Outcome 3')]='Outcome 1'
dd2$Outcome2[which(dd2$Outcome=='Outcome 2')]='Outcome 2'
dd2$Outcome=dd2$Outcome2
dd2 = dd2[,-ncol(dd2)]
dd2$id = dd2$id+100

ddT = rbind(dd1,dd2)
samp = sample(1:200,100, replace = FALSE)
ddS = subset(ddT,id%in%samp) 
grace.simulation.fits <- with(ddS, graceT(t, Y, Outcome, id, plots = TRUE))
grace.simulation.fits1 <- with(dd1, graceT(t, Y, Outcome, id,group = group, plots = TRUE))
grace.simulation.fits2 <- with(dd2, graceT(t, Y, Outcome, id,group = group, plots = TRUE))

fit1 = grace.simulation.fits$fits$`Outcome 1`$subset
fit2 = grace.simulation.fits$fits$`Outcome 2`$subset
fit3 = grace.simulation.fits$fits$`Outcome 3`$subset


fit11 = grace.simulation.fits1$fits$`Outcome 1`$subset
fit12 = grace.simulation.fits1$fits$`Outcome 2`$subset
fit13 = grace.simulation.fits1$fits$`Outcome 3`$subset

fit21 = grace.simulation.fits2$fits$`Outcome 1`$subset
fit22 = grace.simulation.fits2$fits$`Outcome 2`$subset
fit23 = grace.simulation.fits2$fits$`Outcome 3`$subset

#clustB
setwd("~/Documents/Davis2015-2016/AAIC")

gT = with(m1Sim,graceT(argvals = oTime,y  = ve, outcome = variable,id = RID,group = factor(clustB),plots = FALSE))
gT1 = with(subset(m1Sim,clustB==1),graceT(argvals = oTime,y  = ve, outcome = variable,id = RID,group = factor(clustB),plots = FALSE))
gT2 = with(subset(m1Sim,clustB==2),graceT(argvals = oTime,y  = ve, outcome = variable,id = RID,group = factor(clustB),plots = FALSE))
gT3 = with(subset(m1Sim,clustB==2),graceT(argvals = oTime,y  = ve, outcome = variable,id = RID,group = factor(clustB),plots = FALSE))

fitG = with(gT,rbind(fits$M1$subset,fits$M2$subset,fits$M3$subset))
fitG1 = with(gT1,rbind(fits$M1$subset,fits$M2$subset,fits$M3$subset))
fitG2 = with(gT2,rbind(fits$M1$subset,fits$M2$subset,fits$M3$subset))
fitG3 = with(gT3,rbind(fits$M1$subset,fits$M2$subset,fits$M3$subset))

c2Tsim = ddply(fitG,'outcome',function(x){sum(x$smooth.resids^2)})
c2g1g2sim =ddply(rbind(fitG1,fitG2),'outcome',function(x){sum(x$smooth.resids^2)})



pF1 = function(f){
  dd1 = f$fits$M1$subset
  dd2 = f$fits$M2$subset
  dd3 = f$fits$M3$subset
  d1 = ggplot() 
  d1+ 
    geom_line(data = dd1,aes(argvals +gamma0, ghat, group = NULL,color = 'M1'), size = 1.5)+
    geom_line(data = dd2,aes(argvals +gamma0, ghat, group = NULL, color = 'M2'),size = 1.5)+
    geom_line(data = dd3,aes(argvals +gamma0, ghat, group = NULL, color = 'M3'),size = 1.5)+
    theme(legend.position='bottom',plot.title = element_text(size = 20),legend.title=element_blank(),
          legend.text = element_text(size = 20),legend.key.size = unit(1.45, "cm"))+
    theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
    theme(axis.text.x = element_text(size = 15,face = 'bold',color = 'black'))+
    theme(axis.text.y = element_text(size = 15,face = 'bold',color = 'black'))+
    ylab("Biomarker Severity")+xlab("Adjusted Time")+scale_color_manual(values=myColors)
}


g6a = pF1(gT1)
g6b = pF1(gT2)  
g6c = pF1(gT)
fitG = with(gT,rbind(fits$adas$subset,fits$av45$subset,fits$fdg$subset,fits$hoc$subset))
fitG1 = with(gT1,rbind(fits$adas$subset,fits$av45$subset,fits$fdg$subset,fits$hoc$subset))
fitG2 = with(gT2,rbind(fits$adas$subset,fits$av45$subset,fits$fdg$subset,fits$hoc$subset))
fitG3 = with(gT3,rbind(fits$adas$subset,fits$av45$subset,fits$fdg$subset,fits$hoc$subset))

### HC

distSim = dist(cData1[,-1],method = 'euclidean')
clusters <- hclust(distSim)
plot(clusters)
