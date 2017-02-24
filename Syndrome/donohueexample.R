## Not run: 
library(grace)
library(fda)
setwd("/Users/Teresa/Documents/Dissertation/R_Dissertation/Syndrome")

set.seed(10)
x <- seq(-6, 6, by = 0.1)
f1 <- function(t){ 1 / (1 + exp(-t)) }
f2 <- function(t){ t / 12 + 0.5 }
f3 <- function(t){ (t+6)^2/144/2 }

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

grace.simulation.fits <- with(dd, grace(t, Y, Outcome, id, plots = TRUE))


g2 = with(subset(set4,Syndrome==1),grace(argvals = oTime,y = value,outcome = variable,id = RID,plots=TRUE))
g3 = with(subset(set4,Syndrome==2),grace(argvals = oTime,y = value,outcome = variable,id = RID,plots=TRUE))
g4 = with(set4,grace(argvals = oTime,y = value,outcome = variable,id = RID,group = factor(Syndrome),plots=TRUE))

m1 = with(m1,m1[order(RID,type,time1),])

gA1 = with(m1,grace(argvals = time1,y  = quant, outcome = type,id = RID,group = factor(clustB),plots = TRUE))
gA = with(m1,grace(argvals = time,y  = quant, outcome = type,id = RID,plots = TRUE))
gAS1 = with(subset(m1,clustB==1),grace(argvals = time1,y  = quant, outcome = type,id = RID,plots = TRUE))
gAS2 = with(subset(m1,clustB==2),grace(argvals = time1,y  = quant, outcome = type,id = RID,plots = TRUE))
gA3 = with(m1,grace(argvals = time,y  = quant, outcome = type,id = RID,group = factor(clustB),plots = TRUE))
gAS1 = with(subset(m1,clustB==1),grace(argvals = time1,y  = quant, outcome = type,id = RID,plots = TRUE))
gAS2 = with(subset(m1,clustB==2),grace(argvals = time1,y  = quant, outcome = type,id = RID,plots = TRUE))
gAS3 = with(subset(m1,clustB==3),grace(argvals = time1,y  = quant, outcome = type,id = RID,plots = TRUE))

#clustB
setwd("~/Documents/Davis2015-2016/AAIC")
gT = with(m1,graceT(argvals = time,y  = quant, outcome = type,id = RID,group = factor(clustB),plots = TRUE))
gT1 = with(subset(m1,clustB==1),graceT(argvals = time,y  = quant, outcome = type,id = RID,group = factor(clustB),plots = TRUE))
gT2 = with(subset(m1,clustB==2),graceT(argvals = time,y  = quant, outcome = type,id = RID,group = factor(clustB),plots = TRUE))
gT3 = with(subset(m1,clustB==3),graceT(argvals = time,y  = quant, outcome = type,id = RID,group = factor(clustB),plots = TRUE))

m1s = subset(m1, RID%in%c('2002','2007','2010'))
gTs = with(m1s,graceT(argvals = time,y  = quant, outcome = RID,id = type,group = factor(clustB),plots = TRUE))
gT1s = with(subset(m1,clustB==1),graceT(argvals = time,y  = quant, outcome = type,id = RID,group = factor(clustB),plots = TRUE))
gT2s = with(subset(m1,clustB==2),graceT(argvals = time,y  = quant, outcome = type,id = RID,group = factor(clustB),plots = TRUE))
gT3s = with(subset(m1,clustB==3),graceT(argvals = time,y  = quant, outcome = type,id = RID,group = factor(clustB),plots = TRUE))
