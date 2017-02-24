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
s80 = function(t){pnorm(b0+b1*(t-80)*I(t>=80))}

setwd("~/Documents/Dissertation/SyndromeClassification/AnnalsOfStat/figures")
lty1 = 2
lwd1 = 3
png('p1.png',width = 800,height = 800)
par(mfrow = c(2,2))
curve(s60,from = 60, to  = 80,n = 100,col = 'darkred', main = 'Person 1',xlab = 'Age',ylab = "Severity",lwd = lwd1)
curve(s65,from = 60, to  = 80,n = 100,col = 'darkred',lty = lty1,add = TRUE,lwd = lwd1)

curve(s75,from = 70, to  = 90,n = 100,col = 'darkblue', main = 'Person 2',xlab = 'Age',ylab = "Severity",lwd = lwd1)
curve(s80,from = 70, to  = 90,n = 100,col = 'darkblue',lty = lty1,add = TRUE,lwd = lwd1)

curve(s66,from = 60, to  = 80,n = 100,col = 'mediumseagreen',lty = lty1, main = 'Person 3',xlab = 'Age',ylab = "Severity",lwd = lwd1)
curve(s71,from = 60, to  = 80,n = 100,col = 'mediumseagreen',add = TRUE,lwd = lwd1)

curve(s75,from = 70, to  = 90,n = 100,col = 'darkorange', main = 'Person 4',xlab = 'Age',ylab = "Severity",lwd = lwd1)
curve(s70,from = 70, to  = 90,n = 100,col = 'darkorange',lty = lty1,add = TRUE,lwd = lwd1)
dev.off()

png('p2.png',width = 800,height = 800)
par(mfrow = c(2,1))
curve(s60,from = 60, to  = 90,n = 100,col = 'darkred', main = 'Person 1 and Person 2',xlab = 'Age',ylab = "Severity",lwd = lwd1)
curve(s65,from = 60, to  = 90,n = 100,col = 'darkred',lty = lty1,add = TRUE,lwd = lwd1)
curve(s75,from = 60, to  = 90,n = 100,col = 'darkblue', add = TRUE,xlab = 'Age',ylab = "Severity",lwd = lwd1)
curve(s80,from = 60, to  = 90,n = 100,col = 'darkblue',lty = lty1,add = TRUE,lwd = lwd1)

curve(s66,from = 60, to  = 90,n = 100,col = 'mediumseagreen',lty = lty1, main = 'Person 3 and Person 4',xlab = 'Age',ylab = "Severity",lwd = lwd1)
curve(s71,from = 60, to  = 90,n = 100,col = 'mediumseagreen',add = TRUE,lwd = lwd1)
curve(s75,from = 60, to  = 90,n = 100,col = 'darkorange', add = TRUE ,xlab = 'Age',ylab = "Severity",lwd = lwd1)
curve(s70,from = 60, to  = 90,n = 100,col = 'darkorange',lty = lty1,add = TRUE,lwd = lwd1)
dev.off()
