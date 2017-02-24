################################################################################
######    Syndrome Classifiation  		
######	  By: Teresa Filshtein
######	  Date: March 14, 2016			                      
######	  Modified: 
######    EDITED BY: 
################################################################################
setwd("~/Documents/Dissertation/R_Dissertation/Syndrome")
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


f1 = function(t){pnorm(b0+b1*(t-0)*I(t>=0))}
f2 = function(t){pnorm(b0+b1*(t-5)*I(t>=5))}
f3 = function(t){pnorm(b0+b1*(t-10)*I(t>=10))}
f4 = function(t){pnorm(b0+b1*(t-15)*I(t>=15))}

setwd("~/Documents/Davis2015-2016/AAIC")
##ggplotly

df = data.frame(DiseaseTime = seq(-5,30,by = .1))
png('syn1.png',height = 500, width = 800)
base = ggplot(df,aes(DiseaseTime))
base+stat_function(fun = f1, aes(color = 'Marker 1'),size =1.2)+
  stat_function(fun = f2,aes(color = 'Marker 2'),size =1.2)+
  stat_function(fun = f3,aes(color = 'Marker 3'),size =1.2)+
 theme(legend.title=element_blank())+ylab("Biomarker Severity")+xlab("Disease Time")+
  ggtitle("Syndrome I")
dev.off()

png('syn2.png',height = 500, width = 800)

base+stat_function(fun = f3, aes(color = 'Marker 1'),size =1.2)+
  stat_function(fun = f2,aes(color = 'Marker 2'),size =1.2)+
  stat_function(fun = f1,aes(color = 'Marker 3'),size =1.2)+
  theme(legend.title=element_blank())+ylab("Biomarker Severity")+xlab("Disease Time")+
  ggtitle("Syndrome II")
dev.off()
ggplotly()
