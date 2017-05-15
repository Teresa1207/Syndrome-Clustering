## Error Analysis
results1 = readRDS('~/Documents/Dissertation/R_Dissertation/Syndrome/Server/ErrorEstServer/res1.RDA')
results2 = readRDS('~/Documents/Dissertation/R_Dissertation/Syndrome/Server/ErrorEstServer/res2.RDA')
results3 = readRDS('~/Documents/Dissertation/R_Dissertation/Syndrome/Server/ErrorEstServer/res3.RDA')

head(results3[[1]])

tab_res =function(resultsX){data.frame(
  SD = seq(from = 0, to = .5, by = 0.01),
  Syn = unlist(lapply(resultsX, function(x){
  mean(x$TRand,na.rm = TRUE
       )
})),
  M1 = unlist(lapply(resultsX, function(x){
    mean(x$S1Rand,na.rm = TRUE
    )
  })),
  M2 = unlist(lapply(resultsX, function(x){
  mean(x$S2Rand,na.rm = TRUE
  )
})),
  M3 = unlist(lapply(resultsX, function(x){
  mean(x$S3Rand,na.rm = TRUE
  )
})))
}

r1 = tab_res(results1)
r2 = tab_res(results2)
r3  = tab_res(results3)

png('~/Documents/Dissertation/Dissertation_2017/figures/sdPlot.png',height = 1000, width = 1500)
ggplot(r1, aes(x = SD))+ coord_cartesian(xlim = c(0,.5),ylim = c(.7,1))+theme_classic() + 
  geom_line(data = r1,aes(x = SD,y = Syn,linetype = '1:1:1'),size = 3)+
  geom_line(data = r2,aes(x = SD,y = Syn,linetype = '2:1:1'),size = 3)+
  geom_line(data = r3,aes(x = SD,y = Syn,linetype = '5:4:1'),size = 3)+
  scale_linetype_manual(values = c('solid','twodash','dotted'), guide_legend(title = '')) +
    xlab('Standard Deviation') + ylab('Average Rand Index')+
   guides(linetype = guide_legend(override.aes = list(size=1)))+
    theme(legend.text = element_text(size = 24),
          legend.position = 'bottom',
          axis.text.x = element_text(size = 24,face = 'bold'),
          axis.text.y = element_text(size = 24,face = 'bold'),
          axis.title = element_text(size = 24,face = 'bold')
          )
dev.off()

  
  
