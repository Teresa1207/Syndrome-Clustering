## Simulation run, Gap Statistic
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(Rmisc)
library(gridExtra)
library(grid)
res18 = readRDS('~/Documents/Dissertation/R_Dissertation/Syndrome/Server/ServerGap/resGap18.RDA')
res11 = readRDS('~/Documents/Dissertation/R_Dissertation/Syndrome/Server/ServerGap/resGap11.RDA')
res31 = readRDS('~/Documents/Dissertation/R_Dissertation/Syndrome/Server/ServerGap/resGap31.RDA')
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
colP = cbbPalette[c(2,4,6)]
pRes = c(3,3,3,2,2,2,1)

res18df = do.call(rbind, lapply(res18, data.frame, stringsAsFactors=FALSE))
res11df = do.call(rbind, lapply(res11, data.frame, stringsAsFactors=FALSE))
res31df = do.call(rbind, lapply(res31, data.frame, stringsAsFactors=FALSE))

#colnames(res18df)=colnames(res31df)=colnames(res11df) =paste('G',1:7,sep="")
dfFunc = function(df){
df$type1 = c(rep('3 groups',3*500), rep('2 groups',3*500),rep('1 group',500))
df$type2 = c(rep('1:1:1',500),rep('2:1:1',500),rep('5:4:1',500),
                  rep('1:1:0',500),rep('2:1:0',500),rep('9:1:0',500),
                  rep('1:0:0',500))
df$type2 = factor(df$type2, levels = c('1:0:0','1:1:0','2:1:0','9:1:0','1:1:1','2:1:1','5:4:1'), labels = c('1 Group - 1:0:0','2 Groups - 1:1:0','2 Groups - 2:1:0','2 Groups - 9:1:0','3 Groups - 1:1:1','3 Groups - 2:1:1','3 Groups - 5:4:1'))
df$true = c(rep(3,3*500), rep(2,3*500),rep(1,500))
df}

pFunc = function(df){ggplot(df,aes(clus))+theme_classic()+geom_bar(aes(fill=type2),position="dodge")+
    facet_grid(~type1)+ ylab('Count')+xlab("Number of groups found")+
    scale_fill_brewer(palette = "Set1")+theme(legend.position = 'top',legend.title = element_blank(),axis.text.x = element_text(size = 14,face = 'bold'),
                                              axis.text.y = element_text(size = 14,face = 'bold'),
                                              axis.title.x = element_text(size = 16, face = 'bold'),
                                              axis.title.y = element_text(size  = 16, face = 'bold'),
                                              legend.spacing = unit(1,"cm"),
                                              legend.key.size = unit(1,"cm"),
                                              legend.text = element_text(size = 16),
                                              strip.text.x = element_text(size = 16, face = 'bold'))
                      }
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

p1 = pFunc(dfFunc(res18df))
p2 = pFunc(dfFunc(res11df))
p3 = pFunc(dfFunc(res31df))

g_leg = g_legend(p1)
png('~/Documents/Dissertation/Dissertation_2017/figures/groups.png',height = 800, width = 1000)
grid.arrange(arrangeGrob(p2 + theme(legend.position="none")+xlab("")+ggtitle('Low between/Low within error')+theme(title = element_text(size = 16)),
                         p3 + theme(legend.position="none")+xlab("")+ggtitle('Low between/high within error')+theme(title = element_text(size = 16)),
                         p1 + theme(legend.position="none")+ggtitle('Moderate between/Moderate within error')+theme(title = element_text(size = 16)),
                         nrow=3),
            g_leg,nrow = 2,heights=c(15,3),widths = c(15))
dev.off()


lFunc  = function(df){lapply(1:7,function(x){
  length(which(df[,x]==pRes[x]))/500
})
}

lFunc(res18df)
