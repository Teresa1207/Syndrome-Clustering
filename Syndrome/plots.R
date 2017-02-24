### General Plots for Syndrome Classification Paper ###
t= seq(0,30,by = .001)
sig = .01^2
sL = .1 #start level
on = .5 #onset
iT = 10 #time for first marker
l1 = 5
l2 = 3
B0 = qnorm(sL)
B1 = (qnorm(on)-qnorm(sL))/iT
f1 = function(t){pnorm(B0+B1*t)}
f2 = function(t){pnorm(B0+B1*(t-l1))}
f3 = function(t){pnorm(B0+B1*(t-l1-l2))}
f1I = function(t){B0+B1*t}
f2I = function(t){B0+B1*(t-l1)}
f3I = function(t){B0+B1*(t-l1-l2)}

h1 = function(t){pnorm(B0+B1*(t+2))}
h2 = function(t){pnorm(B0+B1*(t-l1))}
h3 = function(t){pnorm(B0+B1*(t))}
h1I = function(t){B0+B1*(t+2)}
h2I = function(t){B0+B1*(t-l1)}
h3I = function(t){B0+B1*(t)}


par(mfrow = c(1,2))
curve(f1,from = -5,to = 30, col = 'darkred',lwd = 3, main = "Syndrome I",ylab = 'Severity',xlab = 'Disease Age')
curve(f2,from = -5,to = 30, add = TRUE,col = 'darkgreen',lwd = 3)
curve(f3,from = -5,to = 30, add = TRUE,col = 'darkblue',lwd = 3)
curve(f1I,from = -5,to = 30, col = 'darkred',lwd = 3, main = "Syndrome I",ylab = 'Linear Predictor',xlab = 'Disease Age')
curve(f2I,from = -5,to = 30, add = TRUE,col = 'darkgreen',lwd = 3)
curve(f3I,from = -5,to = 30, add = TRUE,col = 'darkblue',lwd = 3)

dev.off()
library(gridExtra)
g = ggplot(data.frame(Age = c(-5,30)),aes(Age))+ylab('Severity')+xlab("Disease Time")
ga = ggplot(data.frame(Age = c(-5,30)),aes(Age))+ylab('Linear Predictor')+xlab("Disease Time")

g1 = g+stat_function(fun = f1, colour = "darkred") +
  stat_function(fun = f2, colour = "darkgreen")+
  stat_function(fun = f3, colour = 'darkblue')
# build plot object for rendering 
gg1 <- ggplot_build(g1)

# extract data for the loess lines from the 'data' slot
df1 <- data.frame(x = gg1$data[[1]]$x,
                  ymin = gg1$data[[1]]$y,
                  ymax = gg1$data[[2]]$y) 
df1a <- data.frame(x = gg1$data[[1]]$x,
                  ymin = gg1$data[[2]]$y,
                  ymax = gg1$data[[3]]$y)
g1a = g1 +
  geom_ribbon(data = df1[43:59,], aes(x = x, ymin = ymin, ymax = ymax),
              fill = "#333333", alpha = 0.4)+
  geom_ribbon(data = df1a[43:59,], aes(x = x, ymin = ymin, ymax = ymax),
              fill = "#CC9933", alpha = 0.4)
g2 = ga+stat_function(fun = f1I, colour = "darkred") +
  stat_function(fun = f2I, colour = "darkgreen")+
  stat_function(fun = f3I, colour = 'darkblue')
# build plot object for rendering 
gg2 <- ggplot_build(g2)

# extract data for the loess lines from the 'data' slot
df2 <- data.frame(x = gg2$data[[1]]$x,
                  ymin = gg2$data[[1]]$y,
                  ymax = gg2$data[[2]]$y) 
df2a <- data.frame(x = gg2$data[[1]]$x,
                  ymin = gg2$data[[2]]$y,
                  ymax = gg2$data[[3]]$y)
g2a = g2 +
  geom_ribbon(data = df2[43:59,], aes(x = x, ymin = ymin, ymax = ymax),
              fill = "#333333", alpha = 0.4)+
  geom_ribbon(data = df2a[43:59,], aes(x = x, ymin = ymin, ymax = ymax),
              fill = "#CC9933", alpha = 0.4)

g3 = g+stat_function(fun = h1, colour = "darkred") +
  stat_function(fun = h2, colour = "darkgreen")+
  stat_function(fun = h3, colour = 'darkblue')
# build plot object for rendering 
gg3 <- ggplot_build(g3)

# extract data for the loess lines from the 'data' slot
df3 <- data.frame(x = gg3$data[[1]]$x,
                  ymin = gg3$data[[1]]$y,
                  ymax = gg3$data[[2]]$y) 
df3a <- data.frame(x = gg3$data[[1]]$x,
                  ymin = gg3$data[[2]]$y,
                  ymax = gg3$data[[3]]$y) 
g3a = g3 +
  geom_ribbon(data = df3[43:59,], aes(x = x, ymin = ymin, ymax = ymax),
              fill = "#333333", alpha = 0.4)+
  geom_ribbon(data = df3a[43:59,], aes(x = x, ymin = ymin, ymax = ymax),
              fill = "#CC9933", alpha = 0.4)



g4 = ga+stat_function(fun = h1I, colour = "darkred") +
  stat_function(fun = h2I, colour = "darkgreen")+
  stat_function(fun = h3I, colour = 'darkblue')
# build plot object for rendering 
gg4 <- ggplot_build(g4)

# extract data for the loess lines from the 'data' slot
df4 <- data.frame(x = gg4$data[[1]]$x,
                  ymin = gg4$data[[1]]$y,
                  ymax = gg4$data[[2]]$y) 
df4a <- data.frame(x = gg4$data[[1]]$x,
                  ymin = gg4$data[[2]]$y,
                  ymax = gg4$data[[3]]$y)
g4a = g4 +
  geom_ribbon(data = df4[43:59,], aes(x = x, ymin = ymin, ymax = ymax),
              fill = "#333333", alpha = 0.4)+
  geom_ribbon(data = df4a[43:59,], aes(x = x, ymin = ymin, ymax = ymax),
              fill = "#CC9933", alpha = 0.4)

setwd("~/Documents/Dissertation/SyndromeClassification/AnnalsOfStat")
png('plot1.png',height = 800, width = 800)
grid.arrange(g1a,g2a,g3a,g4a,top=textGrob("Area between Marker trajectories, Top: Syndrome 1, Bottom: Syndrome 2",gp=gpar(fontsize=10,font=3)))
dev.off()
