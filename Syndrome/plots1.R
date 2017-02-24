
setwd("/Users/Teresa/Documents/Dissertation/R_Dissertation/Syndrome")

library(plotly)
source('ModelSetup3.R')
source('DataSetup2.R')
t112 = set1$t2 - set1$t1
t113 = set1$t3 - set1$t1
t123 = set1$t3 - set1$t2

t212 = set2$t2 - set2$t1
t213 = set2$t3 - set2$t1
t223 = set2$t3 - set2$t2

names_ly = c('Lag Between M1 and M2','Lag Between M1 and M3','Lag Between M2 and M3')
names_ly1 = c('Onset M1','Onset M2','Onset M3')

colors_ly = c('blue','orange','green','red','violet')
colors_ly = c('blue','orange','green','red','violet')

axis_template = list(showgrid = F,
                     zeroline = F,
                     nticks = 20,
                     showline = T,
                     title = 'Estimated Slope',
                     mirror = 'all')

p1 = plot_ly(x = t112, opacity = 0.6, type = "histogram",name = names_ly[1],marker = list(color = colors_ly[1])) %>%
  add_trace(x = t113, opacity = 0.6, type = "histogram",name = names_ly[2],marker = list(color = colors_ly[2])) %>%
  add_trace(x = t123, opacity = 0.6, type = "histogram",name = names_ly[3],marker = list(color = colors_ly[3])) %>%
  add_trace(x = t212, opacity = 0.6, type = "histogram",name = names_ly[1],marker = list(color = colors_ly[1]),showlegend = FALSE) %>%
  add_trace(x = t213, opacity = 0.6, type = "histogram",name = names_ly[2],marker = list(color = colors_ly[2]),showlegend = FALSE) %>%
  add_trace(x = t223, opacity = 0.6, type = "histogram",name = names_ly[3],marker = list(color = colors_ly[3]),showlegend = FALSE) %>%
  layout(title = 'Histogram of Lag Times for Syndrome I and II',barmode="overlay",xaxis = axis_template,legend = list(x = 0,y = 100))
p1


p2 = plot_ly(x = set1$t1, opacity = 0.6, type = "histogram",name = names_ly1[1],marker = list(color = colors_ly[1])) %>%
  add_trace(x = set1$t2, opacity = 0.6, type = "histogram",name = names_ly1[2],marker = list(color = colors_ly[2])) %>%
  add_trace(x = set1$t3, opacity = 0.6, type = "histogram",name = names_ly1[3],marker = list(color = colors_ly[3])) %>%
  add_trace(x = set2$t1, opacity = 0.6, type = "histogram",name = names_ly1[1],marker = list(color = colors_ly[1]),showlegend = FALSE) %>%
  add_trace(x = set2$t2, opacity = 0.6, type = "histogram",name = names_ly1[2],marker = list(color = colors_ly[2]),showlegend = FALSE) %>%
  add_trace(x = set2$t3, opacity = 0.6, type = "histogram",name = names_ly1[3],marker = list(color = colors_ly[3]),showlegend = FALSE) %>%
  layout(title = 'Histogram of Times for Syndrome I and II',barmode="overlay",xaxis = axis_template,legend = list(x = 0,y = 100))
p2




png('p1.png',height = 800,width = 800)
p1
dev.off()
