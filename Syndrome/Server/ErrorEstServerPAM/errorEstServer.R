## Error Analysis

source('ClusterSetup_Serv_ERRORest.R')

res = lapply(1:2,function(x){
  varScen = 1
  errorFunc(ps = 1, sd = x)})

saveRDS(res,file = 'res1.RDA')