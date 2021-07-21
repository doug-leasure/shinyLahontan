rm(list=ls())
gc()

setwd('C:/RESEARCH/2015 NASA Trout/wd/lahontan_v2.01')

source('./global.R')

nsim <- 100
forecastyear=2045

pops <- popnames

extinctTable <- data.frame(row.names=pops,
                           GMU=rep(NA, length(pops)),
                           PopulationName=pops, 
                           extinct=rep(NA, length(pops)), 
                           extinctup=rep(NA, length(pops)), 
                           extinctlow=rep(NA, length(pops)))

for(popname in row.names(extinctTable)) {
  extinctTable[popname, 'GMU'] <- metapop[popname,'GMU']
}

for (popname in pops){
  i <- which(popnames==popname)
  
  covs <- simcov(popname, nsim=nsim, forecastyear=forecastyear)
  simN <- simpop(popname, covs, forecastyear=forecastyear)
  sim <- simstat(popname, simN, forecastyear=forecastyear)
  
  extinctTable[popname, c('extinct','extinctlow','extinctup')] <- sim[c('extinct','extinctlow','extinctup'),ncol(sim)]
  print(paste(popname,': ', round(extinctTable[popname,'extinct']*100, 1), '% (', round(extinctTable[popname,'extinctlow']*100, 1), '-', round(extinctTable[popname,'extinctup']*100, 1), '%)', sep=''))
  rm(covs, simN, sim)
  gc()
}

write.csv(extinctTable, 'data/extinct.csv')





