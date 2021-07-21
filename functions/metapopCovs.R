rm(list=ls())

setwd('C:/RESEARCH/2015 NASA Trout/wd/lahontan')

source('./global.R')
scale.factors[scale.factors$cov=='bkt', 'sd'] <- scale.factors[scale.factors$cov=='bkt', 'sd'] / 1000


wtavg <- function(x, weights){
  sum(x * (weights/sum(weights)))
}

for (popname in streams.reconnect$PopulationName){
  metapopname <- metapop[popname, 'MetaPop_Name']
  if(!metapopname==''){
    subpopnames <- metapop[metapop$MetaPop_Name==metapopname, 'PopulationName']
    subpopnames <- subpopnames[!grepl('Reconnect', subpopnames)]
    
    # Extent
    extent.dat[popname,'extent'] <- sum(extent.dat[subpopnames,'extent'])
    
    # Temp
    templag.dat[popname,2:ncol(templag.dat)] <- as.data.frame(t(apply(templag.dat[subpopnames, 2:ncol(templag.dat)], 2, wtavg, weights=extent.dat[subpopnames,'extent'])))
    
    # Non-native
    bkt.dat[popname,2:ncol(bkt.dat)] <- as.data.frame(t(apply(bkt.dat[subpopnames, 2:ncol(bkt.dat)], 2, max)))
    
    # NDVI
    dtemp <- ndvi.dat[subpopnames, 2:ncol(ndvi.dat)]
    for(subpopname in subpopnames){
      mu <- scale.factors[scale.factors$cov=='pndvi' & scale.factors$PopulationName==subpopname,'mu']
      sigma <- scale.factors[scale.factors$cov=='pndvi' & scale.factors$PopulationName==subpopname,'sd']
      dtemp[subpopname,] <- dtemp[subpopname,] * sigma + mu
    }
    dtemp <- apply(dtemp, 2, wtavg, weights=extent.dat[subpopnames,'extent'])
    
    scale.factors[scale.factors$PopulationName==popname & scale.factors$cov=='pndvi','mu'] <- mean(dtemp)
    scale.factors[scale.factors$PopulationName==popname & scale.factors$cov=='pndvi','sd'] <- sd(dtemp)
    
    ndvi.dat[popname,2:ncol(ndvi.dat)] <- (dtemp - mean(dtemp))/sd(dtemp)
    
    # High flow
    dtemp <- hflowlag.dat[subpopnames, 2:ncol(hflowlag.dat)]
    for(subpopname in subpopnames){
      mu <- scale.factors[scale.factors$cov=='qmax3' & scale.factors$PopulationName==subpopname,'mu']
      sigma <- scale.factors[scale.factors$cov=='qmax3' & scale.factors$PopulationName==subpopname,'sd']
      dtemp[subpopname,] <- dtemp[subpopname,] * sigma + mu
    }
    dtemp <- exp(dtemp)
    dtemp <- apply(dtemp, 2, wtavg, weights=extent.dat[subpopnames,'extent'])
    dtemp <- log(dtemp)
    
    scale.factors[scale.factors$PopulationName==popname & scale.factors$cov=='qmax3','mu'] <- mean(dtemp)
    scale.factors[scale.factors$PopulationName==popname & scale.factors$cov=='qmax3','sd'] <- sd(dtemp)
    
    hflowlag.dat[popname,2:ncol(hflowlag.dat)] <- (dtemp - mean(dtemp))/sd(dtemp)  
  }
}

write.csv(extent.dat, file='./data/trimcovs/extent.dat.csv', row.names=F)
write.csv(bkt.dat, file='./data/trimcovs/bkt.dat.csv', row.names=F)
write.csv(templag.dat, file='./data/trimcovs/templag.dat.csv', row.names=F)
write.csv(ndvi.dat, file='./data/trimcovs/ndvi.dat.csv', row.names=F)
write.csv(hflowlag.dat, file='./data/trimcovs/hflowlag.dat.csv', row.names=F)
write.csv(scale.factors, file='./data/scale.factors.csv', row.names=F)


