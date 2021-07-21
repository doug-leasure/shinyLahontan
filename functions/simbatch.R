simbatch <- function(primepop, pops, forecastyear=2045, nsim=10, nmcmc=NA, natvar=FALSE, futrscale=1, hist.bkt=FALSE, extent=0,
                     cut=list(templag=c(0,1), hflowlag=c(0,1), bkt=c(0,1), ndvi=c(0,1)), 
                     const=list(reintro=0, templag=NA, hflowlag=NA, bkt=NA, ndvi=NA),
                     demsto=NA, ci=0.95, Nmax=1e9, model=d,
                     add1=0, add1t=1, add2=0, add2t=1){
  # pops=allpopnames[1:5];
  # forecastyear=2045; nsim=10; natvar=FALSE; hist.bkt=FALSE; extent=NA;
  # cut=list(templag=c(0,1), hflowlag=c(0,1), bkt=c(0,1), ndvi=c(0,1));
  # const=list(reintro=0, templag=NA, hflowlag=NA, bkt=NA, ndvi=NA);
  # demsto=NA; ci=0.95; Nmax=1e9; model=d
  
  bkt.setting <- const[['bkt']]
  if(is.na(bkt.setting) | round(const[['bkt']]) == round(bkt.dat[primepop, as.character(lastyear)]*scale.factors[scale.factors$cov=='bkt','sd'])){
    bkt.setting <- 'default'
  }
  
  if(length(pops)==0){
    output <- data.frame(NoSelection=c('You have not selected any populations','Scroll to the "Batch Processing" settings','Dont forget to click "Run Simulations..."')) 
    return(output)
  }else {
    output <- data.frame(row.names=pops, PopulationName=pops, extinct=rep(NA,length(pops)), extinctup=rep(NA,length(pops)), extinctlow=rep(NA,length(pops)))
    
    probs <- c((1-ci)/2, 0.5, 1-(1-ci)/2)
    
    for(popname in pops){
      
      if(bkt.setting=='default'){
        const[['bkt']] <- bkt.dat[popname, as.character(lastyear)]*scale.factors[scale.factors$cov=='bkt','sd']
      }
      
      simN <- simpop(popname=popname, 
                     covs=simcov(popname=popname, forecastyear=forecastyear, nsim=nsim, natvar=natvar, futrscale=futrscale, hist.bkt=hist.bkt, extent=extent, cut=cut, const=const), 
                     demsto=demsto, forecastyear=forecastyear, ci=ci, model=model)
      
      if (nsim == 1) {
        output[popname, c('extinctlow','extinct','extinctup')] <- paste(round(mean(simN$N[1,,nyears]==0)*100, 1), '%')
      } else {
        output[popname, c('extinctlow','extinct','extinctup')] <- paste(round(quantile(apply(simN$N[,,dim(simN$N)[3]]==0, 1, mean), probs=probs)*100, 1), '%')
      }
    }
    row.names(output) <- NULL #metapop[pops,'PopID']
    names(output) <- c('Population','Extinction Risk','Upper','Lower')
    output <- output[,c('Population','Extinction Risk','Lower','Upper')]
    
    for(col in names(output)) output[,col] <- as.character(output[,col])
    
    output <- rbind(output, rep(NA,4))
    output <- rbind(output, rep(NA,4))
    output <- rbind(output, rep(NA,4))
    output <- rbind(output, rep(NA,4))
    output <- rbind(output, rep(NA,4))
    
    output <- rbind(output, c('SIMULATION SETTING', 'VALUE', NA, NA))
    output <- rbind(output, c('Forecast Year', forecastyear, NA, NA))
    if(extent >= 0) output <- rbind(output, c('Change in extent', paste('+',extent,'km'), NA, NA))
    if(extent < 0) output <- rbind(output, c('Change in extent', paste('-',extent,'km'), NA, NA))
    if(hist.bkt) { output <- rbind(output, c('Non-native density', 'historical', NA, NA))
    } else { output <- rbind(output, c('Non-native density', bkt.setting, NA, NA)) }
    output <- rbind(output, c('Environmental stochasticity', 'default', NA, NA))
    output <- rbind(output, c('Reintro 1 count', add1, NA, NA))
    output <- rbind(output, c('Reintro 1 years', add1t, NA, NA))
    output <- rbind(output, c('Reintro 2 count', add2, NA, NA))
    output <- rbind(output, c('Reintro 2 years', add2t, NA, NA))
    output <- rbind(output, c('Temperature min', cut[['templag']][1], NA, NA))
    output <- rbind(output, c('Temperature max', cut[['templag']][2], NA, NA))
    output <- rbind(output, c('Flow min', cut[['hflowlag']][1], NA, NA))
    output <- rbind(output, c('Flow max', cut[['hflowlag']][2], NA, NA))
    output <- rbind(output, c('NDVI min', cut[['ndvi']][1], NA, NA))
    output <- rbind(output, c('NDVI max', cut[['ndvi']][2], NA, NA))
    output <- rbind(output, c('Credible intervals', ci, NA, NA))
    output <- rbind(output, c('Number of simulations', nsim, NA, NA))
    output <- rbind(output, c('Number of MCMC samples', nmcmc, NA, NA))
    
    return(output)
  }
}