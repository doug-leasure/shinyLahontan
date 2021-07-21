simstat <- function(popname, simN, forecastyear=2045, ci=0.95, model=d){
  # popname='Abel'; simN=simN; forecastyear=2045; ci=0.95; model=d
  
  # Get popnum (corresponds with JAGS indexing)
  popnum <- which(popnames==popname)
  
  # popnum=0 means the population has no data
  if(length(popnum)==0) popnum <- 0
  
  # Determine the first year of the simulation (the last year with data, or 2015)
  if(popnum==0) { 
    year1.sim <- lastyear
  } else { 
    year1.sim <- firstyear + jags.dat$t_ti[popnum,jags.dat$nt[popnum]] - 1 
  }
  
  # Populations that had only data from 2015 can be treated as popnum=0
  if(year1.sim == lastyear) popnum <- 0
  
  # Number of simulation years
  nyears <- forecastyear - year1.sim + 1
  
  # Setup output table (from last year with data to end of forecast)
  rn <- c('Nmean','Nup','Nlow','rmean','rup','rlow','kmean','kup','klow','phimean','phiup','philow','stomean','stoup','stolow','extinct')
  out.df <- data.frame(row.names=rn, matrix(NA, nrow=length(rn), ncol=nyears))
  names(out.df) <-  as.character(year1.sim:forecastyear) #names(templag)
  
  # Calculate summary statistics for simulation results and write them to the table
  probs <- c((1-ci)/2, 0.5, 1-(1-ci)/2)
  
  out.df[c('Nlow','Nmean','Nup'),1] <- quantile(as.vector(simN$N[,,1]), probs=probs)

  for(t in 2:nyears){
    out.df[c('Nlow','Nmean','Nup'),t] <- quantile(as.vector(simN$N[,,t]), probs=probs)
    out.df[c('rlow','rmean','rup'),t] <- quantile(as.vector(simN$r[,,t]), probs=probs)
    out.df[c('klow','kmean','kup'),t] <- quantile(as.vector(simN$k[,,t]), probs=probs)
    out.df[c('philow','phimean','phiup'),t] <- quantile(as.vector(simN$phi[,,t]), probs=probs)

    # out.df['Nmean',t] <- mean(as.vector(simN$N[1,,t]))
    # out.df['rmean',t] <- mean(as.vector(simN$r[1,,t]))
    # out.df['kmean',t] <- mean(as.vector(simN$k[1,,t]))
    # out.df['phimean',t] <- mean(as.vector(simN$phi[1,,t]))
  }
  
  # Calculate extinction probability based on end states of population simulations
  if (dim(simN$N)[1] == 1) {
    out.df[c('extinctlow','extinct','extinctup'),nyears] <- mean(simN$N[1,,nyears]==0)
  } else {
    x <- apply(simN$N[,,nyears]==0, 1, mean)
    out.df['extinctlow',nyears] <- min(x)
    out.df['extinct',nyears] <- median(x)
    out.df['extinctup',nyears] <- max(x)
    rm(x)
  }

  
  #### Add on model-based results between the first and last year with field data ####
  
  if(popnum>0){
    
    year1.wdat <- (firstyear-1) + jags.dat$t_ti[popnum,1] 
    
    nyears.presim <- year1.sim - year1.wdat
    
    # Setup data frame
    out.df <- cbind(data.frame(matrix(NA, nrow=nrow(out.df), ncol=nyears.presim)), out.df)
    names(out.df) <- year1.wdat:forecastyear
    
    out.df[c('stolow','stomean','stoup'),] <- quantile(as.vector(as.matrix(simN$sigmaR)), probs=probs)
    
    if(elevation[popname]==1) b2phi <- model[,'b2phi[1]']
    if(elevation[popname]==2) b2phi <- model[,'b2phi[2]']
    
    for (t in 1:(nyears.presim+1)){
      year <- year1.wdat + t - 1
      t_N <- year - firstyear + 1
      
      r <- model[,'b0r'] + model[,'b1r'] * jags.dat$templag[popnum,t_N] + model[,'b2r'] * jags.dat$hflowlag[popnum,t_N]
      
      phi <- model[,'b0phi'] + model[,'b1phi'] * jags.dat$bkt[popnum,t_N] + b2phi * jags.dat$ndvi[popnum,t_N]
      phi <- sapply(phi, min, -1e-6)
      
      k <- sapply(r/-phi, max, 0)# * covs[[1]]$extent

      out.df[c('Nlow','Nmean','Nup'),t] <- quantile(model[,paste('N[',popnum,',',t_N,']',sep='')], probs=probs)
      out.df[c('rlow','rmean','rup'),t] <- quantile(r, probs=probs)
      out.df[c('klow','kmean','kup'),t] <- quantile(k, probs=probs)
      out.df[c('philow','phimean','phiup'),t] <- quantile(phi, probs=probs)
    }
  }
  return(out.df)
}