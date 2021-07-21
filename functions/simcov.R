simcov <- function(popname, forecastyear=2045, nsim=10, natvar=F, futrscale=1, hist.bkt=FALSE, extent=0,
                    cut=list(templag=c(0,1), hflowlag=c(0,1), bkt=c(0,1), ndvi=c(0,1)), 
                    const=list(reintro=0, templag=NA, hflowlag=NA, bkt=NA, ndvi=NA)
                    ){
  # popname='Abel'; forecastyear=2045; nsim=10; natvar=FALSE; futrscale=1; hist.bkt=TRUE; extent=NA;
  # cut=list(templag=c(0,1), hflowlag=c(0,1), bkt=c(0,1), ndvi=c(0,1))
  # const=list(reintro=0, templag=NA, hflowlag=NA, bkt=NA, ndvi=NA)

  # "Attach" objects from jags.dat
  t_ti <- jags.dat$t_ti
  nt <- jags.dat$nt
  
  # Identify popnum that corresponds to JAGS indexing
  popnum <- which(popnames==popname)
  if (length(popnum)==0) popnum <- 0
  
  # If non-native density is not provided, set it to last observed density.
  if(is.na(const['bkt'])) const['bkt'] <- bkt.dat[popname, ncol(bkt.dat)] * scale.factors[scale.factors$cov=='bkt','sd']
    
  # If resampling historic non-native densities, set other non-native setting to NA
  if(hist.bkt) const['bkt'] <- NA
  
  # Build list with each covariate time series from 1984-2015
  covs <- list(extent = extent.dat[popname, 'extent'], 
               templag = templag.dat[popname, 2:ncol(templag.dat)], 
               hflowlag = hflowlag.dat[popname, 2:ncol(hflowlag.dat)], 
               bkt = bkt.dat[popname, 2:ncol(bkt.dat)], 
               ndvi = ndvi.dat[popname, 2:ncol(ndvi.dat)],
               reintro = add.dat[popname, 2:ncol(add.dat)] - rem.dat[popname, 2:ncol(rem.dat)]
               )  

  # Get covariate names
  covnames <- names(covs)[which(!names(covs)=='extent')]
  
  # Define first year of simulation
  if(popnum==0) { 
    year1.sim <- lastyear
  } else { 
    year1.sim <- (firstyear - 1) + t_ti[popnum,nt[popnum]] 
  }
  
  # Identify years with data for all covariates
  resamp.years <- names(covs[[covnames[1]]])
  for(cov in covnames) resamp.years <- intersect(resamp.years, names(covs[[cov]]))
  resamp.years <- sort(resamp.years)
  
  # Identify years when covariates are within constraints of cut[]
  for(cov in names(cut)){
    q <- quantile(as.numeric(covs[[cov]]), probs=cut[[cov]])
    resamp.years <- resamp.years[ resamp.years %in% names(covs[[cov]])[as.numeric(covs[[cov]]) >= q[1] & as.numeric(covs[[cov]]) <= q[2]] ]
  }
  
  # Setup covariate time series for each simulation
  covs.sims <- array(dim=c(nsim, length(covs), forecastyear-year1.sim+1))
  dimnames(covs.sims) <- list(c(), names(covs), as.character(year1.sim:forecastyear))
  
  # Setup extent time series (a static time series at this point)
  if(is.na(extent)){
    covs.sims[, 'extent',] <- covs[['extent']]
  } else{
    i2 <- which(dimnames(covs.sims)[[3]]==lastyear)
    covs.sims[, 'extent',1:i2] <- covs[['extent']]
    
    i1 <- which(dimnames(covs.sims)[[3]]==lastyear+1)
    i2 <- which(dimnames(covs.sims)[[3]]==as.character(forecastyear))
    covs.sims[, 'extent',i1:i2] <- max(0, covs[['extent']] + extent)
  }
  
  for (s in 1:nsim){
    
    # Sample with replacement from resamp.years to create a vector of forecast years
    years.resamp <- sample(x=resamp.years, size=forecastyear-lastyear, replace=T)
    
    for (cov in covnames) {
      
      # Entire covariate time series (1984-2015)
      x <- covs[[cov]]
      
      # Covariates from last year with LCT data to 2015
      x <- x[,which(names(x)==year1.sim):ncol(x)]
      
      # Setup future covariate time series
      if (cov %in% names(const)[which(!is.na(const))]) { 
        
        # Add constant covariate values for future
        xfuture <- matrix(rep(const[[cov]],length=length(years.resamp)), nrow=1, ncol=length(years.resamp))
        if(cov=='bkt') xfuture <- xfuture / scale.factors[scale.factors$cov=='bkt','sd']
        x <- cbind(x, xfuture)
      
      } else { 
        
        # Add resampled covariate values for future
        x <- cbind(x, covs[[cov]][years.resamp]) 
      }
      
      # Name resampled years consecutively
      names(x) <- as.character(year1.sim:forecastyear)
      
      # Rescale to maintain natural variance
      if(natvar==T & cov %in% c('templag','hflowlag','ndvi')) {
        xall <- as.matrix(covs[[cov]])
        xcut <- as.numeric(x[as.character((lastyear+1):forecastyear)])
        x[which(names(x)==lastyear+1):length(x)] <- (x[which(names(x)==lastyear+1):length(x)] - mean(xcut)) / sd(xcut) * (sd(xall) * futrscale) + mean(xcut) 
        
        # Ensure that hflow isn't negative
        if(cov=='hflowlag'){ 
          mu <- scale.factors[scale.factors$cov=='qmax3' & scale.factors$PopulationName==popname,'mu']
          std <- scale.factors[scale.factors$cov=='qmax3' & scale.factors$PopulationName==popname,'sd']
          zeroflow <- (0-mu)/std
          x[x<zeroflow] <- zeroflow
        }
      }
      covs.sims[s, cov, ] <- as.vector(as.matrix(x))
    }
  }
  return(covs.sims)
}