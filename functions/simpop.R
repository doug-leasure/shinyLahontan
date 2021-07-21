simpop <- function(popname, covs, demsto=NA, forecastyear=2045, 
                   ci=0.95,
                   model=d){
  
  # demsto=NA; forecastyear=2045;
  # ci=0.95; 
  # model=d

  # Identify number of simulations from covariates object
  nsim <- dim(covs)[1]
  
  # Identify number of MCMC iterations from JAGS model
  mcmc <- nrow(model)
  
  # Identify population number that corresponds to JAGS indexing
  popnum <- which(popnames==popname)
  if(length(popnum)==0) popnum <- 0
  
  # If reconnect habitat, identify sub-population names
  if(popname %in% recnames) {
    subpopnames <- metapop[metapop$MetaPop_Name==metapop[popname, 'MetaPop_Name'], 'PopulationName']
    subpopnames <- subpopnames[!grepl('Reconnect', subpopnames) & !is.na(subpopnames)]
    subpopnums <- which(popnames %in% subpopnames)
  } else {subpopnums <- c()}
  
  # "Attach" objects from jags.dat
  t_ti <- jags.dat$t_ti
  nt <- jags.dat$nt
  
  # Get demographic stochasticity (sigmaR) from model results
  if(popnum > 0){
    # Populations with data
    sigmaR <- model[,paste('sigmaR[',popnum,']',sep='')]  
  } else {
    # Populations without data
    if(popname %in% recnames & length(subpopnums)>0) {
      # Reconnect habitat
      sigmaR <- apply(as.matrix(model[,paste('sigmaR[',subpopnums,']',sep='')]), 1, mean)
    } else {
      # Populations without data (Normal-Gamma mixture to represent the half-cauchy; from Martin Plummer JAGS manual)
      sigmaR <- abs( rnorm(n=mcmc, mean=model[,'musig'], sd=model[,'sigsig']) / sqrt(rgamma(n=mcmc, shape=1/2, rate=1/2)) )
    }
  }
  
  # Rescale demographic stochasticity based on user input
  if(!is.na(demsto)) sigmaR <- sigmaR + ( demsto - median(sigmaR) )
  sigmaR[sigmaR < 0] <- 0
    
  # Get regression coefficients from JAGS model
  b0r <- model$b0r; b1r <- model$b1r; b2r <- model$b2r
  b0phi <- model$b0phi; b1phi <- model$b1phi 
  if(elevation[popname]==1) b2phi <- model[,'b2phi[1]']
  if(elevation[popname]==2) b2phi <- model[,'b2phi[2]']
  
  # Define first year of simulation (2015, or last year with data)
  if(popnum==0) { 
    year1.sim <- lastyear
  } else { 
    year1.sim <- firstyear + t_ti[popnum,nt[popnum]] - 1 
  }
  
  # Calculate total years in simulation run
  nyears <- forecastyear - year1.sim + 1
  
  # Setup N
  N <- array(NA, dim=c(nsim, mcmc, nyears))
  
  # Initial N for simulation
  if(popnum==0){ 
    # If no data, initialize at 0
    N[,,1] <- 0
    # If reconnect habitat, sum of all sub-pops N[lastyear]
    if(popname %in% recnames) {
      for (subpopnum in subpopnums){
        N[,,1] <- N[,,1] + round(median(model[,paste('N[',subpopnum,',',jags.dat$t_ti[subpopnum,jags.dat$nt[subpopnum]],']',sep='')]))
      }
    }
  } else { 
    # If pop with data, last year N[lastyear]
    N[,,1] <- model[,paste('N[',popnum,',',t_ti[popnum,nt[popnum]],']',sep='')] 
  }
  
  # Setup arrays to monitor simulation parameters
  k <- phi  <- r <- array(NA, dim=c(nsim, mcmc, nyears))

  # Run simulations
  for (s in 1:nsim){
    
    # For each year after last year with data
    for (t in 2:nyears){
      
      # Adjust N[t-1] based on translocation data
      N[s,,t-1] <-  N[s,,t-1] + covs[s,'reintro',t-1]
      
      # Ensure no N[t-1] < 0
      N[s,,t-1][N[s,,t-1] < 0] <-  0
      
      # Intrinsic growth rate
      r[s,,t] <- b0r + b1r*covs[s,'templag',t] + b2r*covs[s,'hflowlag',t]

      # Density-dependence
      phi[s,,t] <- b0phi + b1phi*covs[s,'bkt',t] + b2phi*covs[s,'ndvi',t]
      phi[s,,t] <- sapply(phi[s,,t], min, -1e-6)

      # Carrying capacity
      k[s,,t] <- sapply(r[s,,t] / -phi[s,,t], max, 0) # * extent
      
      # Realized growth rate
      Rmean <- r[s,,t] + phi[s,,t] * N[s,,t-1] / max(0.001, covs[s,'extent',t])
      
      # Stochastic realized growth rates
      R <- rnorm(mcmc, Rmean, sigmaR)
      
      # Prevent realized growth rate (R) from exceeding intrinsic growth rate (r)
      cap <- max(r[s,,t])
      needCaps <- sum(R > cap) > 0
      while(needCaps){
        x <- R > cap
        R[x] <- rnorm(sum(x), Rmean[x], sigmaR[x])
        needCaps <- sum(R > cap) > 0
      }
      
      # Expected value of N
      muN <- N[s,,t-1] * exp(R)
      
      # Stochastic N
      N[s,,t] <- rpois(mcmc, muN)
    }
  }
  
  return(list(N=N, r=r, phi=phi, k=k, sigmaR=sigmaR))
}


