plotcov <- function(popname, covs, sim, forecastyear=2045, ci=0.95){
  
  # forecastyear=2045; ci=0.95;
  
  popnum <- which(popnames==popname)
  if(length(popnum)==0) popnum <- 0
  
  if(popnum==0){
    year1.sim <- lastyear
    years <- lastyear:forecastyear
  } else {
    attach(jags.dat[c('t_ti','nt','templag','hflowlag','bkt','ndvi','burn')], warn.conflicts=F)
    year1.sim <- firstyear + t_ti[popnum,nt[popnum]] -1
    years <- (firstyear + t_ti[popnum,1] - 1):forecastyear
  }
  
  # Plot Function
  myplot <- function(covname, ylab, scalename, lag=F, scalepop=NA, main=NA, xaxis=F, mar=c(1.5, 7, 1.5, 1), plot1=F){
    # covname='bkt'; ylab='Non-native Density\n(fish per km)'; scalename='bkt';
    # lag=F; scalepop=NA; main=NA; xaxis=F
    # mar=c(1.5, 7, 1.5, 1); plot1=T
    
    
    # Setup data
    dat <- data.frame(row.names=c('low','mid','up'), matrix(NA, nrow=3, ncol=length(years)))
    names(dat) <- years
    
    # Pre 2015
    if(popnum>0) dat['mid', as.character(years[1]:lastyear)] <- dat['up', as.character(years[1]:lastyear)] <- dat['low', as.character(years[1]:lastyear)] <- jags.dat[[covname]][popnum,paste('X',years[1]:lastyear,sep='')]
    
    # Post 2015
    for (yr in as.character((lastyear+1):forecastyear)){
      dat[c('low','mid','up'), yr] <- quantile(covs[,covname,yr], probs=c((1-ci)/2, 0.5, 1-(1-ci)/2))
      dat['mid', yr] <- covs[1,covname,yr]
    }
    
    # Rescale
    if(is.na(scalepop)){
      mu <- scale.factors[scale.factors$cov==scalename,'mu']
      sd <- scale.factors[scale.factors$cov==scalename,'sd']
    } else {
      mu <- scale.factors[scale.factors$cov==scalename & scale.factors$PopulationName==popname,'mu']
      sd <- scale.factors[scale.factors$cov==scalename & scale.factors$PopulationName==popname,'sd']
    }
    
    dat <- dat * sd + mu
    
    # if(covname=='hflowlag') dat <- exp(dat) - 1
    
    # Axis limits
    xlim <- range(as.numeric(names(dat)))
    
    ylim <- range(dat, na.rm=T)
    
    if (covname=='bkt'){
      ylim <- round(ylim)
      ylim[1] <- min(ylim[1], as.matrix(nonnatives[popname, as.character(years[years<=lastyear])]), na.rm=T)
      ylim[2] <- max(ylim[2], as.matrix(nonnatives[popname, as.character(years[years<=lastyear])]), na.rm=T)
    }
    
    if(diff(ylim)==0) {
      ylim[1] <- ylim[1] * 0.9
      ylim[2] <- ylim[2] * 1.1
    }
    
    if(plot1) ylim[2] <- ylim[2] * 1.2
    
    if (diff(ylim)==0) ylim <- c(0,1) 

    #Plot
    par(mar=mar)
    
    if(xaxis){
      plot(NA, xlim=xlim, ylim=ylim, xlab='Year', xaxt='n', ylab=ylab, main=main, cex.main=2, cex.lab=2, cex.axis=2)
      xpad <- diff(range(xlim))*0.04; ypad <- diff(range(ylim))*0.04
      rect(xleft=lastyear+0.5, xright=xlim[2]+xpad, ybottom=ylim[1]-ypad, ytop=ylim[2]+ypad*0.80, col=rgb(0.9,0.9,0.9,0.5), border=NA)
      axis(side=1, at=seq(1990, max(xlim), 10), lwd.ticks=2, cex.axis=2)
      axis(side=1, at=seq(1985, max(xlim), 10), labels=NA)
    } else {
      plot(NA, xlim=xlim, ylim=ylim, xlab=NA, xaxt='n', ylab=ylab, main=main, cex.main=2, cex.lab=2, cex.axis=2)
      xpad <- diff(range(xlim))*0.04; ypad <- diff(range(ylim))*0.04
      rect(xleft=lastyear+0.5, xright=xlim[2]+xpad, ybottom=ylim[1]-ypad, ytop=ylim[2]+ypad*0.80, col=rgb(0.9,0.9,0.9,0.5), border=NA)
      axis(side=1, at=seq(1990, max(xlim), 10), labels=NA, lwd.ticks=3)
      axis(side=1, at=seq(1985, max(xlim), 10), labels=NA)
    }

    if(lag) {
      dat <- data.frame(cbind(dat[,2:ncol(dat)], matrix(NA, ncol=1, nrow=3)))
      names(dat) <- years
      polygon(x=c(years[-length(years)], rev(years[-length(years)])), 
              y=c(dat['up',as.character(years[-length(years)])], rev(dat['low',as.character(years[-length(years)])])), 
              col=rgb(0.8,0.8,0.8, 0.3), border=NA)
    } else{
      polygon(x=c(years, rev(years)), 
              y=c(dat['up',as.character(years)], rev(dat['low',as.character(years)])), 
              col=rgb(0.8,0.8,0.8, 0.3), border=NA)
    }
    
    lines(x=years, y=dat['up', as.character(years)], type='l', lty=2, lwd=1, col=rgb(0.4, 0.4, 0.4))
    lines(x=years, y=dat['low', as.character(years)], type='l', lty=2, lwd=1, col=rgb(0.4, 0.4, 0.4))
    lines(x=years, y=dat['mid', as.character(years)], type='l', lty=1, lwd=2, col=rgb(0.4, 0.4, 0.4))
    
    if(plot1) {
      text(x=lastyear+1, y=ylim[2]-ypad*2, labels='FORECASTS >>', pos=4, cex=2, font=4, col=rgb(0.6, 0.6, 0.6))
    }
    if(covname=='bkt'){
      if(popname %in% row.names(nonnatives)) {
        y <- as.matrix(nonnatives[popname, as.character(years[years<=lastyear])])
        y <- c(y, rep(NA, length(years) - length(years[years<=lastyear])))
        points(x=years, y=y, pch=16, col=rgb(0.4, 0.4, 0.4), cex=1.75)
      }
    }
  }
  
  
  # Setup Plot
  layout(matrix(c(1,2,3,4), nrow=4, ncol=1), heights=c(1.05, 1, 1, 1.2))
  par(oma=c(0,0,6,0))
  
  # Non-natives
  myplot(covname='bkt', ylab='Non-native Density\n(fish per km)', scalename='bkt', mar=c(1.5, 7, 1.5, 1), plot1=T)

  # Temperature
  myplot(covname='templag', ylab='Water Temp\n(C)', scalename='temp', lag=T)
  
  # High flow
  myplot(covname='hflowlag', ylab='High Flow\n(m3/h)', scalename='qmax3', scalepop=popname, lag=T)
  
  # NDVI
  myplot(covname='ndvi', ylab='Riparian Vegetation\n(NDVI)', scalename='pndvi', scalepop=popname, mar=c(4.5, 7, 1.5, 1), xaxis=T)
  
  # Title
  mtext(paste('Population: ',popname,
              '\nExtinction Risk: ', round(sim['extinct', ncol(sim)] * 100, 1),'% (', round(sim['extinctlow', ncol(sim)] * 100, 1),' - ',round(sim['extinctup', ncol(sim)] * 100, 1),'%)',
              sep=''), 
        outer=TRUE, cex=1.5)
}