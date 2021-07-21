ploteffect <- function(popname, covs, sim, forecastyear=2045, ci=0.95, model=d){
  # forecastyear=2045; ci=0.95; model=d
  
  popnum <- which(popnames==popname)
  if(length(popnum)==0) popnum <- 0
  
  # Parameter estimates
  b0r <- model$b0r; b1r <- model$b1r; b2r <- model$b2r; 
  b0phi <- model$b0phi; b1phi <- model$b1phi
  
  if(elevation[popname]==1) b2phi <- model[,'b2phi[1]']
  if(elevation[popname]==2) b2phi <- model[,'b2phi[2]']
  
  # Plot function
  myplot <- function(covname, scalename, b0,b1,b2=0, scalepop=NA, ylab=NA, xlab=NA, mar=c(4.5, 4, 2, 1)){
    # covname='hflowlag'; scalename='qmax3'; b0=b0phi; b1=b1phi;b2=0
    # scalepop=NA; ylab=NA; xlab=NA; mar=c(4.5, 4, 2, 1)

    rng <- range(jags.dat[[covname]])
    x <- seq(rng[1], rng[2], length=100)
    
    # Calculate y-value
    p <- data.frame(matrix(NA, nrow=length(b0), ncol=length(x)))
    
    for (i in 1:length(x)) { p[,i] <- b0 + b1*x[i] + b2*x[i]^2 }
    
    if(covname == 'templag') { p <- p + b2r*mean(covs[,'hflowlag',]) }
    if(covname == 'hflowlag') { p <- p + b1r*mean(covs[,'templag',]) }
    if(covname %in% c('templag','hflowlag')) { p <- (exp(p) - 1) * 100 }
    
    if(covname == 'bkt') { p <- p + b2phi*mean(covs[,'ndvi',]) }
    if(covname == 'ndvi') { p <- p + b1phi*mean(covs[,'bkt',]) }
    if(covname %in% c('bkt','ndvi')) { 
      r <- b0r + b1r*mean(covs[,'templag',]) + b2r*mean(covs[,'hflowlag',])
      r[r<0] <- 0
      p[p>0] <- -1e-3
      
      p <- r / -p
    }
    
    pmed <- apply(p, 2, median)
    pup <- apply(p, 2, quantile, probs=c(1-(1-ci)/2))
    plow <- apply(p, 2, quantile, probs=c((1-ci)/2))
    ypad <- diff(range(plow, pup))*0.01
    
    
    # Plot
    par(mar=mar)
    
    plot(NA, xlim=c(min(x), max(x)), ylim=c(min(plow)-8*ypad, max(pup)+8*ypad), xlab=xlab, ylab=ylab, xaxt='n', cex.lab=1.5)
    
    xpoly <- c(x, rev(x)); ypoly <- c(pup, rev(plow))
    polygon(x=xpoly, y=ypoly, col=rgb(0.8,0.8,0.8, 0.3), border=NA)
    
    # Past data
    if(popnum>0){
      j <- jags.dat[[covname]][popnum, paste('X',1985:2015,sep='')]
      i <- which(xpoly >= min(j) & xpoly <= max(j))
      
      if(diff(range(j))==0) { i <- which.min(abs(xpoly-mean(j)))
      } 
      if(length(i) > 2) { polygon(x=xpoly[i], y=ypoly[i], col=rgb(0.6,0.6,0.6, 0.5), border=NA) }
      
      if(diff(range(j))==0 | length(i)==0){
        y0 <- plow[which.min(abs(x-mean(j)))]
        y1 <- pup[which.min(abs(x-mean(j)))]
        segments(x0=mean(j), y0=y0, y1=y1)
        points(x=mean(j), y=y0, pch=16, cex=1.5)
        points(x=mean(j), y=y1, pch=16, cex=1.5)
        text('Past', x=mean(j), y=y0-8*ypad, cex=1.5, pos=4)
      } else{
        segments(x0=min(j), x1=max(j), y0=min(ypoly[i])-2*ypad, y1=min(ypoly[i])-2*ypad)
        points(x=min(xpoly[i]), y=min(ypoly[i])-2*ypad, pch=16, cex=1.5)
        points(x=max(xpoly[i]), y=min(ypoly[i])-2*ypad, pch=16, cex=1.5)
        text('Past', x=mean(xpoly[i]), y=min(ypoly[i])-8*ypad, cex=1.5)
      }  
    }

    # Future data
    
    j <- covs[, covname, as.character((lastyear+1):forecastyear)]
    i <- which(xpoly >= min(j) & xpoly <= max(j))
    
    if(diff(range(j))==0 | length(i)==0) { i <- which.min(abs(xpoly-mean(j)))  } 
    if(length(i) > 2){ polygon(x=xpoly[i], y=ypoly[i], col=rgb(0.6,0.6,0.6, 0.5), border=NA) }
    
    if(diff(range(j))==0 | length(i)==0){
      y0 <- plow[which.min(abs(x-mean(j)))]
      y1 <- pup[which.min(abs(x-mean(j)))]
      segments(x0=mean(j), y0=y0, y1=y1)
      points(x=mean(j), y=y0, pch=16, cex=1.5)
      points(x=mean(j), y=y1, pch=16, cex=1.5)
      text('Future', x=mean(j), y=y1+8*ypad, cex=1.5, pos=4)
    } else{
      segments(x0=min(xpoly[i]), x1=max(xpoly[i]), y0=max(ypoly[i])+2*ypad, y1=max(ypoly[i])+2*ypad)
      points(x=min(xpoly[i]), y=max(ypoly[i])+2*ypad, pch=16, cex=1.5)
      points(x=max(xpoly[i]), y=max(ypoly[i])+2*ypad, pch=16, cex=1.5)
      text(as.character('Future'), x=mean(xpoly[i]), y=max(ypoly[i])+8*ypad, cex=1.5)
    }
    
    lines(x=x, y=pup, lty=2, col=rgb(0.4, 0.4, 0.4))
    lines(x=x, y=plow, lty=2, col=rgb(0.4, 0.4, 0.4))
    lines(x=x, y=pmed, lty=1, lwd=2, col=rgb(0.4, 0.4, 0.4))
    
    # X-axis
    at <- seq(round(min(x)), round(max(x)))
    if(is.na(scalepop)) { 
      labels <- round(at * scale.factors[scale.factors$cov==scalename, 'sd'] + scale.factors[scale.factors$cov==scalename, 'mu'])
    } else { 
      labels <- round(at * scale.factors[scale.factors$cov==scalename & scale.factors$PopulationName==scalepop, 'sd'] + scale.factors[scale.factors$cov==scalename & scale.factors$PopulationName==scalepop, 'mu'], 2)
      # if(covname=='hflowlag') labels <- round(exp(labels) - 1)
      if(covname=='hflowlag') labels <- round(labels)
    }
    axis(1, at=at, labels=labels, cex=1.5)
  }
  
  layout(matrix(c(1,2,3,4), nrow=2, ncol=2, byrow=T), heights=c(1, 1), widths=c(1, 1))
  par(oma=c(0,0,6,0))
  
  # Temperature
  myplot(covname='templag', scalename='temp', b0=b0r, b1=b1r, 
         xlab='Water Temperature (C)', ylab='Growth Potential (%)', mar=c(5, 5, 2, 1))
  
  # High flow
  myplot(covname='hflowlag', scalename='qmax3', scalepop=popname, b0=b0r, b1=b2r,
         xlab='High Flow (m3/h)', ylab=NA, mar=c(5, 4, 2, 2))
  
  # BKT
  myplot(covname='bkt', scalename='bkt', b0=b0phi, b1=b1phi, 
         xlab='Non-native Trout (fish / km)', ylab='LCT Carrying Capacity (fish / km)', mar=c(5, 5, 1, 1))
  
  # NDVI
  myplot(covname='ndvi', scalename='pndvi', scalepop=popname, b0=b0phi, b1=b2phi, 
         xlab='Riparian Vegetation (NDVI)', ylab=NA, mar=c(5, 4, 1, 2))
  
  # Title
  mtext(paste('Population: ',popname,
              '\nExtinction Risk: ', round(sim['extinct', ncol(sim)] * 100, 1),'% (', round(sim['extinctlow', ncol(sim)] * 100, 1),' - ',round(sim['extinctup', ncol(sim)] * 100, 1),'%)',
              sep=''), outer=TRUE, cex=1.5)
  
}
