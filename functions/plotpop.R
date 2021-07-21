plotpop <- function(popname, dat, Nmax=1e9){
  
  popnum <- which(popnames==popname)
  if(length(popnum)==0) popnum <- 0
  
  if (popnum > 0)  {
    t_ti <- jags.dat$t_ti
    nt <- jags.dat$nt
    Yi <- jags.dat$Yi
  }
  
  # Data
  yrs <- 1:ncol(dat)
  x.yrs <- as.numeric(names(dat))

  # Setup Plot
  layout(matrix(c(1,2,3,4), nrow=4, ncol=1), heights=c(0.20, 0.3, 0.25, 0.30))
  par(oma=c(0,0,6,0))
  
  # Top panel: Survey data
  par(mar=c(1.5,7,1.5,1))
  
  if(popnum==0){
    ytop <- 0
    xmid <- 0
  } else {
    ytop <- c()
    xmid <- c()
    for(ti in 1:nt[popnum]){
      ytop <- c(ytop, Yi[popnum,ti])
      xmid <- c(xmid, (firstyear-1)+t_ti[popnum,ti])
    }
  }
  
  xlim <- range(as.numeric(x.yrs))
  ylim <- c(0, max(5, max(ytop)))
  ylim[2] <- ylim[2] + diff(ylim)*0.04

  plot(NA, xlab=NA, xaxt='n', xlim=xlim, ylim=ylim,
       cex.lab=2, cex.axis=2,
       ylab='Fish Observed\n(see Data tab)',
       main=NA
  )
  if(popnum>0){
    for(ti in 1:nt[popnum]){
      rect(xleft=xmid[ti]-0.5, xright=xmid[ti]+0.5, ybottom=0, ytop=ytop[ti], 
           col=rgb(0.8,0.8,0.8), border=rgb(0.4,0.4,0.4), lwd=2)
    }  
  }
  mtext(paste('Population: ', popname, 
              '\nExtinction Risk: ', round(dat['extinct', ncol(dat)] * 100, 1),'% (', round(dat['extinctlow', ncol(dat)] * 100, 1),' - ',round(dat['extinctup', ncol(dat)] * 100, 1),'%)',
              sep=''), 
        cex=1.5, outer=TRUE)
  axis(side=1, at=seq(1990, max(x.yrs), 10), labels=NA, lwd.ticks=3)
  axis(side=1, at=seq(1985, max(x.yrs), 10), labels=NA)
  
  
  # Population Size
  par(mar=c(1.5,7,1.5,1))

  xlim=range(as.numeric(x.yrs))
  ylim=c(0, 1.20*min(Nmax, max(10, max(dat['Nup',yrs], na.rm=T))))

  plot(NA, xlab=NA, xaxt='n', xlim=xlim, ylim=ylim,
       cex.lab=2, cex.axis=2,
       ylab='Total\nPopulation Size'
  )
  xpad <- diff(range(xlim))*0.04; ypad <- diff(range(ylim))*0.04
  rect(xleft=lastyear+0.5, xright=xlim[2]+xpad, ybottom=ylim[1]-ypad, ytop=ylim[2]+ypad*0.89, col=rgb(0.9,0.9,0.9,0.5), border=NA)
  
  polygon(x=c(x.yrs, rev(x.yrs)), y=c(dat['Nup',yrs], rev(dat['Nlow',yrs])), col=rgb(0.8,0.8,0.8, 0.3), border=NA)
  lines(x=x.yrs, y=dat['Nup',yrs], lty=2, lwd=1, col=rgb(0.4, 0.4, 0.4))
  lines(x=x.yrs, y=dat['Nlow',yrs], lty=2, lwd=1, col=rgb(0.4, 0.4, 0.4))
  lines(x=x.yrs, y=dat['Nmean',yrs], lty=1, lwd=2, col=rgb(0.4, 0.4, 0.4))

  text(x=lastyear+1, y=ylim[2]-ypad*2, labels='FORECASTS >>', pos=4, cex=2, font=4, col=rgb(0.6,0.6, 0.6))
  
  axis(side=1, at=seq(1990, max(x.yrs), 10), labels=NA, lwd.ticks=3)
  axis(side=1, at=seq(1985, max(x.yrs), 10), labels=NA)
  
  # Carrying Capacity
  par(mar=c(1.5,7,1.5,1))

  xlim <- range(as.numeric(x.yrs))
  ylim <- c(0, max(10, min(1000, max(dat['kup',yrs],na.rm=T))))

  plot(NA, xlim=xlim, ylim=ylim,
       ylab='Carrying Capacity\n(fish per km)', xlab=NA,
       xaxt='n', cex.lab=2, cex.axis=2, lwd=2)
  
  xpad <- diff(range(xlim))*0.04; ypad <- diff(range(ylim))*0.04
  rect(xleft=lastyear+0.5, xright=xlim[2]+xpad, ybottom=ylim[1]-ypad, ytop=ylim[2]+ypad*0.85, col=rgb(0.9,0.9,0.9,0.5), border=NA)

  polygon(x=c(x.yrs, rev(x.yrs)), y=c(dat['kup',yrs], rev(dat['klow',yrs])), col=rgb(0.8,0.8,0.8, 0.3), border=NA)
  lines(x=x.yrs, y=dat['kup',yrs], lty=2, lwd=1, col=rgb(0.4, 0.4, 0.4))
  lines(x=x.yrs, y=dat['klow',yrs], lty=2, lwd=1, col=rgb(0.4, 0.4, 0.4))
  lines(x=x.yrs, y=dat['kmean',yrs], lty=1, lwd=2, col=rgb(0.4, 0.4, 0.4))

  axis(side=1, at=seq(1990, max(x.yrs), 10), labels=NA, lwd.ticks=3)
  axis(side=1, at=seq(1985, max(x.yrs), 10), labels=NA)
  
  # Growth Rate
  par(mar=c(4.5,7,1.5,1))

  dat['rup',yrs] <- (exp(dat['rup',yrs]) - 1) * 100
  dat['rmean',yrs] <- (exp(dat['rmean',yrs]) - 1) * 100
  dat['rlow',yrs] <- (exp(dat['rlow',yrs]) - 1) * 100

  xlim <- range(as.numeric(x.yrs))
  ylim <- c(min(dat['rlow',yrs],na.rm=T), max(dat['rup',yrs],na.rm=T))

  plot(NA, xlim=xlim, ylim=ylim,
       ylab='Growth Potential\n(%)', xlab='Year', xaxt='n',
       cex.lab=2, cex.axis=2, lwd=2)
  
  xpad <- diff(range(xlim))*0.04; ypad <- diff(range(ylim))*0.04
  rect(xleft=lastyear+0.5, xright=xlim[2]+xpad, ybottom=ylim[1]-ypad, ytop=ylim[2]+ypad*0.85, col=rgb(0.9,0.9,0.9,0.5), border=NA)

  polygon(x=c(x.yrs, rev(x.yrs)), y=c(dat['rup',yrs], rev(dat['rlow',yrs])), col=rgb(0.8,0.8,0.8, 0.3), border=NA)
  lines(x=x.yrs, y=dat['rup',yrs], lty=2, lwd=1, col=rgb(0.4, 0.4, 0.4))
  lines(x=x.yrs, y=dat['rlow',yrs], lty=2, lwd=1, col=rgb(0.4, 0.4, 0.4))
  lines(x=x.yrs, y=dat['rmean',yrs], lty=1, lwd=2, col=rgb(0.4, 0.4, 0.4))
  abline(h=0)
  
  axis(side=1, at=seq(1990, max(x.yrs), 10), lwd.ticks=3, cex.axis=2)
  axis(side=1, at=seq(1985, max(x.yrs), 10), labels=NA)
}