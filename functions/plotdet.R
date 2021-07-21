plotdet <- function(popname, model=d, ci=0.95){
  # popname='Ackler'; model=d; ci=0.95
  
  probs <- c((1-ci)/2, 0.5, 1-(1-ci)/2)
  
  popnum <- which(popnames==popname)
  if(length(popnum)==0) popnum <- 0
  
  x <- 1:max(jags.dat$npasses, na.rm=T)
  
  par(mar=c(4,4,1,1)+4.5)
  
  plot(NA, xlim=range(x), ylim=c(0,1), ylab='Detection Probability', xlab='Electrofishing Pass', xaxt='n', cex.lab=2, cex.axis=1.5)
  axis(1, at=x, labels=x, cex.lab=2, cex.axis=1.5)
  
  obsmod <- function(m, b0, b1, x, delta){
    if(!is.finite(x)) x <- 0
    boot::inv.logit(b0 + b1*x) * exp(delta*(m-1))
  }
  
  if(popnum==0){
    # Mean drainage area
    dat <- data.frame(row.names=x, matrix(NA, ncol=3, nrow=max(x)))
    names(dat) <- c('low','mid','up')
    for (m in 1:max(x)){
      dat[m,] <- quantile(obsmod(m, b0=model$b0p, b1=model$b1p, 
                                 x=0,
                                 delta=model$delta),
                          probs=probs)
    }
    polygon(x=c(x, rev(x)), y=c(dat[,'up'], rev(dat[,'low'])), col=rgb(0.8,0.8,0.8, 0.3), border=NA)
    lines(x=x, y=dat[,'up'], lty=2, lwd=1, col=rgb(0.4, 0.4, 0.4))
    lines(x=x, y=dat[,'low'], lty=2, lwd=1, col=rgb(0.4, 0.4, 0.4))
    lines(x=x, y=dat[,'mid'], lty=1, lwd=2, col=rgb(0.4, 0.4, 0.4))
    
  } else {
    # Max drainage area
    dat <- data.frame(row.names=x, matrix(NA, ncol=3, nrow=max(x)))
    names(dat) <- c('low','mid','up')
    for (m in 1:max(x)){
      dat[m,] <- quantile(obsmod(m, b0=model$b0p, b1=model$b1p, 
                                 x=max(jags.dat$drain[popnum,,], na.rm=T),
                                 delta=model$delta),
                          probs=probs)
    }
    polygon(x=c(x, rev(x)), y=c(dat[,'up'], rev(dat[,'low'])), col=rgb(0,0,1, 0.3), border=NA)
    lines(x=x, y=dat[,'up'], lty=2, lwd=1, col=rgb(0, 0, 1))
    lines(x=x, y=dat[,'low'], lty=2, lwd=1, col=rgb(0, 0, 1))
    lines(x=x, y=dat[,'mid'], lty=1, lwd=2, col=rgb(0, 0, 1))
    
    # Min drainage area
    dat <- data.frame(row.names=x, matrix(NA, ncol=3, nrow=max(x)))
    names(dat) <- c('low','mid','up')
    for (m in 1:max(x)){
      dat[m,] <- quantile(obsmod(m, b0=model$b0p, b1=model$b1p, 
                                 x=min(jags.dat$drain[popnum,,], na.rm=T),
                                 delta=model$delta),
                          probs=probs)
    }
    polygon(x=c(x, rev(x)), y=c(dat[,'up'], rev(dat[,'low'])), col=rgb(1,0,0, 0.3), border=NA)
    lines(x=x, y=dat[,'up'], lty=2, lwd=1, col=rgb(1,0,0))
    lines(x=x, y=dat[,'low'], lty=2, lwd=1, col=rgb(1,0,0))
    lines(x=x, y=dat[,'mid'], lty=1, lwd=2, col=rgb(1,0,0))
    
    bigkm <- round(max(jags.dat$drain[popnum,,], na.rm=T) * scale.factors[scale.factors$cov=='drain','sd'] + scale.factors[scale.factors$cov=='drain','mu'], 1)
    smallkm <- round(min(jags.dat$drain[popnum,,], na.rm=T) * scale.factors[scale.factors$cov=='drain','sd'] + scale.factors[scale.factors$cov=='drain','mu'], 1)
    legend('topright', legend=c(paste('Smallest sampled site (catchment =',smallkm,'sq km)'), paste('Largest sampled site (catchment =',bigkm,'sq km)')), 
           border=c('red','blue'), fill=c(rgb(1,0,0, 0.3),rgb(0,0,1, 0.3)),
           cex=1.5, box.lwd=0, bg=rgb(1,1,1,0.5), xjust=0.5)  
  }
  
  
}