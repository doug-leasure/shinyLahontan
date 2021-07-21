rawdata <- function(popname, jd=jags.dat, poplist=popnames){
  # popname='Abel'; jd=jags.dat; poplist=popnames;
  
  if (popname %in% popnames){
    popnum <- which(popnames==popname)
    
    nrows <- sum(!is.na(jd$y[popnum,,,]))
    
    dat <- data.frame(row.names=1:nrows, Year=rep(NA,nrows), Site=rep(NA,nrows), Length_m=rep(NA,nrows), Pass=rep(NA,nrows), Fish=rep(NA,nrows))
    row <- 0
    
    for(ti in 1:jd$nt[popnum]){
      year <- firstyear + jd$t_ti[popnum, ti] - 1
      
      for(site in 1:jd$nj[popnum,ti]){
        for(pass in 1:jd$npasses[popnum, ti, site]){
          row <- row + 1
          dat[row,] <- c(Year = year, 
                         Site = site, 
                         Length_m = round(jd$extent[popnum]*(jd$sitelength[popnum,ti,site]/jd$extent[popnum])*1000), 
                         Pass = pass, 
                         Fish = jd$y[popnum,ti,site,pass]
                         )
    }}}
    return(dat)
  }
}