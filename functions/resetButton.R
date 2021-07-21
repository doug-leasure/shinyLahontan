resetButton <- function(session, button, popname, dthin=NA){
  if(!is.null(button)) {
    
    updateNumericInput(session, 'extent', value=0)
    updateSliderInput(session, 'const.bkt', value=bkt.dat[popname, as.character(lastyear)]*scale.factors[scale.factors$cov=='bkt','sd'])
    
    if(popname %in% recnames){
      subpopnames <- metapop[metapop$MetaPop_Name==metapop[popname, 'MetaPop_Name'], 'PopulationName']
      subpopnames <- subpopnames[!grepl('Reconnect', subpopnames) & !is.na(subpopnames)]
      subpopnums <- which(popnames %in% subpopnames)
    } else {subpopnums <- c()}
    
    if(popname %in% popnames) { 
      updateSliderInput(session, 'demsto', value=mean(d[,paste('sigmaR[',which(popnames==popname),']',sep='')]))
    } 
    else if(popname %in% recnames & length(subpopnums)>0){
      sigmaR <- apply(as.matrix(dthin[,paste('sigmaR[',subpopnums,']',sep='')]), 1, mean)
      updateSliderInput(session, 'demsto', value=mean(sigmaR))
    } 
    else {
      updateSliderInput(session, 'demsto', value=mean(d[,'musig']))  
    }

    updateRadioButtons(session, 'hist.bkt', selected=FALSE)
    updateSliderInput(session, 'cut.temp', value=c(0,1))
    updateSliderInput(session, 'cut.hflow', value=c(0,1))
    updateSliderInput(session, 'cut.ndvi', value=c(0,1))
    # updateCheckboxInput(session, 'natvar', value=FALSE)
    updateSliderInput(session, 'futrscale', value=1)
    updateNumericInput(session, 'add1', value=0)
    updateNumericInput(session, 'add1t', value=1)
    updateNumericInput(session, 'add2', value=0)
    updateNumericInput(session, 'add2t', value=99)
    
    # updateSliderInput(session, 'nsim', value=10)
    # updateSliderInput(session, 'mcmc', value=500)
    # updateSliderInput(session, 'ci', value=0.95)
  }
  return(session)
}