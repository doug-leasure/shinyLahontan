shinyServer(
function(input, output, session){
  
  # Thin MCMC
  dthin <- reactive( thin(d, input$mcmc) )
  
  # Simulate covariates
  covs <- reactive( simcov(input$popname, forecastyear=input$forecastyear, nsim=input$nsim,
                           hist.bkt=input$hist.bkt, extent=input$extent,
                           cut=list(templag=input$cut.temp, hflowlag=input$cut.hflow, ndvi=input$cut.ndvi, bkt=c(0,1)), 
                           const=list(bkt=input$const.bkt, reintro=c(rep(input$add1, input$add1t), rep(input$add2, input$add2t)), templag=NA, hflowlag=NA, ndvi=NA)) )
  
  # Simulate population
  simN <- reactive( simpop(input$popname, covs=covs(), forecastyear=input$forecastyear, 
                          demsto=input$demsto, ci=input$ci,
                          model=dthin()) )
  
  # Stats from simulations for plots
  sim <- reactive( simstat(input$popname, simN=simN(), forecastyear=input$forecastyear, ci=input$ci) )
  
  # Plot population
  output$simPlot <- renderPlot( plotpop(input$popname, sim()) )
  
  # Plot covariates
  output$covPlot <- renderPlot( plotcov(input$popname, forecastyear=input$forecastyear, covs=covs(), sim=sim(), ci=input$ci) )
  
  # Plot effects
  output$effectPlot <- renderPlot( ploteffect(input$popname, covs=covs(), sim=sim(), 
                                              forecastyear=input$forecastyear, ci=input$ci,
                                              model=dthin()) )
  # Plot detection
  output$detectPlot <- renderPlot( plotdet(input$popname, model=dthin(), ci=input$ci) )
  
  # Map
  output$map <- renderLeaflet(map())
  
  observeEvent(input$popname, proxy(input$popname))
  
  observeEvent(input$map_shape_click$id, {
    if(input$map_shape_click$id %in% allpopnames) updateSelectInput(session, 'popname', selected=input$map_shape_click$id) 
  })
  
  # Raw data
  output$dataTable <- renderTable( rawdata(input$popname), striped=T, digits=0, na="", align='c' )
  
  # Batch
  output$batchTable <- renderTable({
    input$batchButton
    isolate(simbatch(primepop=input$popname, pops=input$batchnames, forecastyear=input$forecastyear, nsim=input$nsim, nmcmc=input$mcmc,
                     add1=input$add1, add1t=input$add1t, add2=input$add2, add2t=input$add2t,                      
                     hist.bkt=input$hist.bkt, extent=input$extent, 
                     cut=list(templag=input$cut.temp, hflowlag=input$cut.hflow, ndvi=input$cut.ndvi, bkt=c(0,1)), 
                     const=list(bkt=input$const.bkt, reintro=c(rep(input$add1, input$add1t), rep(input$add2, input$add2t)), templag=NA, hflowlag=NA, ndvi=NA),
                     demsto=NA, ci=input$ci, model=dthin()))
  },striped=T, digits=3, na="", align='c')
  
  observeEvent(input$batchGMU, {
    if(!input$batchGMU==""){
      updateSelectizeInput(session, 'batchState', selected="")
      updateSelectInput(session, 'batchnames', selected=metapop[metapop$GMU==input$batchGMU, 'PopulationName'])
      updateSelectizeInput(session, 'batchGMU', selected="")
    }
  })
  
  observeEvent(input$batchState, {
    if(!input$batchState==""){
      updateSelectizeInput(session, 'batchGMU', selected="")
      updateSelectInput(session, 'batchnames', selected=metapop[metapop$US_State==input$batchState, 'PopulationName'])
      updateSelectizeInput(session, 'batchState', selected="")
    }
  })
  
  observeEvent(input$batchRemoveNoData, {
    updateSelectInput(session, 'batchnames', selected=input$batchnames[input$batchnames %in% streams.lahontan$PopulationName])
  })
  
  observeEvent(input$batchRemovewData, {
    updateSelectInput(session, 'batchnames', selected=input$batchnames[!input$batchnames %in% streams.lahontan$PopulationName])
  })
  
  observeEvent(input$batchClear, {
    updateSelectInput(session, 'batchnames', selected="")
  })
  
  # Reset button
  observe( resetButton(session, input$reset, input$popname, dthin()) )
  
  # Dynamic default values
  observeEvent(input$popname, {
    updateSliderInput(session, 'const.bkt', value=bkt.dat[input$popname, as.character(lastyear)]*scale.factors[scale.factors$cov=='bkt','sd'])
    
    if(input$popname %in% recnames) {
      subpopnames <- metapop[metapop$MetaPop_Name==metapop[input$popname, 'MetaPop_Name'], 'PopulationName']
      subpopnames <- subpopnames[!grepl('Reconnect', subpopnames) & !is.na(subpopnames)]
      subpopnums <- which(popnames %in% subpopnames)
    } else {subpopnums <- c()}
    
    if(input$popname %in% popnames){
      updateSliderInput(session, 'demsto', value=median(dthin()[,paste('sigmaR[',which(popnames==input$popname),']',sep='')]))
    } 
    else if(input$popname %in% recnames & length(subpopnums)>0) {
      sigmaR <- apply(as.matrix(dthin()[,paste('sigmaR[',subpopnums,']',sep='')]), 1, mean)
      updateSliderInput(session, 'demsto', value=median(sigmaR))
    } 
    else {
      updateSliderInput(session, 'demsto', value=median(dthin()[,'musig']))
    }
  })
})









