map <- function() {
  
  leaflet(options = leafletOptions(minZoom=5, maxZoom=18)) %>%
      
    addProviderTiles(provider='Esri.WorldImagery', group='Satellite') %>%
      
    addProviderTiles(provider= 'Esri.NatGeoWorldMap', group='Map') %>%
     
    addPolylines(weight=2, color=~pal(streams.lahontan$extinct), opacity=1, data=streams.lahontan, group='Lahontan streams', 
                 layerId=streams.lahontan$PopulationName, label=streams.lahontan$PopulationName, 
                 labelOptions=labelOptions(noHide=TRUE, textsize=12) ) %>%
      
    addPolylines(weight=2, color='white', opacity=1, data=streams.other, group='Potential streams', 
                 layerId=streams.other$PopulationName, label=streams.other$PopulationName, 
                 labelOptions=labelOptions(noHide=TRUE, textsize=12) ) %>%
    
    addPolylines(weight=3, color='black', dashArray=c(5,5), opacity=1, data=streams.reconnect, group='Reconnect corridors', 
                 layerId=streams.reconnect$PopulationName, label=streams.reconnect$PopulationName, 
                 labelOptions=labelOptions(noHide=TRUE, textsize=12) ) %>%
    
    addLayersControl(baseGroups=c('Satellite','Map'), 
                     overlayGroups=c('Lahontan streams','Potential streams','Reconnect corridors'), 
                     options=layersControlOptions(collapsed=FALSE)) %>%
    
    addLegend(position='bottomright', pal=pal, values=extinctTable$extinct, 
              title='Extinction Risk (30 yr)',
              labFormat=labelFormat(suffix='%')
              ) %>%
      
    setMaxBounds(lng1=streams@bbox['x','min']-0.5*(streams@bbox['x','max']-streams@bbox['x','min']),
                lng2=streams@bbox['x','max']+0.5*(streams@bbox['x','max']-streams@bbox['x','min']),
                lat1=streams@bbox['y','min']-0.25*(streams@bbox['y','max']-streams@bbox['y','min']),
                lat2=streams@bbox['y','max']+0.25*(streams@bbox['y','max']-streams@bbox['y','min']) ) %>%
      
    addEasyButton(easyButton(
        icon="fa-globe", title="Zoom Out",
        onClick=JS("function(btn, map){ map.setZoom(7); }"))) %>%
      
    addScaleBar(position='topleft')
}

proxy <- function(popname){
  
  stream <- subset(streams, streams$PopulationName==popname)

  if (popname %in% streams.reconnect$PopulationName) {
    col <- 'black'
    dash <- c(5,5)
  } else if(popname %in% streams.other$PopulationName){
    col <- 'white'
    dash <- NULL
  } else {
    col <- pal(stream$extinct)
    dash <- NULL
  }
    
  leafletProxy('map') %>% 
    
    clearGroup('selected') %>%
    
    addPolylines(weight=5, color=col, dashArray=dash, opacity=1, data=stream, label=popname, group='selected' ) %>%
    
    addMarkers(lng = loc[popname,'POINT_X'], lat = loc[popname,'POINT_Y'], 
               label = paste(popname, '(', extent.dat[popname,'extent'],' km)'),
               labelOptions = list(noHide=TRUE, textsize=14), group='selected') %>%
    
    fitBounds(lng1=stream@bbox['x','min'],
              lng2=stream@bbox['x','max'],
              lat1=stream@bbox['y','min'],
              lat2=stream@bbox['y','max'] )
}
