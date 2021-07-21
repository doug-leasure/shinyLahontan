shinyUI(
fluidPage(
  
  titlePanel("Lahontan Cutthroat Trout Population Simulator (v2.01 beta)"),
  
  sidebarLayout(position='left',
    sidebarPanel(style='background:white; border:none',
      
      wellPanel(style='margin:0;',
        selectInput('popname', 'Population:', choices=allpopnames),
        sliderInput('forecastyear', 'Forecast year:', min=2020, max=2115, value=2045, step=1, sep='')
      ),
      
      br(),
      
      wellPanel(style='margin:0; overflow-y:scroll; max-height:500px',
        div(style='display:inline-block', h4('Future Conditions')),
        div(style='display:inline-block; float:right', actionButton('reset','Reset')),
        p(),
        
        numericInput('extent', 'Change in population extent (km):', value=0, min=-100, max=100, width=250, step=0.1),
        br(),
        
        radioButtons('hist.bkt', "Non-native trout:", 
                     choices=c('Use historic densities'=TRUE, 'Set a constant density (slider below)'=FALSE),
                     selected=FALSE),
        sliderInput('const.bkt', NULL, min=0, max=2000, step=10, value=0, width='100%', post=' per km'),
        br(),
        
        sliderInput('demsto', 'Environmental stochasticity:', value=1, min=0.25, max=2.0, round=2),
        br(),
        
        strong('Stream Habitat:'),br(),
        'Sliders define the range of historical values used for forecasting',br(),
        '(0 = min, 0.5 = median, 1 = max). Open the Habitat tab to watch changes.',br(),
        br(),
        'Temperature',
        sliderInput('cut.temp', NULL, min=0, max=1, value=c(0,1)),
        'High flow',
        sliderInput('cut.hflow', NULL, min=0, max=1, value=c(0,1)),
        'Riparian vegetation (NDVI)',
        sliderInput('cut.ndvi', NULL, min=0, max=1, value=c(0,1)),
        
        # br(),
        # strong(checkboxInput('natvar', label='Manipulate habitat variability (slider below)', value=F)),
        # sliderInput('futrscale', NULL, min=0, max=2, step=0.1,  value=1, width='100%'),
        
        hr(),br(),
        h4('Reintroductions'),
        div(style='display:inline-block', numericInput('add1', 'Add x fish/year', value=0, min=0, max=999)),
        div(style='display:inline-block; margin-left:10px', numericInput('add1t', 'for t years,', value=1, min=0, max=99)),
        br(),'followed by...', br(),
        div(style='display:inline-block', numericInput('add2', 'x fish/year', value=0, min=0, max=999)),
        div(style='display:inline-block; margin-left:10px', numericInput('add2t', 'for t years', value=99, min=0, max=99)),
        
        br(),'* Assumes introduced fish are reproductive adults that all survive.',
        br(),'* Reintroduction schedule will be repeated throughout simulation.',
        br(),'* Add negative numbers of fish to simulate removals.',
        
        hr(),br(),
        h4('Batch Processing'),
        '* View results in "Batch" tab.',br(),
        '* See "Help" tab for help.',br(),
        div(style='display:inline-block; float:right', actionButton('batchButton','Run simulations...')),
        br(),br(),
        selectInput('batchnames', 'Manually select multiple populations:', choices=allpopnames, multiple=T),
        div(style='display:inline-block', selectizeInput('batchGMU', 'Select by GMU', choices=gmu.list, multiple=F, width=200)),
        div(style='display:inline-block; margin-left:10px', selectizeInput('batchState', 'Select by State', choices=state.list, multiple=F, width=200)),
        br(),
        strong('Remove from list:'),br(),
        div(style='display:inline-block', actionButton('batchRemoveNoData', 'Unoccupied streams')), 
        div(style='display:inline-block', actionButton('batchRemovewData', 'Occupied streams')),
        div(style='display:inline-block', actionButton('batchClear','Clear list')),
        
        
        hr(),br(),
        h4('Simulation Settings'),
        sliderInput('ci', 'Credible intervals:', min=0.50, max=0.99, value=0.95, step=0.05),
        sliderInput('nsim', 'Number of simulations (more is slower):', min=1, max=100, value=10, step=1),
        sliderInput('mcmc', 'Number of MCMC samples (more is slower):', min=1000, max=nrow(d), value=1000, step=100),
        
        br(),br()
      )
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel('Map', leafletOutput('map', height=700)),
        tabPanel('Population', plotOutput('simPlot', height=700)),
        tabPanel('Habitat', plotOutput('covPlot', height=700)),
        tabPanel('Effects', plotOutput('effectPlot', height=700)),
        tabPanel('Detection', plotOutput('detectPlot', height=700)),
        tabPanel('Data', tableOutput('dataTable'), style='overflow-y:scroll; max-height:700px'),
        tabPanel('Batch', tableOutput('batchTable'), style='overflow-y:scroll; max-height:700px'),
        tabPanel('Help', includeHTML('help.html'), style='overflow-y:scroll; max-height:700px'),
        tabPanel('About', includeHTML('about.html'), style='overflow-y:scroll; max-height:700px')
      ),
      hr(),
      helpText(paste('Last Update:', lastupdate))
    )
  )
))