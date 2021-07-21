rm(list=ls())
gc()
library(leaflet);library(shiny);

# Constants
firstyear <- 1984
lastyear <- 2015

# Load R objects
load('./data/popnames.R')
load('./data/dthin.R')
load('./data/jags.dat.R')
load('./data/lastupdate.R')

thin <- function(d, mcmc){
  thin <- round(nrow(d)/mcmc, 0)
  dthin <- d[seq(1,nrow(d), by=thin),]
  return(dthin)
}

d <- thin(d, 5e3)
save(d, file='./data/dthin.R')

# Scale factors
scale.factors <- read.csv('./data/scale.factors.csv', stringsAsFactors=F)
scale.factors[scale.factors$cov=='bkt', 'mu'] <- 0
scale.factors[scale.factors$cov=='bkt', 'sd'] <- scale.factors[scale.factors$cov=='bkt', 'sd'] * 1000

# Covariates
extent.dat <- read.csv('./data/trimcovs/extent.dat.csv', stringsAsFactors=F)
templag.dat <- read.csv('./data/trimcovs/templag.dat.csv', stringsAsFactors=F)
hflowlag.dat <- read.csv('./data/trimcovs/hflowlag.dat.csv', stringsAsFactors=F)
bkt.dat <- read.csv('./data/trimcovs/bkt.dat.csv', stringsAsFactors=F)
ndvi.dat <- read.csv('./data/trimcovs/ndvi.dat.csv', stringsAsFactors=F)
add.dat <- read.csv('./data/trimcovs/add.dat.csv', stringsAsFactors=F)
rem.dat <- read.csv('./data/trimcovs/rem.dat.csv', stringsAsFactors=F)

row.names(extent.dat) <- extent.dat$PopulationName
row.names(templag.dat) <- templag.dat$PopulationName
row.names(hflowlag.dat) <- hflowlag.dat$PopulationName
row.names(bkt.dat) <- bkt.dat$PopulationName
row.names(ndvi.dat) <- ndvi.dat$PopulationName
row.names(add.dat) <- add.dat$PopulationName
row.names(rem.dat) <- rem.dat$PopulationName

names(templag.dat) <- sub('X','',names(templag.dat))
names(hflowlag.dat) <- sub('X','',names(hflowlag.dat))
names(bkt.dat) <- sub('X','',names(bkt.dat))
names(ndvi.dat) <- sub('X','',names(ndvi.dat))
names(add.dat) <- sub('X','',names(add.dat))
names(rem.dat) <- sub('X','',names(rem.dat))

# Non-natives
nonnatives <- read.csv('./data/nonnatives.csv', check.names=F)
rownames(nonnatives) <- nonnatives$PopulationName
nonnatives <- nonnatives[,-c(1,2)] * 1000
 
# Pop names
allpopnames <- extent.dat$PopulationName
recnames <- allpopnames[grepl('Reconnect', allpopnames, fixed=T)]

# Metapop names
metapop <- read.csv('./data/metapop.csv', stringsAsFactors=F)
row.names(metapop) <- metapop$PopulationName
metapop[metapop$MetaPop_Name=='','MetaPop_Name'] <- NA

for(popname in recnames){
  subpopnames <- metapop[metapop$MetaPop_Name==metapop[popname, 'MetaPop_Name'], 'PopulationName']
  subpopnames <- subpopnames[!grepl('Reconnect', subpopnames) & !is.na(subpopnames)]
  if(length(subpopnames) < 2) {
    allpopnames <- allpopnames[-which(allpopnames==popname)]
    recnames <- recnames[-which(recnames==popname)]
  }
}

gmu.list <- c("", unique(metapop$GMU))
state.list <- c("", unique(metapop$US_State))

# elevation
elevation <- read.csv('./data/elevation.csv')
for(row in row.names(elevation)){
  elevation[row,'PopulationName'] <- metapop[metapop$PopID==elevation[row,'PopID'], 'PopulationName']
}
elev <- (elevation$elev_m > 2000) + 1
names(elev) <- elevation$PopulationName
elevation <- elev
rm(elev)

# Locations
loc <- read.csv('./data/shapefiles/locations.csv')
row.names(loc) <- loc$PopulationName
loc <- loc[,c('POINT_X','POINT_Y')]

# extinctTable
extinctTable <- read.csv('./data/extinct.csv', stringsAsFactors=F)
row.names(extinctTable) <- extinctTable$PopulationName
extinctTable[,c('extinctlow','extinct','extinctup')] <- extinctTable[,c('extinctlow','extinct','extinctup')] * 100

# Streams
streams <- rgdal::readOGR('./data/shapefiles', layer='streams', verbose=F, stringsAsFactors=F, pointDropZ=T)
names(streams) <- c('PopID', 'Extent', 'PopulationName')
streams <- streams[c('PopID', 'PopulationName', 'Extent')]

streams.reconnect <- streams[streams$PopulationName %in% recnames,]

streams <- streams[streams$PopulationName %in% allpopnames,]

for(popname in streams$PopulationName[!streams$PopulationName %in% streams.reconnect$PopulationName]){
  streams[streams$PopulationName==popname, 'LCT'] <- sum(jags.dat$Yi[which(popnames==popname),], na.rm=T)
  streams[streams$PopulationName==popname, 'extinct'] <- extinctTable[popname,'extinct']
}
streams.lahontan <- streams[streams$LCT>0,]
streams.other <- streams[streams$LCT==0,]

pal <- colorBin(palette=c('blue','green','yellow','orange','red'), bins=c(0.0, 0.1, 0.3, 0.7, 0.9, 1.0)*100)

# Load functions
source('./functions/simpop.R')
source('./functions/simstat.R')
source('./functions/simcov.R')
source('./functions/simbatch.R')
source('./functions/plotpop.R')
source('./functions/plotcov.R')
source('./functions/ploteffect.R')
source('./functions/plotdet.R')
source('./functions/map.R')
source('./functions/rawdata.R')
source('./functions/resetButton.R')