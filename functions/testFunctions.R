rm(list=ls())

setwd('C:/RESEARCH/2015 NASA Trout/wd/lahontan')

source('./global.R')

# simbatch(primepop='Abel', pops=c('Abel','Andorno'))

popname <- 'Abel'

covs <- simcov(popname, nsim=10)
simN <- simpop(popname, covs)
sim <- simstat(popname, simN, covs)

plotcov(popname, covs, sim)
plotpop(popname, sim)
ploteffect(popname, covs, sim)



####################33
plot(y=apply(simN$r[1,,], 2, mean), x=1:dim(simN$r)[3], type='l')
plot(y=sim['rmean',], x=1:ncol(sim), type='l')
