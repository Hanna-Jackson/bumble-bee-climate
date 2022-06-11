## Written by: Hanna M. Jackson - hmj2@sfu.ca

## This file is sourced in the Main.R file just before a model is run 

## Dependencies
library(coda)
library(rjags)
library(R2jags)
library(hms) ## This is ONLY for getting the run time, not necesary 

## Some required functions 
expit <- function(x) 1/(1+exp(-x))
logit <- function(x) log(x/(1-x))  

## Load files
load(sprintf("data_prep/saved/for_analysis/JAGS.arr-%s.RData",   site.resolution), verbose=TRUE) 
load(sprintf('data_prep/saved/for_analysis/JAGS.visit-%s.RData', site.resolution), verbose=TRUE)
load(sprintf("data_prep/saved/for_analysis/JAGS.site-%s.RData",  site.resolution), verbose=TRUE)
load(sprintf("data_prep/saved/for_analysis/JAGS.era-%s.RData",   site.resolution), verbose=TRUE)

load(sprintf("data_prep/saved/for_analysis/species.ranges-%s.RData", site.resolution), verbose=TRUE)
load(sprintf("data_prep/saved/for_analysis/master.index-%s.RData",   site.resolution), verbose=TRUE)
load(sprintf("data_prep/saved/for_analysis/visit.arr-%s.RData",      site.resolution), verbose=TRUE) 


## Initial Values for the MCMC to start at
my.inits<-function(){
   tmp <- array(as.matrix(species.ranges*1), dim = c(nspecies, nsite, nera)) # takes species x site 1s and 0s (1= yes in species range according to 9.Make_Range.R) and then just duplicates it nera times
   tmp2 <- aperm(tmp, c(2,3,1))  # makes it so the dimensions are the same as in JAGS.arr 
   list(Z=tmp2)    
}

## Scaling function we will use for all of our variables
my.scale <- function(x){
  return((x - mean(x,na.rm=T)) / sd(x,na.rm=T))
}

## Modification to the data - scaling variables
JAGS.site[,"sitearea"]   <- my.scale(JAGS.site[,"sitearea"])                                 
JAGS.era[,,"meanmaxt"]   <- my.scale(JAGS.era[,,"meanmaxt"])
JAGS.era[,,"meanprec"]   <- my.scale(JAGS.era[,,"meanprec"])
JAGS.era[,,"floral.all"] <- my.scale(JAGS.era[,,"floral.all"])
                                    
## Making the final object that will go into JAGS for the run! 
my.data <- list(X              = JAGS.arr[master.index], 
                nsite          = dim(JAGS.arr)[1],
                nvisit         = dim(JAGS.arr)[2],
                nera           = dim(JAGS.arr)[3],
                nspecies       = dim(JAGS.arr)[4],
                nind           = nrow(master.index),
                master.index   = master.index,
                
                site           = site,
                sitevec        = site.vec,
                visit          = visit,
                visitvec       = visit.vec,
                era            = era,
                eravec         = era.vec,
                species        = species,
                
                sitearea       = JAGS.site[,"sitearea"],
                
                meanmaxt       = JAGS.era[,,"meanmaxt"],
                meanprec       = JAGS.era[,,"meanprec"],
                floral         = JAGS.era[,,"floral.all"],     
                
                species.ranges = species.ranges,               
                visit.arr      = visit.arr
                )


