
## This model uses environmental covariates as predictors of
## occupancy.
##    meanmaxt = temperature
##    meanprec = precipitation
##    floral   = floral resources 

my.model <- function() {
  
  ## ~~~~~~~~~~~~~~ Priors ~~~~~~~~~~~~~~~
  psi.0              ~ dnorm(0,0.001) ## all species share this occupancy intercept 
  psi.sitearea       ~ dnorm(0,0.001) ## prior for the effect of site area on occupancy
  p.0                ~ dnorm(0,0.001) ## all species share gets same detection intercept
  
  ## ~~~ Random effect of species on occupancy intercept ~~~
  sigma.psi.sp       ~ dunif(0,10)
  tau.psi.sp        <- 1/(sigma.psi.sp*sigma.psi.sp)
  for(species in 1:nspecies) { 
    psi.sp[species]  ~ dnorm(0, tau.psi.sp) 
  } 
  
  ## ~~~ Random intercept for each site-era's detection probability ~~~
  ##     This means that each site gets it's own detection probability
  ##     in each of the 6 eras. This adds a high degree of flexibility
  ##     (and a large number of parameters) to the detection
  ##     probability model. 
  sigma.p.site           ~ dunif(0,10)
  tau.p.site            <- 1/(sigma.p.site*sigma.p.site)
  for(site in 1:nsite){
    for (era in 1:nera){
      p.site[site,era]   ~ dnorm(0, tau.p.site)
    }
  }

  ## ~~~ ALL species share a  meanmaxt.sq parameter ~~~
  psi.meanmaxt.sq ~ dnorm(0,0.001)
  
  ## ~~~ Random slopes for the effects of our environmental variables on each species' occupancy ~~~
  ##     This means that each species' occupancy can be affected by
  ##     each of our variables in different ways! So for temperature as
  ##     an example, some species might be increasing in occupancy with
  ##     increasing temperature (positive value of psi.meanmaxt for
  ##     that species and some might be decreasing (negative value of
  ##     psi.meanmaxt for that species)

  ## Random slopes for the effect of temperature on each species' occupancy
  mu.psi.meanmaxt    ~ dnorm(0,0.001)
  sigma.psi.meanmaxt ~ dunif(0,100)
  tau.psi.meanmaxt <- 1/(sigma.psi.meanmaxt*sigma.psi.meanmaxt)
  for(species in 1:nspecies) {
    psi.meanmaxt[species] ~ dnorm(mu.psi.meanmaxt, tau.psi.meanmaxt)
  }
  
 ## Random slopes for the effect of precipitation on  each species' occupancy
  mu.psi.meanprec    ~ dnorm(0,0.001)
  sigma.psi.meanprec ~ dunif(0,10)
  tau.psi.meanprec    <- 1/(sigma.psi.meanprec *sigma.psi.meanprec )
  for(species in 1:nspecies) {  
    psi.meanprec[species] ~ dnorm(mu.psi.meanprec, tau.psi.meanprec ) 
  }

  ## Random slopes for the effect of floral resources for each species' occupancy
  mu.psi.floral     ~ dnorm(0,0.001)
  sigma.psi.floral  ~ dunif(0,10)
  tau.psi.floral <- 1/(sigma.psi.floral *sigma.psi.floral )
  for(species in 1:nspecies) {  
    psi.floral[species] ~ dnorm(mu.psi.floral, tau.psi.floral ) 
  } 

  
  ## ~~~~~~~~~~~~~~ Model ~~~~~~~~~~~~~~~
  for(site in 1:nsite) {
    for(species in 1:nspecies){
      for(era in 1:nera){ 
        logit(psi[site, era, species]) <- psi.0          + ## The occupancy probability for a given species at a given site in a given era is equal to an intercept...
          psi.sp[species]                                + ## Plus a species-specific intercept
          psi.meanmaxt[species]   * meanmaxt[site,era]   + ## Plus a species-specific effect of temperature
          psi.meanmaxt.sq         * meanmaxt[site,era]^2 + ## Plus a species-specific effect of precipitation
          psi.meanprec[species]   * meanprec[site,era]   + ## Plus a species-specific effect of floral resources
          psi.floral[species]     * floral  [site,era]   + ## Plus a species-specific effect of floral resources
          psi.sitearea            * sitearea[site]         ## Plus the effect of the area of that site (accounts for fact that bigger sites have higher occupancy)
       
        
        for(visit in 1:nvisit){
          logit(p[site, visit, era, species]) <- p.0 +      ## The detection probability for a given species at a given site in a given visit interval in a given era interval is equal to an intercept...
            p.site[site, era]                               ## Plus a random effect of site-era (so each site in each of the 6 eras can have its own detection probability)
          
        }
      } 
    } 
  } 

  
  ## ~~~~~~~~~~~~~~ Likelihood ~~~~~~~~~~~~~~~
   ## This section describes how occupancy (psi) and detection probability (p) are related to our data (X) via the true, unknown occupancy state (Z)
  for(site in 1:nsite) {
    for(era in 1:nera) {
      for(species in 1:nspecies){
        Z[site, era, species] ~ dbern(psi[site, era, species]) ## True occupancy (0 or 1) is a Bernouli draw (a single weighted coin flip) with probability of a 1 equal to the occupancy probability of that species at that site in that era
      } 
    }
  }
  for(ind in 1:nind){
    p.eff[site[ind],  visit[ind], era[ind], species[ind]] <- Z[site[ind], era[ind], species[ind]] * p[site[ind], visit[ind], era[ind], species[ind]] ## "Realized" detection probability is equal to true occupancy (0 or 1, is the site occupuied) times detection probability
    X[ind] ~ dbern(p.eff[site[ind],  visit[ind], era[ind], species[ind]])  ## Our observations (0 or 1, did we see that species at that site in that era) are a Bernoulli draw from "Realized" detection probability)
  } 
   
} ## End of my.model 


## ~~~~~~~~~~~~~~ Parameters to Track  ~~~~~~~~~~~~~~~
## This section just tells JAGS which parameter values we want it to keep track of
my.params <- c('psi.0', 
               
               'psi.sp',
               'sigma.psi.sp',
               
               'psi.meanmaxt',
               'mu.psi.meanmaxt',
               'sigma.psi.meanmaxt',

               'psi.meanmaxt.sq',
               
               'psi.meanprec',
               'mu.psi.meanprec',
               'sigma.psi.meanprec',
               
               'psi.floral',
               'mu.psi.floral',
               'sigma.psi.floral',
               
               'psi.sitearea',
               
               'p.0', 
               'p.site',
               'sigma.p.site'
               )


