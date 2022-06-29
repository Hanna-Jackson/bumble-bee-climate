
## This model has the same structure as both era_model.R and env_model.R,
## the only difference is that is includes the predictors from both of
## them. So occupancy is a function of era, temperature, precipitation
## and floral resources 

## NOTE: Please see era_model.R and env_model.R for detailed comments
## describing what each section of code does

my.model <- function() {
 ## ~~~~~~~~~~~~~~ Priors ~~~~~~~~~~~~~~~
  psi.0              ~ dnorm(0,0.001)
  psi.sitearea       ~ dnorm(0,0.001) 
  p.0                ~ dnorm(0,0.001) 
  
  ## ~~~ Random effect of species on occupancy intercept ~~~
  sigma.psi.sp       ~ dunif(0,10)
  tau.psi.sp         <- 1/(sigma.psi.sp*sigma.psi.sp)
  for(species in 1:nspecies) {  
    psi.sp[species] ~ dnorm(0, tau.psi.sp) 
  } 
  
  ## ~~~ Random intercept for each site-era's detection probability ~~~
  sigma.p.site         ~ dunif(0,10)
  tau.p.site           <- 1/(sigma.p.site*sigma.p.site)
  for(site in 1:nsite){
    for (era in 1:nera){
      p.site[site,era]   ~ dnorm(0, tau.p.site)
    }
  } 
  
  ## ~~~ Random slopes for the effect of temperature on each species'occupancy ~~~
  mu.psi.meanmaxt ~ dnorm(0,0.001)
  sigma.psi.meanmaxt ~ dunif(0,100)
  tau.psi.meanmaxt <- 1/(sigma.psi.meanmaxt*sigma.psi.meanmaxt)
  for(species in 1:nspecies) {
    psi.meanmaxt[species] ~ dnorm(mu.psi.meanmaxt, tau.psi.meanmaxt)
  }

  ## ~~~ ALL species meanmaxt.sq parameter ~~~
  psi.meanmaxt.sq ~  dnorm(0,0.001)  

  ## ~~~ Random slopes for the effect of precipitation on each species' occupancy ~~~
  mu.psi.meanprec           ~ dnorm(0,0.001)
  sigma.psi.meanprec        ~ dunif(0,10)
  tau.psi.meanprec          <- 1/(sigma.psi.meanprec *sigma.psi.meanprec )
  for(species in 1:nspecies) {  
    psi.meanprec[species] ~ dnorm(mu.psi.meanprec, tau.psi.meanprec ) 
  }

  ## ~~~ Random slopes for the effect of floral resources on each species' occupancy ~~~
  mu.psi.floral         ~ dnorm(0,0.001)
  sigma.psi.floral      ~ dunif(0,10)
  tau.psi.floral         <- 1/(sigma.psi.floral *sigma.psi.floral )
  for(species in 1:nspecies) {  
    psi.floral[species] ~ dnorm(mu.psi.floral, tau.psi.floral ) 
  } 

  ## ~~~ Random slopes for the effect of era on each species' occupancy ~~~
  mu.psi.era           ~ dnorm(0,0.001)
  sigma.psi.era        ~ dunif(0,10)
  tau.psi.era           <- 1/(sigma.psi.era *sigma.psi.era )
  for(species in 1:nspecies) {  
    psi.era[species] ~ dnorm(mu.psi.era , tau.psi.era) 
  } 

  
  ## ~~~~~~~~~~~~~~ Model ~~~~~~~~~~~~~~~
  for(site in 1:nsite) {
    for(species in 1:nspecies){
      for(era in 1:nera){ 
        logit(psi[site, era, species]) <- psi.0          +
          psi.sp[species]                                +
          psi.era[species]        * era                  + 
          psi.meanmaxt[species]   * meanmaxt[site,era]   +
          psi.meanmaxt.sq         * meanmaxt[site,era]^2 +                                       
          psi.meanprec[species]   * meanprec[site,era]   +
          psi.floral[species]     * floral[site,era]     +
          psi.sitearea            * sitearea[site]      
        
        for(visit in 1:nvisit){
          logit(p[site, visit, era, species]) <- p.0 + 
            p.site[site, era]
          
        }
      } 
    } 
  } 

  
 ## ~~~~~~~~~~~~~~ Likelihood ~~~~~~~~~~~~~~~
  for(site in 1:nsite) {
    for(era in 1:nera) {
      for(species in 1:nspecies){
        Z[site, era, species] ~ dbern(psi[site, era, species]) 
      } 
    }
  }
  for(ind in 1:nind){
    p.eff[site[ind],  visit[ind], era[ind], species[ind]] <- Z[site[ind], era[ind], species[ind]] * p[site[ind], visit[ind], era[ind], species[ind]] 
    X[ind] ~ dbern(p.eff[site[ind],  visit[ind], era[ind], species[ind]])
  } 
  
  
} ## End of my.model 


## ~~~~~~~~~~~~~~ Parameters to Track  ~~~~~~~~~~~~~~~
my.params <- c('psi.0', 
               
               'psi.sp',
               'sigma.psi.sp',

               'psi.era',
               'mu.psi.era',
               'sigma.psi.era',
               
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


