
my.model <- function() {
  ## PRIORS
  psi.0              ~ dnorm(0,0.001) # prior for the intercept 
  psi.sitearea       ~ dnorm(0,0.001) # prior for fixed effect of site area
  p.0                ~ dnorm(0,0.001) #each spp gets same detection - no modification 
  
  ## Random effect of species on occupancy
  sigma.psi.sp       ~ dunif(0,10)
  tau.psi.sp         <- 1/(sigma.psi.sp*sigma.psi.sp)
  for(species in 1:nspecies) { #for loop that creates a psi.sp for each species 
    psi.sp[species] ~ dnorm(0, tau.psi.sp) 
  } 
  
  ## Random effect of site on  detection probability
  sigma.p.site         ~ dunif(0,10)
  tau.p.site           <- 1/(sigma.p.site*sigma.p.site)
  for(site in 1:nsite){
    for (era in 1:nera){
      p.site[site,era]   ~ dnorm(0, tau.p.site)
    }
  } 
  
  ## Random slopes for the effect of meanmaxt for each species
  mu.psi.meanmaxt ~ dnorm(0,0.001)
  sigma.psi.meanmaxt ~ dunif(0,100)
  tau.psi.meanmaxt <- 1/(sigma.psi.meanmaxt*sigma.psi.meanmaxt)
  for(species in 1:nspecies) {
    psi.meanmaxt[species] ~ dnorm(mu.psi.meanmaxt, tau.psi.meanmaxt)
  }

  ## ALL species meanmaxt.sq parameter 
  psi.meanmaxt.sq ~  dnorm(0,0.001)  

  ## Random slopes for the effect of  meanprec for each species
  mu.psi.meanprec           ~ dnorm(0,0.001)
  sigma.psi.meanprec        ~ dunif(0,10)
  tau.psi.meanprec          <- 1/(sigma.psi.meanprec *sigma.psi.meanprec )
  for(species in 1:nspecies) {  
    psi.meanprec[species] ~ dnorm(mu.psi.meanprec, tau.psi.meanprec ) 
  }

  ## Random slopes for the effect of floral resources for each species
  mu.psi.floral         ~ dnorm(0,0.001)
  sigma.psi.floral      ~ dunif(0,10)
  tau.psi.floral         <- 1/(sigma.psi.floral *sigma.psi.floral )
  for(species in 1:nspecies) {  
    psi.floral[species] ~ dnorm(mu.psi.floral, tau.psi.floral ) 
  } 

  ## Random slopes for the effect of era for each species
  mu.psi.era           ~ dnorm(0,0.001)
  sigma.psi.era        ~ dunif(0,10)
  tau.psi.era           <- 1/(sigma.psi.era *sigma.psi.era )
  for(species in 1:nspecies) {  
    psi.era[species] ~ dnorm(mu.psi.era , tau.psi.era) 
  } 

  
  ## Model for Occupancy and detection (mean + spp-specific effect)
  for(site in 1:nsite) {
    for(species in 1:nspecies){
      for(era in 1:nera){ 
        logit(psi[site, era, species]) <- psi.0           +
          psi.sp[species]                                 +
          psi.era[species]        * era                   + 
          psi.meanmaxt[species]   * meanmaxt[site, era]   +
          psi.meanmaxt.sq         * meanmaxt[site,era]^2  +                                       
          psi.meanprec[species]   * meanprec[site, era]   +
          psi.floral[species]     * floral[site,era]      +
          psi.sitearea            * sitearea[site]      
        
        for(visit in 1:nvisit){
          logit(p[site, visit, era, species]) <- p.0 + 
            p.site[site, era]
          
        }
      } 
    } 
  } 
  
  ## LIKELIHOOD
  for(site in 1:nsite) {
    for(era in 1:nera) {
      for(species in 1:nspecies){
        ## Occurrence 
        Z[site, era, species] ~ dbern(psi[site, era, species]) 
      } 
    }
  }
  for(ind in 1:nind){
    p.eff[site[ind],  visit[ind], era[ind], species[ind]] <- Z[site[ind], era[ind], species[ind]] * p[site[ind], visit[ind], era[ind], species[ind]] 
    X[ind] ~ dbern(p.eff[site[ind],  visit[ind], era[ind], species[ind]])
  } 
  
  
}#model 

## parameters to track
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


