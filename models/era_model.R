
my.model <- function() {
  ## PRIORS
  psi.0              ~ dnorm(0,0.001)#each spp gets SAME occupancy but modify with psi.sp[species]
  psi.sitearea       ~ dnorm(0,0.001)
  p.0                ~ dnorm(0,0.001) #each spp gets same detection - no modification 

  ## Random effect of species on occupancy 
  sigma.psi.sp       ~ dunif(0,10)
  tau.psi.sp         <- 1/(sigma.psi.sp*sigma.psi.sp)
  for(species in 1:nspecies) { #for loop that creates a psi.sp for each species 
    psi.sp[species] ~ dnorm(0, tau.psi.sp) 
  } 
  
  ## Random slopes for the effect of era for each species  
  mu.psi.era           ~ dnorm(0,0.001)
  sigma.psi.era        ~ dunif(0,10)
  tau.psi.era           <- 1/(sigma.psi.era *sigma.psi.era )
  for(species in 1:nspecies) {  
    psi.era[species] ~ dnorm(mu.psi.era , tau.psi.era) 
  } 
  
  
  ## Random intercept of site on detection probability
  sigma.p.site         ~ dunif(0,10)
  tau.p.site           <- 1/(sigma.p.site*sigma.p.site)
  for(site in 1:nsite){
    for (era in 1:nera){
      p.site[site,era]   ~ dnorm(0, tau.p.site)
    }
  } 
  
  ## Model for Occupancy and detection (mean + spp-specific effect)
  for(site in 1:nsite) {
    for(species in 1:nspecies){
      for(era in 1:nera){ 
        logit(psi[site, era, species]) <- psi.0 +
          psi.sp[species]                            +
          psi.era[species]         * era             +
          psi.sitearea             * sitearea[site]      
        
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
  
  
}

## parameters to track
my.params <- c('psi.0', 
               
               'psi.sp',
               'sigma.psi.sp',
               
               'psi.era',
               'sigma.psi.era',
               'mu.psi.era',
               
               'psi.sitearea',
               
               'p.0', 
               'p.site',
               'sigma.p.site'
               )
