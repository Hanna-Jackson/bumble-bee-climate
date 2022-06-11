## Written by: Hanna M. Jackson - hmj2@sfu.ca

make.covariates <- function(site.resolution,
                            n.cores.to.use){
  library(parallel)
  library(data.table)

  ## Load Data
  load(sprintf("data_prep/saved/for_analysis/JAGS.arr-%s.RData",         site.resolution), verbose=TRUE)
  load(sprintf("data_prep/saved/visitinfo-%s.RData",                     site.resolution), verbose=TRUE)
  load(sprintf("data_prep/saved/spatial.cleaned.data.filtered-%s.RData", site.resolution), verbose=TRUE)
  load(sprintf("data_prep/saved/climate-all-%s.RData",                   site.resolution), verbose=TRUE)

  floral.all <- readRDS(sprintf('data_prep/data/bee_suit_%s_floral_all.rds', site.resolution))
  
  message("Getting era and visit numbers into the climate summary based on the year")
  time.ref.dt <- data.table(unique(mysdata[,c("eraIDtext",
                                              "era",
                                              "visitIDtext",
                                              "visit",
                                              'year')]))

  ## Changes the class of the vector year to numeric 
  climate[,year:=as.numeric(year)]
  ## What's the column that I want to use to join these two things
  setkey(time.ref.dt,'year') 
  setkey(climate, 'year')
  ## Join the two, so now climate has era visit!
  climate <- climate[time.ref.dt]  

  message("Getting era into land use based on decade")
  time.ref.dt$decade <- paste0(substr(x     = time.ref.dt$year,
                                      start = 1,
                                      stop  = 3), '0')

  ## We'll call it decade instead becasue its in 10 year bins
  floral.all$decade <- floral.all$year 
  
  get.era <- function(de){
    use <- which(time.ref.dt$decade == de)
    if(length(use) == 0) return (NA)
    return(max(unique(time.ref.dt$era[use])))
  }
  
  ## Apply the get.era function to our data 
  floral.all$era <- mapply(get.era, de=floral.all$decade)

  ## Rename
  floral.all$V1 <- floral.all$floral_all  
  
  climate.means.tmax <- climate[variable == 'tmax',mean(V1, na.rm=TRUE), by=.(era,site)] #from data.table
  climate.means.prec <- climate[variable == 'prec',mean(V1, na.rm=TRUE), by=.(era,site)] #from data.table
  
  floral.all.means <- floral.all[, mean(V1, na.rm=TRUE), by=.(era,site)]
  
  era.means <- list(meanmaxt = climate.means.tmax, 
                    meanprec = climate.means.prec,
                    floral.all = floral.all.means#,
                    )
  save(era.means, file = sprintf("era.means-%s.RData", site.resolution))
  
  
  ## ~~~~~~~~~~~~~~ Adding variables to VISITINFO ~~~~~~~~~~~~~~~
  ## this is where we putvisit level info
  ## for example - average julian day
  ## we dont use this but it is avaliable if you want to add visit
  ## level variables
  
  message("Adding variables to visitinfo from mysdata")
  ## create a new object that is every combination of site and decade - regardless of if it was visited or not
  visitinfolong <- expand.grid(siteIDtext  = unique(visitinfo$siteIDtext), 
                               visitIDtext = unique(visitinfo$visitIDtext), 
                               eraIDtext   = unique(visitinfo$eraIDtext)
                               )
  
  ## and now add the other format of writing the values - just so we have it
  visitinfolong$site  <- as.numeric(sub('site.' ,'', visitinfolong$siteIDtext))
  visitinfolong$visit <- as.numeric(sub('visit.','', visitinfolong$visitIDtext))
  visitinfolong$era   <- as.numeric(sub('era.'  ,'', visitinfolong$eraIDtext))
  
  ## Add average julian day 
  get.avgjulianday <- function(site, era, visit){
    keep <- which(mysdata$site==site &
                  mysdata$era==era   &
                  mysdata$visit==visit)
    if (length(keep) == 0) return (NA)
    return(mean(mysdata$julianday[keep], na.rm=TRUE))
  }
  visitinfolong$avgjulianday <- mcmapply(get.avgjulianday,
                                         site     = visitinfolong$site,
                                         era      = visitinfolong$era,
                                         visit    = visitinfolong$visit,
                                         mc.cores = n.cores.to.use)
  summary(visitinfolong$avgjulianday)

  
  ## ~~~~~~~~~~~~~~ Adding variables to SITEINFO ~~~~~~~~~~~~~~~
  message("Adding variables to siteinfo from mysdata")
  
  siteinfo <- unique(mysdata[c("siteIDtext")])
  
  ## Function to add site size to siteinfo 
  get.site.area <- function(x){
    keep <- which(mysdata$siteIDtext == x)
    if (length(keep) == 0) return (NA)
    return(unique(mysdata$sitearea[keep]))
  }
  ## Apply that function
  siteinfo$sitearea <- mapply(get.site.area, x=siteinfo$siteID)

  
  ## ~~~~~~~~~~ Now putting siteinfo and visitinfolong into JAGS format ~~~~~~~~~~~~
  message("Putting siteinfo and visitinfo into array format")
  
  ## Making JAGS.site array 
  index <- match(rownames(JAGS.arr), siteinfo$siteIDtext)
  JAGS.site <- data.frame(siteIDtext = siteinfo$siteIDtext[index],
                          sitearea   = siteinfo$sitearea  [index]
                          )
  row.names(JAGS.site) <- JAGS.site[,'siteIDtext']
  
  ## Making JAGS.era array that we will use to run our model 
  ## Define what variables you want
  era.variables <- list(names(era.means))
  
  ##make an object with the right dimensions
  neravariables <- length(era.variables[[1]])
  JAGS.era <- array(NA, dim = c(site     = nsite, 
                                era      = nera,
                                variable = neravariables))
  
  ##name those dimensions
  dimnames(JAGS.era) <- c(dimnames(JAGS.arr)[c('site',
                                               'era')], 
                          variable = era.variables )

  ##then fill in JAGS.era with values from era.means based on the dimension names of JAGS.era 
  message("Filling in our era-level variables into an array (don't worry if this takes a little time)")
  for(ss in dimnames(JAGS.era)[['site']]){
    for (ee in dimnames(JAGS.era)[['era']]){
      for (var in dimnames(JAGS.era)[['variable']]){
        use <- which(era.means[[var]]$site == ss & paste0('era.',era.means[[var]]$era) == ee)
        JAGS.era[sprintf("%s",ss),  
                 sprintf("%s",ee), 
                 sprintf("%s",var)] <- as.numeric(era.means[[var]][use,'V1'])
      }
    }
  } 

  ## Making JAGS.visit array
  ## Define what variables you want - we dont use these, but the code
  ## is here if you do want to add visit-level variables, average
  ## julian day is an example of one you could use 
  visit.variables <- list(c("avgjulianday"))
  
  ## Make an object with the right dimensions
  nvariables <- length(visit.variables[[1]])
  JAGS.visit <- array(NA, dim=c(site     = nsite, 
                                visit    = nvisit, 
                                era      = nera,
                                variable = nvariables))
  
  ## Name those dimensions
  dimnames(JAGS.visit) <- c(dimnames(JAGS.arr)[c('site',
                                                 'visit',
                                                 'era')], 
                            variable = visit.variables ) 

  message("Filling in our visit-level variables into an array")
  ## Filling in according to the dimensions
  for(ss in dimnames(JAGS.visit)[['site']]){
    for (vv in dimnames(JAGS.visit)[['visit']]){
      for (ee in dimnames(JAGS.visit)[['era']]){
        for (var in dimnames(JAGS.visit)[['variable']]){
          use <- which(visitinfolong$siteIDtext == ss &  visitinfolong$visitIDtext == vv & visitinfolong$eraIDtext == ee)
          JAGS.visit[sprintf("%s",ss),  
                     sprintf("%s",vv), 
                     sprintf("%s",ee), 
                     sprintf("%s",var)] <- visitinfolong[use,sprintf('%s',var)]
        }
      }
    }
  }
  
  ## Vectors of that objct
  visit.vec <- as.numeric(sub("visit.",'', dimnames(JAGS.visit)[["visit"]]))
  era.vec   <- as.numeric(sub("era."  ,'', dimnames(JAGS.visit)[["era"]]))
  site.vec  <- as.numeric(sub("site." ,'', dimnames(JAGS.visit)[['site']]))

  ## Saving
  message("Saving objects")
  mysdata2 <- mysdata 
  save(mysdata2, file = "data_prep/saved/spatial.cleaned.data.2.RData")  
  
  save(siteinfo,      file = sprintf("data_prep/saved/siteinfo-%s.RData",      site.resolution))  
  save(visitinfolong, file = sprintf("data_prep/saved/visitinfolong-%s.RData", site.resolution))    
  
  save(JAGS.era, site.vec , era.vec,             file = sprintf("data_prep/saved/for_analysis/JAGS.era-%s.RData",   site.resolution))
  save(JAGS.visit, era.vec, visit.vec, site.vec, file = sprintf('data_prep/saved/for_analysis/JAGS.visit-%s.RData', site.resolution))
  save(JAGS.site, site.vec,                      file = sprintf("data_prep/saved/for_analysis/JAGS.site-%s.RData",  site.resolution))
  
}
