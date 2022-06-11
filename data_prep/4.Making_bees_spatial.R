## Written by: Hanna M. Jackson - hmj2@sfu.ca

make.bees.spatial<-function(site.era.min,
                            remove.sites.without.certain.number.visits.in.any.era,
                            min.site.visits.in.any.era,
                            remove.sites.without.certain.number.observations.in.all.eras,
                            min.site.obs.in.all.eras,
                            site.resolution,
                            n.cores.to.use){

  ## Dependencies 
  library(maptools)
  library(parallel)
  library(raster)
  library(rgdal)
  library(rgeos)
  library(sp)
  library(spatstat)
  library(stringr)

  ## Load the necessary files
  load(file=        'data_prep/saved/cleaned.data.RData',                  verbose = TRUE) 
  load(file=sprintf('data_prep/sites/%s/dd.sites.RData', site.resolution), verbose = TRUE) 

  ## Making dd.bee - a spatial points data frame with all the bee records, and their details
  ## Select the coordinates
  coords  <- mycdata[,c('longitude', 'latitude')]
  ## The other things we want to go into dd.bee
  dd.attr <- mycdata[,c('year', 'visit', 'era', 'species')]
  ## Actually making it using those two objexts and the CRS we want
  dd.bee <- SpatialPointsDataFrame(coords = coords, 
                                   data = dd.attr,
                                   proj4string = CRS('+proj=longlat +datum=WGS84'))
  dd.bee <- spTransform(dd.bee, CRS(proj4string(dd.sites)))

  ## Now we have a spatial object of the bees, so next we want to know which of our sites each one is at
  message("Figuring out which site each bee is at")
  beesite <- over(dd.bee, dd.sites)
  save(beesite, file = sprintf('data_prep/saved/beesite-%s.RData',site.resolution))
  save(dd.bee,  file = sprintf("data_prep/saved/dd.bee-%s.RData",site.resolution))
  
  ## Add each bee's site to mycdata 
  message("Adding site to bee data (mycdata)")
  mycdata$siteIDtext <- beesite$site 
  mycdata$site <- sub(pattern     = "site.",
                      replacement = "",
                      x           = mycdata$siteIDtext)

  ## Make the notation consistent in the bee data 
  mycdata$visitIDtext   <- paste("visit",    
                                 formatC(mycdata$visit, width = nchar(max(mycdata$visit, na.rm=T)), format = "d", flag = "0"),
                                 sep = ".")
  mycdata$eraIDtext     <- paste("era",      
                                 formatC(mycdata$era, width = nchar(max(mycdata$era, na.rm=T)), format = "d", flag = "0"),
                                 sep = ".")
  mycdata$speciesIDtext <- paste("species", 
                                 mycdata$species, 
                                 sep = ".")
  
  ## ~~~~~~~~~~~~~~~ Filtering Sites based on bumble bee records ~~~~~~~~~~~~~~~~~~~
  
  ## Remove bees not in sites from mycdata
  message("Remove bees not found in any site")
  remove <- which(is.na(mycdata$site))
  if(!is.empty(remove)){ # if not empty (so nothing to remove) then remove them 
    mycdata <- mycdata[-remove,]
  }
  
  ## Add site areas to mycdata from site.areas
  message("Add site areas to mycdata")
  get.site.areas <- function(site){
    keep <- which(site.areas$site == site)
    if (length(keep) == 0) return (NA)
    return(unique(site.areas$area[keep]))
  }
  mycdata$sitearea <- mcmapply(get.site.areas,
                               site     = mycdata$siteIDtext,
                               mc.cores = n.cores.to.use)        

  ## Keep track of the number of records we'll remove 
  initial.nobs <- nrow(mycdata)
  
  ## Exclude sites that dont have enough observtaions of eras
  message("Excluding sites that dont have at least ", site.era.min, " eras observed")
  neras.per.site <- rowSums(table(unique(mycdata[,c("site","era")])))
  sites.keep     <- names(neras.per.site[which(neras.per.site >= site.era.min)])
  sites.remove   <- setdiff(unique(mycdata$site), sites.keep)
  mycdata        <- mycdata[which(mycdata$site %in% sites.keep),]
  message("nsites removed: ", length(sites.remove),
          ", nsites kept: ", length(sites.keep))

  ## Exclude sites that dont have enough repeat visits in any era 
  if(remove.sites.without.certain.number.visits.in.any.era == TRUE){
    message("Removing sites without repeat visits in ",min.site.visits.in.any.era, " eras")
    sitevisitsnew   <- unique(mycdata[,c("site","era","visit")])
    table.sitevisit <- table(sitevisitsnew)
    nvisit.site.era <- apply(table.sitevisit,c(1,2),sum)
    max.visit.site  <- apply(nvisit.site.era, 1, max)
    
    sites.keep   <- names(max.visit.site[which(max.visit.site >= min.site.visits.in.any.era)])
    sites.remove <- setdiff(unique(mycdata$site), sites.keep) #just so we know 
    mycdata      <- mycdata[which(mycdata$site %in% sites.keep),]
    message("nsites removed: ", length(sites.remove), ", nsites kept: ", length(sites.keep))
  }

  ## Exclude sites that dont have enough observations in all eras(from reviewer comment)
  if(remove.sites.without.certain.number.observations.in.all.eras == TRUE){
    message("Removing sites without ", min.site.obs.in.all.eras,
            " observations in every era")
    sites.remove <- vector()
    sites.keep   <- vector()
    for(ii in unique(mycdata$eraIDtext)){
      row.index    <- which(mycdata[,"eraIDtext"] == ii)
      site.era.obs <- data.frame(table(mycdata[row.index,'siteIDtext']))
      index.sites.remove.this.era <- which(site.era.obs$Freq < min.site.obs.in.all.eras)
      ## Here is the condition!
      sites.remove.this.era <- as.vector(site.era.obs[index.sites.remove.this.era,"Var1"])
      sites.remove          <- c(sites.remove,sites.remove.this.era)
    }
    sites.keep <- setdiff(unique(mycdata$siteIDtext), sites.remove)
    mytestdata <- mycdata[which(mycdata$siteIDtext %in% sites.keep),]
    message("Removed: ", length(sites.remove), " sites and kept: ", length(sites.keep))
    mycdata <- mytestdata
  }

  message("Went from ", initial.nobs, " to ", nrow(mycdata), " bumble bee observations due to removing sites with not enough records")
  
  ## Saving
  mysdata <- mycdata
  save(mysdata, file = sprintf("data_prep/saved/spatial.cleaned.data-%s.RData", site.resolution))
}
