## Written by: Hanna M. Jackson - hmj2@sfu.ca

format.data <- function(visit.bin.time,
                        era.bin.time, 
                        start.year 
                        ){
  leifsdata <- read.csv(file = "data_prep/data/Leif_Richardson_BBNA_10-07-2020.csv")
  
  ## Define timesteps
  message("Defining a duration for eras and visits and adding it to bee data")

  ## For calculations later
  zero.year <- start.year-1
  
  ## Decade (No longer used)
  decade.bin.time  <- 10
  leifsdata$decade <- ceiling((leifsdata$year - zero.year+1) / decade.bin.time)
  
  ## Era
  leifsdata$era <- ceiling((leifsdata$year - zero.year) / era.bin.time)
  
  ## Visit
  visits.per.era     <- era.bin.time/visit.bin.time
  leifsdata$visitnum <- ceiling((leifsdata$year - zero.year) / visit.bin.time)
  leifsdata$visit    <- leifsdata$visitnum-(visits.per.era*(leifsdata$era-1))
  
  message("Creating and saving myfdata object (my formatted data)")
  ## Creating myfdata - choosing which variables will be included! 
  myfdata <- data.frame(
    ##Dataset
    source   = leifsdata$data.source,
    BBNAcode = leifsdata$BBNA.code,
    
    ##species 
    species     = leifsdata$species,
    detected.as = leifsdata$det.as,
    
    ##location
    latitude   = leifsdata$latitude,
    longitude  = leifsdata$longitude,
    continent  = leifsdata$continent,
    country    = leifsdata$country,
    state.prov = leifsdata$state.prov,
    county     = leifsdata$county,
    location   = leifsdata$location,
    siteinput  = leifsdata$site,
    
    ##observation
    julianday   = leifsdata$dayno,
    year        = leifsdata$year,
    decade      = leifsdata$decade,
    era         = leifsdata$era,
    visit       = leifsdata$visit,
    observers   = leifsdata$observers,
    det.by      = leifsdata$det.by,
    obs.notes   = leifsdata$obs.notes,
    geo.notes.1 = leifsdata$geo.notes.1,
    geo.notes.2 = leifsdata$geo.notes.2,
    notes.1     = leifsdata$notes.1,
    notes.2     = leifsdata$notes.2,
    notes.3     = leifsdata$notes.3,
    notes.4     = leifsdata$notes.4,
    notes.5     = leifsdata$notes.5
  )
  myfdata$beeID <- paste("bee",1:nrow(leifsdata), sep="")

  ## Saving so it can be loaded later
  save(myfdata ,  file = "data_prep/saved/formatted.data.RData")
  save(leifsdata, file = "data_prep/saved/leifsdata.RData")
}
