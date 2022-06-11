## Written by: Hanna M. Jackson - hmj2@sfu.ca

clean.data <- function(is.subset, 
                       subset.cap, 
                       end.year, 
                       start.year,
                       split.occi,
                       remove.problem.datasets){
  load(file="data_prep/saved/formatted.data.RData", verbose=TRUE)
  
  ## If we are subsetting the data, use subset.cap as number of data points
  if(is.subset == TRUE){
    myfdata <- myfdata[1:subset.cap,]
    message("Data has been subsetted to ", subset.cap, " entries. This feature is for troubleshooting only")
  }

  ## Remove duplicate observations by taking observations with unique entries of all these variables
  message("Removing duplicate observations")
  myfdata <- unique(myfdata[,c("source", 
                               'species', 
                               'year', 
                               'era',
                               'visit',
                               'julianday', 
                               'observers', 
                               'country', 
                               'county',
                               'siteinput',
                               'longitude',
                               'latitude')])
  
  ## Remove entries with no coordintes 
  message("Remove entries with no coordintes")
  myfdata <- myfdata[!is.na(myfdata$longitude) & !is.na(myfdata$latitude),] 
  
  ## Remove entries with no year
  message("Remove entries with no year")
  myfdata <- myfdata[!is.na(myfdata$year),]
  
  ## Remove entries with wrong years - we have no climate data beyond this period of time 
  message("Remove entries with year earlier than start.year = ",start.year)
  myfdata <- myfdata[!myfdata$year < start.year,]
  
  ## Remove entries after the end year (2020)
  message("Remove entries with year after end.year = ",end.year)
  myfdata <- myfdata[which(myfdata$year <= end.year),]
  

  ## Remove species with not enough data
  message("Removing species with not enough data")
  myfdata <- myfdata[!myfdata$species == 'cockerelli',]
  myfdata <- myfdata[!myfdata$species == 'kluanensis',]   
  myfdata <- myfdata[!myfdata$species == 'distinguendus',]  
  
  ## Split Occidentalis into two species: Occidentalis and Mckayi 
  if(split.occi==TRUE){
    mckayi.index <- which(myfdata$species == "occidentalis" & myfdata$latitude > 57)
    myfdata$latitude[mckayi.index]
    myfdata$species[mckayi.index]
    myfdata[mckayi.index,"species"] <- "mckayi"
    myfdata$latitude[which(myfdata$species == "mckayi")]
    message("Occidentalis split into Occi and mckayi")
  } else {
    message("Occidentalis not split")
  }
  
  ## Remove specific outlier observations as observed based on Species range plots
  message("Remove outlier observations based on extreme latitudes/longitudes")
  myfdata <- myfdata[-c(which(myfdata$species == "occidentalis" &  myfdata$longitude > -80),
                        which(myfdata$species == "appositus"    &  myfdata$latitude  <  30),
                        which(myfdata$species == "centralis"    &  myfdata$longitude > -95),
                        which(myfdata$species == "borealis"     &  myfdata$latitude  > 63.9      & myfdata$longitude > (-146)),
                        which(myfdata$species == "citrinus"     &  myfdata$longitude < (-115)),
                        which(myfdata$species == "griseocollis" &  myfdata$latitude  > 64),
                        which(myfdata$species == "polaris"      &  myfdata$longitude < (-100)    & myfdata$latitude < 53),
                        which(myfdata$species == "polaris"      &  myfdata$latitude  > 65        & myfdata$latitude < 75       & myfdata$longitude > (-21) ),
                        which(myfdata$species == "huntii"       &  myfdata$latitude  > 59 )
                        )
                    ,]
  
  if(remove.problem.datasets==TRUE){
    ## Remove observations from datasets with suspected issues
    index.laura <- which(grepl("Laura F"        , myfdata$source))
    index.SCAN  <- which(grepl("SCAN 10-04-2018", myfdata$source))
    myfdata     <- myfdata[-c(index.laura, index.SCAN),]
    message("Removed ",
            length(index.laura)+length(index.SCAN),
            "entries from laura F and SCAN datasets")
  } else(message("Didnt remove entries from Laura F and SCAN databases"))
  
  ## Making and saving 
  mycdata <- myfdata 
  save(mycdata, file="data_prep/saved/cleaned.data.RData") 
}
