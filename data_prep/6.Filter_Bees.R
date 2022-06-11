## Written by: Hanna M. Jackson - hmj2@sfu.ca

filter.bees <- function(site.resolution){

  ## Load files
  load(sprintf("data_prep/saved/climate-all-%s.RData",                    site.resolution), verbose = TRUE)
  load(sprintf("data_prep/saved/spatial.cleaned.data-%s.RData",           site.resolution), verbose = TRUE)
  landuse <- readRDS(sprintf('data_prep/data/bee_suit_%s_floral_all.rds', site.resolution))

  ## Dependencies
  library(data.table)

  ## Rename climate$site to siteIDtext for consistency 
  colnames(climate)[which(colnames(climate) == "site")] <- "siteIDtext"
  
  ## And make one without the text
  climate$site <- as.numeric(sub("site.",
                                 "",
                                 climate$siteIDtext))
  
  ## Getting which sites have  NA or NAN in climate data
  message("Getting list of sites that have NA or NAN values in climate data")
  na.values <- which(!is.finite(climate$V1)) 

  ## Get the site names that those correspond to and provide a helpful message
  if (length(na.values) > 0){ 
    na.sites <- unique(climate$site[na.values]) 
    message(length(na.sites), " sites are have NA values in climate data")
  } else {
    message("No sites have climate data that is NA")
  }
  ## Note: That part did NOT filter this climate data, I just need the
  ## list of sites with NA for the next step where we instead filter
  ## the bees to only sites with associated climate data 

  ## Filtering the bee data by the climatae data 
  message("Removing observations from bee data that are in sites where any value of any climate or land use data is NA or NaN")
  remove <- vector()
  
  if(length(na.values) > 0){
    
    for(i in 1:length(na.sites)){
      remove <- c(remove, which(mysdata$site == na.sites[i]))
    }
    
    if(length(remove) > 0){
      mysdata <- mysdata[-remove,]
      message(length(remove), 'entries were removed from mysdata due to being from sites with NAs in CLIMATE data')
    }
    
  } else {
    message("Climate:No observations were removed from the bee data")
  }
  
  ## Filtering the bee data by the land use data 
  landuse.na.sites <- setdiff(paste0('site.',mysdata$site), unique(landuse$site))
  remove <- vector()
  
  if(length(landuse.na.sites) > 0){
    
    for(i in 1:length(landuse.na.sites)){
      remove <- c(remove, which(paste0('site.',mysdata$site) == landuse.na.sites[i]))
    }
    
    if(length(remove) > 0){
      mysdata <- mysdata[-remove,]
      message("Land Use: ",length(remove), ' entries were removed from mysdata due to being from sites with NAs in LAND USE data')
    } else {
      message("Land Use: Sites that had NAs didnt have bees so no entries were removed")
    }
    
  } else {
    message("Land Use: No observations were removed from the bee data")
  }
  
  ## Saving
  message("Saving mysdata")
  save(mysdata,file=sprintf("data_prep/saved/spatial.cleaned.data.filtered-%s.RData", site.resolution))
  ## this will be used to make JAGS.arr (observation history going into JAGS model) and it WON'T have sites that are NA for climate data now!
}
