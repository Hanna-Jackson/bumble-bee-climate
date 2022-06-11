## Written by: Hanna M. Jackson - hmj2@sfu.ca

prep.climate <- function(site.resolution){
  ## Dependencies
  library(parallel)

  ## Load necessary files
  load(sprintf("data_prep/saved/spatial.cleaned.data-%s.RData", site.resolution),verbose = T)
  
  message("Loading climate data and getting it into long format")
  ## ~~~~~~~~ For the files that need fixing ~~~~~~~~~~
  
  ## IMPORTANT NOTE: When I downloaded my files from CHELSA I messed up on some and
  ## had to go back and get the names and then use those
  ## If you've downloaded the data yourseld and arent useing my
  ## version, you can delete this part and skip to the "files that
  ## dont need fixing" section
  ## Sorry for any confusion! 
  
  climate.files <- list.files(sprintf('data_prep/climate/need_year_month/resolution_%s', site.resolution))
  
  ## Create the structure
  climatefiles <- data.frame(climate.files = climate.files)
  
  for (i in 1:length(climatefiles$climate.files)){
    climatefiles$variable[i] <- sub(".rds",'', strsplit(x = climatefiles$climate.files, split = "-")[[i]][[3]])
    climatefiles$run_id[i]   <-     as.numeric(strsplit(x = climatefiles$climate.files, split = "-")[[i]][[2]])
  }
  
  ## Now a lookup for the correct values 
  load("data_prep/data/files.that.ran.RData", verbose=T) # this has the correct info
  
  prec <- data.frame(precfiles = files.that.ran[which(grepl(pattern="prec", x = files.that.ran))]) # in the order that they ran, so prec.files[1] was run_ID= 1 for prec 
  tmax <- data.frame(tmaxfiles = files.that.ran[which(grepl(pattern='tmax', x = files.that.ran))])
  
  prec$run_id <- 1:length(prec$precfiles)
  for (i in 1:length(prec$precfiles)){
    prec$year[i]  <- strsplit(x = prec$precfiles, split = "_")[[i]][[3]]
    prec$month[i] <- strsplit(x = prec$precfiles, split = "_")[[i]][[4]]
  }
  
  tmax$run_id <- 1:length(tmax$tmaxfiles)
  for (i in 1:length(tmax$tmaxfiles)){
    tmax$year[i]  <- strsplit(x = tmax$tmaxfiles, split = "_")[[i]][[3]]
    tmax$month[i] <- strsplit(x = tmax$tmaxfiles, split = "_")[[i]][[4]]
  }
  
  for (i in 1:nrow(climatefiles)){
    var    <- climatefiles$variable[i]
    run_id <- climatefiles$run_id[i]
    ## Which file name had that run ID for that variable 
    if (var == 'prec'){
      index <- which(prec$run_id == run_id)
      climatefiles$filename[i] <- prec$precfiles[index]
    } 
    if (var == 'tmax'){
      index <- which(tmax$run_id == run_id)
      climatefiles$filename[i] <- tmax$tmaxfiles[index]
    }
  }

  ## Now assign the year by getting it from the file name 
  for (i in 1:nrow(climatefiles)){
    climatefiles$year[i]  <- strsplit(climatefiles$filename,'_')[[i]][[3]]
    climatefiles$month[i] <- strsplit(climatefiles$filename,'_')[[i]][[4]]
  }
  
  ## Alright so we now have our `climatefiles` database, now what are we gonna do to it
  ## Well, for these files, unfortunately the month and year and variable needs to be adjusted - issue when downloaded them 
  message("Add data to climatefiles object from the data that needs to be updated")
  library(data.table)
  climate.list <- list()
  
  for (file in 1:length(climate.files)){ #for each file, load it,edit it, and attatch it to a main dataframe 
    ## Load
    climate.data <- readRDS(sprintf('data_prep/climate/need_year_month/resolution_%s/%s', site.resolution, climatefiles$climate.files[file]))
    climate.data <- data.table(climate.data)
    
    ## HERE WE EDIT IT - wont do this in the for loop for the good files
    climate.data$variable <- climatefiles$variable[file]
    climate.data$year     <- climatefiles$year[file]
    climate.data$month    <- climatefiles$month[file]
    
    ## Attatch new entry to the data frame
    climate.list[[file]] <- climate.data 
    message("climate files: ",file," out of ", length(climate.files))
  }
  length <- length(climate.list)
  
  ## ~~~~~~~~ For the files that don't need fixing ~~~~~~~~~~
  ## For the files that dont need adjusting, all we need is to loop through them and add them to the climate data 
  message("Add data to climatefiles object from the data that DONT need adjusting")
  more.climate.files <- list.files(sprintf('data_prep/climate/good_to_go/resolution_%s', site.resolution))
  
  for (file in (1:length(more.climate.files))){
    climate.data <- readRDS(sprintf('data_prep/climate/good_to_go/resolution_%s/%s', site.resolution, more.climate.files[file]))
    climate.data <- data.table(climate.data)
    climate.list[[file+length]] <- climate.data
    message("more climate files: ", file, " out of ", length(more.climate.files))
  }
  
  ## Now take the list of all climate data and put it into one dataframe 
  message("rbinding climate data together")
  climate <- rbindlist(climate.list)
  
  ## Saving climate data 
  save(climate, file=sprintf("data_prep/saved/climate-all-%s.RData", site.resolution))
}
