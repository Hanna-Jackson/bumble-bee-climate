## Written by: Hanna M. Jackson - hmj2@sfu.ca 

make.species.ranges <- function(site.resolution){

  ## Dependiecies
  library(maptools)
  library(parallel)
  library(raster)
  library(rgdal)
  library(rgeos)
  library(sp)
  library(spatstat)
  library(stringr)

  ## Load required files
  load(sprintf("data_prep/saved/dd.bee-%s.RData",               site.resolution), verbose=TRUE)  
  load(sprintf("data_prep/sites/%s/dd.sites.RData",             site.resolution), verbose=TRUE) 
  load(sprintf("data_prep/saved/for_analysis/JAGS.arr-%s.RData",site.resolution), verbose=TRUE)
  
  ## Ensure the CRS matchs
  dd.sites <- spTransform(dd.sites, CRS(proj4string(dd.bee))) 
  
  ## calculate and save ranges for all bee species
  tab.bees.all <- names(table(dd.bee$species))

  message("Drawing a convex hull around observations of each species, then get spatial object of just those sites")
  get.range <- function(bee) {
    dd.focal.bee <- dd.bee[dd.bee$species==bee,]

    region <- gConvexHull(dd.focal.bee)
    
    ## returns TRUE if that site intersects that region 
    keep <- !is.na(over(dd.sites, region)) 

    ## the subset of dd.sites that's in range
    dd.sites.in.range <- dd.sites[keep,] 

    ## so this thing I make when I apply this function is a list of subsets of dd.sites
    return(dd.sites.in.range) 
  }
  ranges <- sapply(X = tab.bees.all,
                   FUN = get.range) 
  

  ## Making the matrix that will constrain species to their ranges 
  get.species.site.ranges <- function(bee) {
    dd.focal.bee   <- dd.bee[dd.bee$species == bee,]
    region         <- gConvexHull(dd.focal.bee) 
    species.ranges <- !is.na(over(dd.sites, region))
    return(species.ranges) 
  }
  species.ranges <- sapply(X   = tab.bees.all,
                           FUN = get.species.site.ranges) 
  
  ##rename to site name 
  my_list <- strsplit(dimnames(species.ranges)[[1]] , "[.]")
  dimnames(species.ranges)[[1]] <- paste("site",
                                         formatC((lapply(my_list,`[[`, 1)),
                                                 width  = nchar(max(as.numeric(  lapply(my_list, `[[`, 1) ))),
                                                 format = "d",
                                                 flag   = "0"), # this just gets the site names and adds leading zeros 
                                         sep = '.')
  dimnames(species.ranges)[[2]] <- paste("species",
                                         dimnames(species.ranges)[[2]],
                                         sep = '.')

  ## Just checking to make sure that worked
  output   <- list()
  checking <- c(1:length(tab.bees.all))
  
  for(i in 1:length(tab.bees.all)){
    ## These are the sites that bee i is at
    ranges.list <- ranges[[i]][['site']]
    ## This is what we just made and we want to check the naming went correctly and didn't scramble them or something
    species.ranges.list <- dimnames(species.ranges)[[1]][which(species.ranges[,i])]
    ## Are these two the same?
    output[[i]] <- as.vector(ranges.list == species.ranges.list)
    
    ## If theyre not the same print the warning! 
    if(length(which(output[[i]] == F)) != 0){
      print(c("WARNING: for species number:", i, "species.range doesn't contain the true sites this species was found in")) # this lets us know which species messed up
      checking[i] <- 'bad'
    } else {
      checking[i] <- 'all good'
    }
  }
  
  if( length(which(checking == 'all good')) == length(checking) ){
    print("Success: species' ranges were assigned correctly")
  } else {
    print("Warning: at least one species' range was assigned incorrectly")
  }
  
  ## Rearragning species ranges to ensure they are in the same order as in JAGS.arr 
  index.species <- match(dimnames(JAGS.arr)[['species']],dimnames(species.ranges)[[2]])
  index.site    <- match(dimnames(JAGS.arr)[['site']],   dimnames(species.ranges)[[1]]) 

  ## Ensures that species and sites in JAGS.arr and species.ranges are in the same order
  species.ranges <- species.ranges[index.site, index.species] 

  ## Switch dimensions 
  species.ranges <- as.matrix(species.ranges)
  species.ranges <- aperm(a = species.ranges ,
                          perm = c(2,1))
  
  ## saving 
  save(ranges,         file=sprintf("data_prep/saved/ranges-%s.RData",                      site.resolution))
  save(species.ranges, file=sprintf("data_prep/saved/for_analysis/species.ranges-%s.RData", site.resolution))
}
