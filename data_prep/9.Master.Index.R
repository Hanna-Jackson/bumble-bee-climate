## Written by: Hanna M. Jackson - hmj2@sfu.ca

what.to.iterate.over <- function(site.resolution){  

  ## Load files
  load(sprintf("data_prep/saved/for_analysis/JAGS.arr-%s.RData",      site.resolution), verbose=TRUE)
  load(sprintf("data_prep/saved/for_analysis/species.ranges-%s.RData",site.resolution), verbose=TRUE)

  ## Limit the modeling range of species to only their range
  sp.range <- species.ranges

  ## Limit analysis to visits that happened
  ## Sum across species so now 0s 1s 2s etc for this array of site x visit x era
  visit.arr <- apply(JAGS.arr,
                     c(1:3),
                     sum,
                     na.rm = TRUE)  

  ## Now figure out the indencies of the site,  species, era combinations we want to itterate over
  get.indices <- function(sp) {
    nsp.detected                 <- apply(JAGS.arr, 1:3, sum)  # sum across species so its (site, visit, era) - this tells us which sites were visited
    visit.arr[TRUE]              <- 1                          # set all values to 1 
    visit.arr[nsp.detected == 0] <- 0                          # then the values that are supposed to be 0 are set to 0 

    ## So now visit.arr tells us which site visits happened - this is an assumption!
    ## We're assuming that if species A was found at a site in a given
    ## time period, that that means that species B was looked for but
    ## not found and thus has a 0 at that site visit (as opposed to an NA)
    
    ## So now that we have our limits to visits and sites, we combine them! 
    ## This takes visit array, and for each species, if sp.range is FALSE for a given site, it assigns it a zero 
    visit.arr[!sp.range[sp,],,] <- 0       
    
    ## Next we want to find what the index of each of the 1s is 
    tmp <- which(visit.arr == 1, arr.ind=TRUE) # gives the array indicies ie [2,3,6] for all the 1s in visit.arr

    ## The OUTPUT of get.indicies:
    ## takes the thing we just made (tmp)and adds the species number as the first row, then it collumn
    ## binds that together with that same process for all the other species
    cbind(rep(sp,nrow(tmp)),
          tmp)  
  }
  
  master.index <- do.call(rbind,
                          lapply(X = 1:dim(JAGS.arr)[4],
                                 FUN = get.indices)
                          ) # apply get.indicies for each species 

  ## Giving it dimension names
  colnames(master.index) <- c('species','site','visit','era')
  rownames(master.index) <- 1:nrow(master.index)
  
  ## Ensure in the same order as in JAGS.arr for consistency 
  index        <- match(names(dimnames(JAGS.arr)),
                        colnames(master.index)) # take each dim of JAGS.arr and look for which column that is in master index 
  master.index <- master.index[,index] # now in the order of site visit era species 
  
  ## Making the vector versions 
  site    <- as.vector(master.index[,'site'])
  visit   <- as.vector(master.index[,'visit'])
  era     <- as.vector(master.index[,'era'])
  species <- as.vector(master.index[,'species'])
  
  visit.arr[visit.arr > 0] <- 1
  
  save(master.index, site, visit, era, species, file=sprintf("data_prep/saved/for_analysis/master.index-%s.RData", site.resolution))
  save(visit.arr,                               file=sprintf("data_prep/saved/for_analysis/visit.arr-%s.RData",    site.resolution))
  
}

