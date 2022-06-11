## Written by: Hanna M. Jackson - hmj2@sfu.ca

make.occupancy.data <- function(site.resolution){ 
  
  load(file=sprintf("data_prep/saved/spatial.cleaned.data.filtered-%s.RData",
                    site.resolution),
       verbose=TRUE)
  
  ## This is a list of sites that had visits in the given possible 5 year 'visit' chunk of time 
  visitinfo <- unique(mysdata[c("siteIDtext",
                                "site",
                                "visitIDtext",
                                "visit",
                                "eraIDtext",
                                "era")])
  
  nsite    <- length(unique(mysdata$siteIDtext))
  nvisit   <- length(unique(mysdata$visitIDtext))
  nera     <- length(unique(mysdata$eraIDtext))
  nspecies <- length(unique(mysdata$speciesIDtext)) 

  ## make array of correct dimensions and empty fill with NAs
  detarray <- array(NA, dim = c(nsite,
                                nvisit,
                                nera,
                                nspecies))
  
  ## we'll make it a string instead of a number or it won't index properly  
  dimnames(detarray) <- list(site    = sort(unique(mysdata$siteIDtext)),
                             visit   = sort(unique(mysdata$visitIDtext)), 
                             era     = sort(unique(mysdata$eraIDtext)),
                             species = sort(unique(mysdata$speciesIDtext))
                             )
  
  ## Fill in site VISITS that occurred with a 0 
  for(i in 1:nrow(visitinfo)){
    si <- visitinfo$siteIDtext [i]
    vi <- visitinfo$visitIDtext[i]
    er <- visitinfo$eraIDtext  [i]
    detarray[si,vi,er,] <- 0 
  }
  
  ## Adding 1s to where DETECTIONS occurred 
  for(i in 1:nrow(mysdata)){
    si <- mysdata$siteIDtext   [i]
    vi <- mysdata$visitIDtext  [i]
    er <- mysdata$eraIDtext    [i]
    sp <- mysdata$speciesIDtext[i]
    ## This will keep adding one each time a bee was in that site
    detarray[si,vi,er,sp] <- detarray[si,vi,er,sp] + 1 
  }
  
  ## detarray should sum to the number of rows of mysdata - Quick Check 
  if ((sum(detarray, na.rm=TRUE)) == (length(mysdata$site) )) {
    message("SUCCESS: Occ Array sums to length of mysdata")
  } else {
    message("ERROR: Occ Array DOES NOT sum to length of mysdata")
  }
  
  ## detarray is now a count matrix, here we turn count matrix into 0 and 1 
  JAGS.arr               <- detarray
  JAGS.arr[JAGS.arr > 0] <- 1   

  ## saving
  save(detarray,                                file = sprintf("data_prep/saved/detection_array-%s.RData",       site.resolution))
  save(JAGS.arr, nsite, nvisit, nspecies, nera, file = sprintf("data_prep/saved/for_analysis/JAGS.arr-%s.RData", site.resolution))
  save(visitinfo,                               file = sprintf("data_prep/saved/visitinfo-%s.RData",             site.resolution))

}
