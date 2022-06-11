## Written by Hanna M. Jackson - hmj2@sfu.ca

workingdirectory <- '~/Dropbox/Bumble_Bee_Climate/climate-bumblebee/for_github'
setwd(workingdirectory)
message("working directory: ",getwd())

## Dependencies
library(maptools)
library(parallel)
library(raster)
library(rgdal)
library(rgeos)
library(sp)
library(spatstat)
library(RColorBrewer)
library(gridBase)
library(grid)

## Some plotting parameters and a list of labels
cex.main   <- 2.8
cex.axis   <- 1.9 
era.labels <- c("1900-1920",
                "1920-1940",
                "1940-1960",
                "1960-1980",
                "1980-2000",
                "2000-2020")

## Plotting parameters
my.par <- function(){
   par(mar = c(1, 1.5, 1, 0.1)) 
}

## Make a datastructure that will hold some values that we use in
## the manuscript
site.era.info <- array(NA, dim = c(6,3,3))
dimnames(site.era.info) <- list(c('era.1','era.2','era.3','era.4','era.5','era.6'),
                                c("nspecies",'nvisit','ndet'),
                                c("50","100","250"))
save(site.era.info, file = "site.era.info.RData")



## ~~~~~~~~~~~~~ Figure 1: Number of species in each site in each era ~~~~~~~~~~~~~~~~

make.species.figure <- function(site.resolution, cex.main, cex.axis){
  pdf(file = sprintf("~/Desktop/nspecies_site_era_%s.pdf",site.resolution), width = 18, height = 12) 
  layout(matrix(1:6, nrow = 2, ncol = 3, byrow = T))
  
  prj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  dd.box <- bbox2SP(n = 87,
                    e = -40,
                    s = 12,
                    w = -190,
                    proj4string = CRS(prj))
  dd.box <- bbox2SP(n = 4600000,
                    e = 3100000,
                    s = -3000000,
                    w = -5000000,
                    proj4string = crs(dd.AM))
  
  dd.AM <- raster::intersect(dd.AM, dd.box)
  
  ## This is how many VISITS for each species at each site
  JAGS.arr.nvisit <- apply(JAGS.arr,
                           c("site","era","species"),
                           sum,
                           na.rm = TRUE)
  
  ## And now we make it 0 or 1, was that species detectd at that site 
  JAGS.arr.nvisit[which(JAGS.arr.nvisit > 0, arr.ind = TRUE)] <- 1
  
  ## And now sum across species to get the number of species at each site this is the quantity that we will plot on the map!
  JAGS.arr.nspecies <- apply(JAGS.arr.nvisit,
                             c("site","era"),
                             sum,
                             na.rm = TRUE) 
  
  make.pannel <- function(era){
    ## Find the maximum across all eras so that each plot can be on the same scale
    max.nspecies <- max(JAGS.arr.nspecies)
    
    ## Then select only the data from the era that we want for this particular pannel
    JAGS.arr.era.nspecies <- JAGS.arr.nspecies[,era]
    
    ## Now we want to plot that on the map such that sites with higher numbers get shaded
    ## Step one is to match the names of the spatial object and the structure containing the data
    index <- match(dd.sites$site, names(JAGS.arr.era.nspecies))

    ## Make a spatial objext that has the number of species per site in that era 
    dd.sites@data[sprintf('era%snspecies',era)] <- JAGS.arr.era.nspecies[index]

    ## Decide which colours you're going to use!
    cols <- brewer.pal(9, 'Greens')

    ## Take the number of species in each site and turn it into a
    ## number from 1 to 9 (the number of colours) and store that
    ## number in a new object called era1cols, era2cols, ..., era9cols,
    dd.sites@data[[sprintf('era%scols',era)]] <-
      ceiling(dd.sites@data[[sprintf('era%snspecies',era)]] * 9 / max.nspecies)

    ## Subset to only the sites that have finite values (i.e. the ones
    ## we have data for)
    dd.sites.plot <- dd.sites[c(which(is.finite(dd.sites@data[[sprintf('era%snspecies',era)]]))),]
    my.par()

    ## Plot!
    plot(dd.AM,
         main     = NA,
         cex.main = cex.main,
         lwd      = 0.5,
         col      = "lightgrey")
    plot(dd.sites.plot,
         col    = cols[ dd.sites.plot@data[[sprintf('era%scols',era)]] ], 
         add    = TRUE,
         border = NA
         )
    era.labels <- c("1900-1920",
                    "1920-1940",
                    "1940-1960",
                    "1960-1980",
                    "1980-2000",
                    "2000-2020")
    title(main = sprintf("%s", era.labels[era]), 
          line = -1.3,
          cex.main = cex.main)

    ## Adding the inset histogram to the pannel 
    add.small.panel <- function(var, era2) {
      vp <- baseViewports()
      pushViewport(vp$inner,vp$figure,vp$plot)
      ## push viewport that will contain the inset
      pushViewport(viewport(x      = 0.07,
                            y      = 0.4,
                            width  = 0.43,
                            height = 0.3, 
                            just   = c("left","top")))
      
      op <- par(plt = gridPLT(), new = T)
      on.exit(par(op))
      ## plot inset figure
      ## We want a histogram of the values of each site, which we have
      ## stored in the object that we used to plot  
      hist(dd.sites.plot@data[[sprintf('era%s%s',era2,var)]],
           xlim = c(0,max.nspecies) ,
           main = NA,
           ylab = NA,
           xlab = NA,
           yaxt = 'n',
           xaxt = 'n',
           col = "lightgrey",
           breaks = seq(from=0, to=25, by=1))
      axis(side = 1, cex.axis = cex.axis)
      axis(side = 2, cex.axis = cex.axis, las=1)
      abline(v = mean(dd.sites.plot@data[[sprintf('era%s%s',era2,var)]]), col = "red")
      
      load("site.era.info.RData", verbose=T)
      site.era.info[paste0('era.',era),var,sprintf("%s",site.resolution)] <- mean(dd.sites.plot@data[[sprintf('era%s%s',era2,var)]])
      save(site.era.info, file = "site.era.info.RData")
      
      ## pop all viewports from stack
      popViewport(4)
    }
    
    ## Call that histogram function that we just made for the era that we're in! 
    add.small.panel(var = "nspecies",
                    era2 = era)

    ## Now we need one final plot for this pannel where we add the legend
    add.legend.panel <- function(){
      vp <- baseViewports()
      pushViewport(vp$inner,
                   vp$figure,
                   vp$plot)
      ## push viewport that will contain the inset
      pushViewport(viewport(x      = 0.07,
                            y      = 0.05,
                            width  = 0.43,
                            height = 0.04, 
                            just   = c("left","top")))
      op <- par(plt = gridPLT(),new = T)
      on.exit(par(op))
      
      ## plot the figure 
      image(matrix(1:9, nrow = 9), 
            col  = cols, 
            yaxt = "n",
            xaxt = 'n',
            yaxp = c(0,max.nspecies,9))
      
      ## pop all viewports from stack
      popViewport(4)
    }
    ## But we only want the legend on the first pannel, so we'll only
    ## actually call that function for the first era's pannel
    if(era == 1){
      add.legend.panel()
    }
    message("done ",era)
  }

  ## And now we call make.pannel to actually plot that all 
  for (ii in 1:6){
    make.pannel(era = ii)
  }

  ## Write to PDF 
  dev.off()
}






## ~~~~~~~~~~ Figure 2: Number of visits each site recieved in each era ~~~~~~~~~~~~~

make.visit.figure <- function(site.resolution,cex.main,cex.axis){
  ## PUT YOUR FILE PATH HERE! 
   pdf(file = sprintf("~/Desktop/nvisits_site_era_%s.pdf",site.resolution),
       width = 18,
       height = 12)
  layout(matrix(1:6,nrow=2,ncol=3,byrow=T))
  
   prj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
   dd.box <- bbox2SP(n = 87,
                     e = -40,
                     s = 12,
                     w = -190,
                     proj4string=CRS(prj))

   dd.box <- bbox2SP(n = 4600000,
                     e = 3100000,
                     s = -3000000,
                     w = -5000000,
                     proj4string = crs(dd.AM))
   
   dd.AM <- raster::intersect(dd.AM, dd.box)
   
  ## We're looking to get the number of visits that occured at each site in each era

  ## Sum across the detection history array across species 
  JAGS.arr.visit <- apply(JAGS.arr,
                          c("site","era","visit"),
                          sum,
                          na.rm=TRUE)
  ## Anything above 1 assign the value of 1 
  JAGS.arr.visit[which(JAGS.arr.visit>0,arr.ind = TRUE)] <- 1
  ## for each era, now each site visit has 1 or 0, visited or not
  
  ## So we're going to sum that across visits - this tells us the number of visits in each site-era
  JAGS.arr.visit <- apply(JAGS.arr.visit,c("site", "era"),
                          sum,
                          na.rm=TRUE) 

  ## Now we've made our data that we're going to use, lets actually
  ## make code to plot it!
   make.pannel <- function(era){
      
      JAGS.arr.visit <- JAGS.arr.visit[,era]
      ## Now we want to plot that on the map. Sites with higher numbers get shaded

     ## Match the sites in the spatial objext to the dimension names in our data
      index <- match(dd.sites$site, names(JAGS.arr.visit))

     ## Make a spatial object that contains the data 
      dd.sites@data[sprintf('era%snvisit',era)] <- JAGS.arr.visit[index]

     ## Choose colours 
     cols <- brewer.pal(5, 'Reds')

     ## Find the maximum value so we can make all plots on the same scale
     max.nvisit <- max(JAGS.arr.visit)

     ## Make an object that will store the colours 
      dd.sites@data[[sprintf('era%scols',era)]] <- dd.sites@data[[sprintf('era%snvisit',era)]] 

     ## Subset to only the sites that we have data for 
      dd.sites.plot <- dd.sites[c(which(is.finite(dd.sites@data[[sprintf('era%snvisit',era)]]))),]
      
      my.par()
      
      plot(dd.AM,
           main = NA,
           cex.main = cex.main,
           lwd = 0.5,
           col = "lightgrey")
      plot(dd.sites.plot,
           col = cols[ dd.sites.plot@data[[sprintf('era%scols',era)]]+1 ], 
           add = TRUE, 
           border = NA)
      title(main = sprintf("%s", era.labels[era]), 
            line = -1.3,
            cex.main = cex.main)

     ## Add the histogram
      add.small.panel <- function(var,era2) {
         vp <- baseViewports()
         pushViewport(vp$inner,vp$figure,vp$plot)
         
         ## push viewport that will contain the inset
         pushViewport(viewport(x = 0.11,
                               y = 0.4,
                               width = 0.22,
                               height = 0.3, 
                               just = c("left","top")))
         op <- par(plt=gridPLT(),new=T)
         on.exit(par(op))
         
         barplot(table(dd.sites.plot@data[[sprintf('era%s%s',era2,var)]]),
                 main = NA,
                 ylab = NA,
                 xlab = NA,
                 yaxt = 'n',
                 xaxt = 'n',
                 col = "lightgrey")
         
         axis(side   = 1,
              cex.axis = cex.axis,
              at     = c(0.7,1.9,3.1,4.3,5.5),
              labels = c("0","1","2","3","4"))
         
         axis(side = 2,
              cex.axis = cex.axis,
              las  = 1)

         load("site.era.info.RData", verbose = T)
         site.era.info[paste0('era.',era),var,sprintf("%s",site.resolution)] <- mean(dd.sites.plot@data[[sprintf('era%s%s',era2,var)]])
         save(site.era.info, file="site.era.info.RData")
         
         ## pop all viewports from stack
         popViewport(4)
      } 
      add.small.panel(var="nvisit",era2=era)

     ## Make the function that we will call to add the legend
      add.legend.panel <- function(){
         vp <- baseViewports()
         pushViewport(vp$inner,vp$figure,vp$plot)
         ## push viewport that will contain the inset
         pushViewport(viewport(x = 0.11,
                               y = 0.05,
                               width  = 0.22,
                               height = 0.04, 
                               just   = c("left","top")))
         op <- par(plt = gridPLT(),new = T)
         on.exit(par(op))
         ## plot inset figure
         image(matrix(1:5,nrow = 5), 
               col  = cols, 
               yaxt = "n",
               xaxt = 'n',
               yaxp = c(0,max.nvisit,9))
         
         ## pop all viewports from stack
         popViewport(4)
      } 
      if(era == 1){
         add.legend.panel()
      }
      
      message("done ",era)
   }
   make.pannel(era=1)
   make.pannel(era=2)
   make.pannel(era=3)
   make.pannel(era=4)
   make.pannel(era=5)
   make.pannel(era=6)
   
   dev.off()
}



## ~~~~~~~~~~ Figure 3: Number of detections in  each site in each era ~~~~~~~~~~~~~

make.detections.figure <- function(site.resolution,cex.main,cex.axis){
  pdf(file   = sprintf("~/Desktop/ndetections_site_era_%s.pdf", site.resolution),
      width  = 18,
      height = 12)
   layout(matrix(1:6,nrow = 2,ncol = 3,byrow = T))

   prj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
   dd.box <- bbox2SP(n = 87,
                     e = -40,
                     s = 12,
                     w = -190,
                     proj4string = CRS(prj))
   dd.box <- bbox2SP(n = 4600000,
                     e = 3100000,
                     s = -3000000,
                     w = -5000000,
                     proj4string = crs(dd.AM))
   
   dd.AM <- raster::intersect(dd.AM, dd.box)

  detarray <- apply(detarray,c("site","era","visit"),
                    sum,
                    na.rm = TRUE)
  
  JAGS.arr.det <- apply(detarray,c("site", "era"),
                        sum,
                        na.rm = TRUE)
  message("Number of detections by era for site resolution ", site.resolution)
  print(colSums(JAGS.arr.det))
  JAGS.arr.det <- log(JAGS.arr.det+1)

   make.pannel <- function(era){
      
      max.ndet <- max(JAGS.arr.det)
      JAGS.arr.det <- JAGS.arr.det[,era]
      ## Now we want to plot that on the map. Sites with higher numbers get shaded 
      
      index <- match(dd.sites$site, names(JAGS.arr.det))
      
      dd.sites@data[sprintf('era%sndet',era)] <- JAGS.arr.det[index] 
      
      cols <- brewer.pal(9, 'Blues')

      dd.sites@data[[sprintf('era%scols',era)]] <- ceiling(dd.sites@data[[sprintf('era%sndet',era)]] * 9 / max.ndet)
      
      dd.sites.plot <- dd.sites[c(which(is.finite(dd.sites@data[[sprintf('era%sndet',era)]]))),]
      
     my.par()
     
      ## Call the actual plot code: 
      plot(dd.AM,
           main = NA,
           cex.main = cex.main,
           lwd  = 0.5,
           col  = "lightgrey")
      
      plot(dd.sites.plot,
           col  = cols[ dd.sites.plot@data[[sprintf('era%scols',era)]]+1 ], 
           add  = TRUE, 
           border = NA)

       title(main = sprintf("%s", era.labels[era]), 
            line  = -1.3,
            cex.main = cex.main)
      
     add.small.panel <- function(var,era2) {
         vp <- baseViewports()
         pushViewport(vp$inner,vp$figure,vp$plot)
         ## push viewport that will contain the inset
         pushViewport(viewport(x = 0.07,
                               y = 0.4,
                               width  = 0.42,
                               height = 0.3, 
                               just   = c("left","top")))
         op <- par(plt = gridPLT(),new = T)
         on.exit(par(op))
         
         ## plot inset figure
         myvals <- dd.sites.plot@data[[sprintf('era%s%s',era2,var)]]
         if(site.resolution == 50){
           hist(myvals,
                main = NA,
                ylab = NA,
                xlab = NA,
                ylim = c(0,1400),
                xlim = c(0,7),
                yaxt = 'n',
                xaxt = 'n',
                col  = "lightgrey"
                )
           nice.vals <- c(0,5,50,500)

           ## convert nice vals to find locations
           tick.locs <- log(nice.vals+1)
           axis(1,cex.axis = cex.axis, at = tick.locs, labels = nice.vals)
         }

         if(site.resolution == 100){
           hist(myvals,
                main = NA,
                ylab = NA,
                xlab = NA,
                ylim = c(0,500),
                xlim = c(0,8),
                yaxt = 'n',
                xaxt = 'n',
                col = "lightgrey"
                )
           nice.vals <- c(0,5,50,500)
           ## convert nice vals to find locations
           tick.locs <- log(nice.vals+1)
           axis(1,cex.axis = cex.axis, at = tick.locs, labels = nice.vals)
         }
         
         if(site.resolution == 250){
           hist(myvals,
                main = NA,
                ylab = NA,
                xlab = NA,
                ylim = c(0,150),
                xlim = c(0,11),
                yaxt = 'n',
                xaxt = 'n',
                col  = "lightgrey"
                )
           nice.vals <- c(0,5,50,500,5000)
           ## convert nice vals to find locations
           tick.locs <- log(nice.vals+1)
           axis(1, cex.axis = cex.axis, at = tick.locs, labels = nice.vals)
         }
         axis(side = 2, cex.axis = cex.axis,las = 1)
   
         abline(v = mean(dd.sites.plot@data[[sprintf('era%s%s',era2,var)]]),col = "red")
         message("Mean is", mean(dd.sites.plot@data[[sprintf('era%s%s',era2,var)]]))
         
         load("site.era.info.RData",verbose=T)
         site.era.info[paste0('era.',era),var,sprintf("%s",site.resolution)] <- mean(dd.sites.plot@data[[sprintf('era%s%s',era2,var)]])
         save(site.era.info, file = "site.era.info.RData")

         ## pop all viewports from stack
         popViewport(4)
      } 
      
      add.small.panel(var = "ndet",era2 = era)
     
         
     add.legend.panel <- function(){
         vp <- baseViewports()
         pushViewport(vp$inner,vp$figure,vp$plot)
         ## push viewport that will contain the inset
         pushViewport(viewport(x = 0.07, y = 0.05, width = 0.42, height = 0.04,  just = c("left","top")))
         op <- par(plt = gridPLT(),new = T)
         on.exit(par(op))
         
         ## plot inset figure
         image(matrix(1:9, nrow = 9), 
               col  = cols, 
               yaxt = "n",
               xaxt = 'n',
               yaxp = c(0, max.ndet,9))
         
         ## pop all viewpornts from stack
         popViewport(4)
      } 
      
      if(era == 1){
         add.legend.panel()
      }

      message("done ",era)
   }
   
   make.pannel(era = 1)
   make.pannel(era = 2)
   make.pannel(era = 3)
   make.pannel(era = 4)
   make.pannel(era = 5)
   make.pannel(era = 6)
   
   dev.off()
}


for(ii in c(100)){
  message("Starting resolution: ",ii)
  
  load(sprintf("data_prep/saved/detection_array-%s.RData",            ii), verbose=TRUE)
  load(sprintf("data_prep/saved/dd.bee-%s.RData",                     ii), verbose=TRUE)  
  load(sprintf("data_prep/sites/%s/dd.sites.RData",                   ii), verbose=TRUE) 
  load(sprintf("data_prep/saved/for_analysis/JAGS.arr-%s.RData",      ii), verbose=TRUE)
  load(sprintf("data_prep/saved/for_analysis/species.ranges-%s.RData",ii), verbose=TRUE)
  load(        "data_prep/saved/dd.AM.RData",                                            verbose=TRUE)

  message("Making species figure")
  make.species.figure(site.resolution = ii,
                      cex.main        = cex.main,
                       cex.axis        = cex.axis)
  message("Making visit figure")
   make.visit.figure(site.resolution = ii,
                     cex.main        = cex.main,
                     cex.axis        = cex.axis)
  message("Making detections figure")
   make.detections.figure(site.resolution = ii,
                          cex.main        = cex.main,
                          cex.axis        = cex.axis)
}



## Some values for the manuscript
load("site.era.info.RData")

mean(exp(site.era.info[,"ndet","50"]))
mean(exp(site.era.info[,"ndet","100"]))
mean(exp(site.era.info[,"ndet","250"]))

mean(site.era.info[,"nvisit","50"])
mean(site.era.info[,"nvisit","100"])
mean(site.era.info[,"nvisit","250"])

mean(site.era.info[,"nspecies","50"])
mean(site.era.info[,"nspecies","100"])
mean(site.era.info[,"nspecies","250"])
