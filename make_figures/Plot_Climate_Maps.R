## Written by Hanna M. Jackson - hmj2@sfu.ca

## Set your working directory here:
workingdirectory <- 'PUT YOUR WORKING DIRECTORY HERE'
setwd(workingdirectory)
message("working directory: ",getwd())

## Dependencies
library(grid)
library(gridBase)
library(maptools)
library(parallel)
library(raster)
library(rgdal)
library(rgeos)
library(sp)
library(spatstat)
library(RColorBrewer)

## Make a function with all the plotting code that we will then call after

make.climate.figures <- function(site.resolution,
                                 cex.main,
                                 width,
                                 x,
                                 y,
                                 cex.axis,
                                 title.line){

  ## Put the file path for where you want the file to be saved here: 
  pdf(file   = sprintf("climate_fig%s.pdf",site.resolution),
      width  = 18,
      height = 6)

  ## Making it into a multipannel plot 
  layout(m <- matrix(1:3,nrow=1,byrow=T))

  ## Plotting parameters
  my.par <- function(){
    par(mar=c(0, 0.6, 0.4, 0.5),
        mai=c(0.3,0.8,0.3,0)) #oma is for each pannel and mar is for around all 4 pannels 
  }

  ## Loading required files
  load(sprintf("data_prep/sites/%s/dd.sites.RData",             site.resolution), verbose=TRUE) 
  load(        "data_prep/saved/dd.AM.RData",                                     verbose=TRUE)
  load(sprintf("data_prep/saved/for_analysis/JAGS.era-%s.RData",site.resolution), verbose=TRUE)

  ## Cropping the map 
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

  ## Print some calculated values - for the manuscript 
  JAGS.era.delta.t <- JAGS.era[,"era.6","meanmaxt"] - JAGS.era[,"era.1","meanmaxt"]
  message("Mean change in temp: ", mean(JAGS.era.delta.t,na.rm=T))
  
  JAGS.era.delta.p <- JAGS.era[,"era.6","meanprec"] - JAGS.era[,"era.1","meanprec"]
  message("Mean change in prec: ", mean(JAGS.era.delta.p,na.rm=T))
  
  JAGS.era.delta.lu <- JAGS.era[,"era.6","floral.all"] - JAGS.era[,"era.1","floral.all"]
  message("Mean change in floral: ", mean(JAGS.era.delta.lu,na.rm=T))

  ## Now we want to plot that on the map. Sites with higher numbers get shaded 
  index <- match(dd.sites$site, names(JAGS.era.delta.t))
  dd.sites@data["delta.t"]   <- JAGS.era.delta.t[index]
  dd.sites@data[,"delta.p"]  <- JAGS.era.delta.p[index]
  dd.sites@data[,"delta.lu"] <- JAGS.era.delta.lu[index]
  
  ## Function code to make an inset histogram!  - We will call this
      ## for each pannel we plot
  add.small.panel <- function(var,extra.width=0) {
    vp <- baseViewports()
    pushViewport(vp$inner,vp$figure,vp$plot)
    pushViewport(viewport(x      = x,
                          y      = 0.42,
                          width  = width+extra.width,
                          height = 0.3, 
                          just   = c("left","top")
                          )
                 )
    op <- par(plt=gridPLT(),new=T)
    on.exit(par(op))
    if(var == "delta.t"){
      hist(JAGS.era.delta.t,
           xlim = c(-1 * range(JAGS.era.delta.t)[2], range(JAGS.era.delta.t)[2]),
           main = NA,
           ylab = NA,
           xlab = NA,
           yaxt = 'n',
           xaxt = 'n',
           col  = "lightgrey",
           breaks = seq(from = -1 * range(JAGS.era.delta.t)[2],
                        to   = range(JAGS.era.delta.t)[2],
                        by   = (range(JAGS.era.delta.t)[2]--1 * range(JAGS.era.delta.t)[2])/20
                        )
           )
      
      axis(side = 1,cex.axis = cex.axis)
      axis(side = 2,cex.axis = cex.axis,las = 1)
      
      abline(v = mean(JAGS.era.delta.t),col="red")
    }
    if(var == "delta.p"){
      hist(JAGS.era.delta.p,
           xlim = c(range(JAGS.era.delta.p)[1], range(JAGS.era.delta.p)[1]*-1),
           main = NA,
           ylab = NA,
           xlab = NA,
           yaxt = 'n',
           xaxt = 'n',
           col  = "lightgrey",
           breaks=seq(from = range(JAGS.era.delta.p)[1],
                      to   = range(JAGS.era.delta.p)[1]*-1,
                      by   = ((range(JAGS.era.delta.p)[1]*-1)-(range(JAGS.era.delta.p)[1]))/20
                      )
           )

      axis(side = 1,cex.axis = cex.axis)
      axis(side = 2,cex.axis = cex.axis,las = 1)         
      abline(v   = mean(JAGS.era.delta.p),
             col = "red")
    }
    if(var == "delta.lu"){
      hist(JAGS.era.delta.lu,
           xlim = c(range(JAGS.era.delta.lu)[1], range(JAGS.era.delta.lu)[1]*-1),
           main = NA,
           ylab = NA,
           xlab = NA,
           yaxt = 'n',
           xaxt = 'n',
           col  = "lightgrey",
           breaks=seq(from = range(JAGS.era.delta.lu)[1],
                      to   = range(JAGS.era.delta.lu)[1]*-1,
                      by   = ((range(JAGS.era.delta.lu)[1]*-1)-(range(JAGS.era.delta.lu)[1]))/20
                      )
           )
      axis(side = 1,cex.axis = cex.axis)
      axis(side = 2,cex.axis = cex.axis,las = 1)         
      abline(v = mean(JAGS.era.delta.lu),
             col="red")
    }
    ## pop all viewports from stack
    popViewport(4)
  } 
  
  
  ## ~~~~~~~~ PANNEL 1: Temperature ~~~~~~~~~~~ 
  lower <- range(dd.sites@data[['delta.t']],na.rm = T)[[1]]
  upper <- range(dd.sites@data[['delta.t']],na.rm = T)[[2]]
  real.range.t <- upper-lower
  my.range.t   <- 2*upper
  my.lower.t   <- upper*-1
  
  cols  <- rev(brewer.pal(11,"RdBu"))
  ncols <- length(cols)
  
  bin.size <- my.range.t/ncols
  
  dd.sites.plot <- dd.sites[c(which(is.finite(dd.sites@data[['delta.t']]))),]
  
  ## Create new data column
  dd.sites.plot@data[,'delta.t.cols'] <- rep(NA ,nrow(dd.sites.plot@data))

  ## Add colours based on values
  allvals <- vector()
  for (i in 0:(ncols-1)){
    values <- which(dd.sites.plot@data[['delta.t']]  > (my.lower.t+bin.size*i) & dd.sites.plot@data['delta.t'] <= (my.lower.t+bin.size*(i+1)))
    dd.sites.plot@data[values,'delta.t.cols']<-cols[i+1] 
    allvals <- c(values,allvals)
  }
  
  my.par()
  ## Plot grey everywhere   
  plot(dd.AM,
       main     = NA,
       cex.main = cex.main,
       lwd      = 0.5,
       col      = "lightgrey")
  ## Then plot our sites with their corresponding correct colour on top
  plot(dd.sites.plot,
       col    = dd.sites.plot@data[['delta.t.cols']], 
       add    = TRUE, 
       border = NA)
  title(main     = "Temperature Change", 
        line     = title.line,
        cex.main =cex.main)
  
  add.small.panel(var = "delta.t",extra.width = 0.01) # Call the function that we made earlier
  
  add.legend.panel <- function(){
    vp <- baseViewports()
    pushViewport(vp$inner,vp$figure,vp$plot)
    pushViewport(viewport(x      = x,
                          y      = y,
                          width  = width,
                          height = 0.04, 
                          just   = c("left","top")))
    op <- par(plt = gridPLT(),new = T)
    on.exit(par(op))
    
    image(matrix(1:9, nrow = 9), 
          col  = rev(brewer.pal(11,"RdBu")), 
          yaxt = "n",
          xaxt = 'n',
          yaxp = c(my.lower.t,upper,9)
          ) 
    
    ## pop all viewports from stack
    popViewport(4)
  } 
  add.legend.panel()
  
  

  
  ## ~~~~~~~~ PANNEL 2: Precipitation ~~~~~~~~~~~ 
  range(dd.sites@data[['delta.p']],na.rm = T)
  lower <- range(dd.sites@data[['delta.p']],na.rm = T)[[1]]
  upper <- range(dd.sites@data[['delta.p']],na.rm = T)[[2]]
  real.range.p <- upper-lower
  my.range.p <- 2*abs(lower)
  my.upper.p <- lower*-1
  my.lower.p <- lower
  
  cols <- brewer.pal(11,"PuOr")
  
  ncols <- length(cols)
  
  bin.size <- my.range.p/ncols
 
  dd.sites.plot <- dd.sites[c(which(is.finite(dd.sites@data[['delta.p']]))),]
  
  ## Create new data column
  dd.sites.plot@data[,'delta.p.cols'] <- rep(NA,nrow(dd.sites.plot@data))
  
  ## Add colours based on values
  allvals <- vector()
  for (i in 0:(ncols-1)){
    values <- which(dd.sites.plot@data[['delta.p']]  >= (my.lower.p+bin.size*i) & dd.sites.plot@data['delta.p'] <= (my.lower.p+bin.size*(i+1)))
    dd.sites.plot@data[values,'delta.p.cols'] <- cols[i+1] 
    allvals <- c(values,allvals)
  }
  my.par()
  
  plot(dd.AM,
       main     = NA,
       cex.main = cex.main,
       lwd      = 0.5,
       col      = "lightgrey")
  
  plot(dd.sites.plot,
       col    = dd.sites.plot@data[['delta.p.cols']], 
       add    = TRUE, 
       border = NA)
  
  title(main = "Precipitation Change", 
        line = title.line,
        cex.main = cex.main)
  
  add.small.panel(var="delta.p")
  
  add.legend.panel <- function(){
    vp <- baseViewports()
    pushViewport(vp$inner,vp$figure,vp$plot)
    pushViewport(viewport(x      = x,
                          y      = y,
                          width  = width,
                          height = 0.04, 
                          just   = c("left","top")))
    op <- par(plt = gridPLT(),new = T)
    on.exit(par(op))
    ## plot inset figure - The legend for the histogram 
    image(matrix(1:9,nrow=9), 
          col  = brewer.pal(11,"PuOr"), 
          yaxt = "n",
          xaxt = 'n',
          yaxp = c(my.lower.t,upper,9)
          ) 
    
    ## pop all viewports from stack
    popViewport(4)
  } 
  add.legend.panel()


  
  
  ## ~~~~~~~~ PANNEL 3: Land Use ~~~~~~~~~~~ 
  range(dd.sites@data[['delta.lu']],na.rm=T)
  lower <- range(dd.sites@data[['delta.lu']],na.rm=T)[[1]]
  upper <- range(dd.sites@data[['delta.lu']],na.rm=T)[[2]]
  real.range.p <- upper-lower
  my.range.lu  <- 2*abs(lower)
  my.upper.lu  <- lower*-1
  my.lower.lu  <- lower
  
  cols <- c("#BE2758", "#D56680", "#E899A9", "#F6CCD3", "#FFFFFF", "#D0D5F3", "#A1ACE6", "#7085DA","#2F5DD4")
  ncols <- length(cols)
  
  bin.size <- my.range.lu/ncols
  
  dd.sites.plot <- dd.sites[c(which(is.finite(dd.sites@data[['delta.lu']]))),]
  
  ## Create new data column
  dd.sites.plot@data[,'delta.lu.cols'] <- rep(NA ,nrow(dd.sites.plot@data))
  ## Add colours based on values
  allvals <- vector()
  for (i in 0:(ncols-1)){
    values <- which(dd.sites.plot@data[['delta.lu']]  >= (my.lower.lu+bin.size*i) & dd.sites.plot@data['delta.lu'] <= (my.lower.lu+bin.size*(i+1)))
    dd.sites.plot@data[values,'delta.lu.cols'] <- cols[i+1] 
    allvals <- c(values,allvals)
  }
  
  my.par()
  
  plot(dd.AM,
       main     = NA,
       cex.main = cex.main,
       lwd      = 0.5,
       col      = "lightgrey")
  
  plot(dd.sites.plot,
       col    = dd.sites.plot@data[['delta.lu.cols']], 
       add    = TRUE, 
       border = NA)
  
  title(main="Floral Resources Change", 
        line=title.line,
        cex.main=cex.main)
  
  add.small.panel(var="delta.lu")
  
  add.legend.panel <- function(){
    vp <- baseViewports()
    pushViewport(vp$inner,vp$figure,vp$plot)
    pushViewport(viewport(x      = x,
                          y      = y,
                          width  = width,
                          height = 0.04, 
                          just   = c("left","top")))
    op <- par(plt=gridPLT(),new=T)
    on.exit(par(op))
    ## plot inset figure
    image(matrix(1:9, nrow=9), 
          col  = cols,
          yaxt = "n",
          xaxt = 'n',
          yaxp = c(my.lower.t,upper,9)
          )
    
    ## pop all viewports from stack
    popViewport(4)
  } 
  add.legend.panel()
    
  dev.off()
}

## This is the code that actually makes the figure!
## It calls the huge function that we just made above
## Put in the site resolutions that you want to make the figure for in
## the first line here
for (i in c(50, 100, 250)){
  make.climate.figures(site.resolution = i,
                       cex.main        = 2.2,
                       width           = 0.51, 
                       x               = -0.065,
                       y               = 0.040, 
                       cex.axis        = 1.7,
                       title.line      = -1.2
                       )
}


