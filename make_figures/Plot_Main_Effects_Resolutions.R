## Written by Hanna M. Jackson - hmj2@sfu.ca

## Put your working directory here! 
workingdirectory <- 'PUT YOUR WORKING DIRECTORY HERE'
setwd(workingdirectory)
message("working directory: ",getwd())

## Specify which resolutions and what their colours and line types will be! 
resolutions <- c(50, 100, 250)
cols        <- c("black","black","black")
linetypes   <- c("dotted","dashed","solid")

## Function that makes a colour transparent - we will call this later
makeTransparent <- function(..., alpha=0.5) {
  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
  alpha = floor(255*alpha)  
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  apply(newColor, 2, .makeTransparent, alpha=alpha)
}

## We will make a huge function that contains all the code needed to
## make this plot. Then once that function is loaded, we will call it
## with the options we want to actually plot (see code for calling it
## at the end of the code)

make.resolutions.figure <- function(full=FALSE, cex.lab, cex.axis, cex.main){
  ## Your file path for the figure you want to make here 
  pdf(file   = sprintf("MainEffects_resolutions.pdf"),
      width  = 9,
      height = 9)
  
  m <- matrix(c(1,2,3,4), ncol=2, nrow=2,byrow=TRUE)
  layout(m)
  
  library(abind)
  
  par(oma = c(0.2, 0.2,   0, 0.2), 
      mar = c(4.5,   4, 0.1, 0.1),
      mgp = c(2,0.2,0),
      tcl = 0, cex.axis = cex.axis, cex.main = cex.main, pty = 's')
  
  ## ~~~~~~~~~~~ Initiation - Loading files and sourcing ~~~~~~~~~~~
  ## Some required functions 
  expit <- function(x) 1/(1+exp(-x))
  logit <- function(x) log(x/(1-x))  

  covariate <- list(SR50  = 'output/modeloutput-env-50.RData',
                    SR100 = 'output/modeloutput-env-100.RData' ,
                    SR250 = 'output/modeloutput-env-250.RData' )
  era <- list(SR50  =       'output/modeloutput-era-50.RData',
              SR100 =       'output/modeloutput-era-100.RData' ,
              SR250 =       'output/modeloutput-era-250.RData' )
  files <- list(era = era,
                covariate = covariate)
  
  ## Load the outputs from both the era and covariate model for all 3 resolutions
  load(files[["era"]][["SR50"]])
  era50 <- bugs$BUGSoutput$sims.matrix
  
  load(files[["era"]][["SR100"]])
  era100 <- bugs$BUGSoutput$sims.matrix
  
  load(files[["era"]][["SR250"]])
  era250 <- bugs$BUGSoutput$sims.matrix
  
  load(files[["covariate"]][["SR50"]])
  covariate50 <- bugs$BUGSoutput$sims.matrix
  
  load(files[["covariate"]][["SR100"]])
  covariate100 <- bugs$BUGSoutput$sims.matrix
  
  load(files[["covariate"]][["SR250"]])
  covariate250 <- bugs$BUGSoutput$sims.matrix

  ## Now put those all into a list so we can easily access any of
  ## those model outputs
  chains.list <- list(era50  = era50,
                      era100 = era100,
                      era250 = era250,
                      covariate50  = covariate50,
                      covariate100 = covariate100,
                      covariate250 = covariate250)
  
  load("data_prep/saved/for_analysis/JAGS.era-50.RData",verbose=T)
  JAGS.era.unscaled <- JAGS.era

  ## Function that we will use to make pannel labels 
  fig.labs <- c("a","b","c","d")
  add.fig.label <- function(fig.num,extra){
    text(x = par('usr')[2]-par('usr')[2]*(0.08)+extra, 
         y = par('usr')[4]*0.95, 
         fig.labs[fig.num],
         cex = 2, 
         pos = 4
         )
  }

  
#####~~~~~~~~~~~ PANNEL 1: Occupancy vs Era ~~~~~~~~~~~~~#####
  for (ii in 1:length(resolutions)){
    message("Loop Number ",ii)
    sims.mat <- chains.list[[sprintf("era%s",resolutions[ii])]]
    ## figure out x-values for our figure
    x.vals <- seq(from   = 1, 
                  to     = 6,
                  length = 6)

    ## function to extract the chains of y-values for a given species
    get.for.sp <- function(ss) {
      get.y.val <- function(xx) {
        range <- my.data$species.ranges[ss,]
        chains <- sims.mat[,'psi.0'] +
          sims.mat[,sprintf('psi.sp[%d]', ss)]                 +       
          sims.mat[,sprintf('psi.era[%d]', ss)]  * xx         +
          sims.mat[,        'psi.sitearea']      * mean(my.data$sitearea)         
      }
      sapply(x.vals, get.y.val)
    }
    message("Extracting chains for era, site resolution: ", resolutions[ii])    ## extract chains by species
    chains.by.sp <- lapply(1:my.data$nspecies, get.for.sp)
    ## and then stitch together into a multidimensional array
    chains.arr <- abind(chains.by.sp, along=3)
                                        # browser()
    ## Set up this pannel on the first itteration
    if(ii==1){
      plot(NA,
           xlim = range(x.vals),
           ylim = c(0,1),
           xlab = '',
           ylab = 'Occupancy Probability',
           las  = 1,
           xaxt = "n",
           cex.lab  = cex.lab,
           cex.axis = cex.axis)
    }

    add.line <- function(xx, bci=FALSE, cc, lwd=1, ...) {
      to.plot <- apply(xx, 2, function(ii)
        expit(c(mean = mean(ii), quantile(ii, probs = c(0.025,0.975)))))
      ## add mean line
      lines(x = x.vals, y = to.plot["mean",], type = "l", col = cc,lwd = lwd,lty = linetypes[ii])
      ## add BCI, if desired
      if(bci) {
        polygon(x   = c(x.vals, rev(x.vals)),
                y   = c(to.plot["2.5%",],rev(to.plot["97.5%",])),
                col = makeTransparent(sprintf("%s",cc), alpha=0.25),
                border = NA)
      }
    }

    ## add the mean line across species
    xx <- apply(chains.arr, 1:2, mean)
    add.line(xx, bci = TRUE, cc = cols[ii], lwd = 2, lty = 2)
    
    ## After plotting the last line on this pannel, run some final plotting options 
    if(ii == length(resolutions)){
      op <- par(new=TRUE) #does this command, but returns the previous par
      labs <- c("1900-1920",
                "1920-1940",
                "1940-1960",
                "1960-1980",
                "1980-2000",
                "2000-2020")
      text(x   = seq_along(labs)+0.1, 
           y   = par('usr')[3]-0.01,
           labels = labs,
           srt = 20,
           adj = c(1.1,1.1),
           xpd = NA,
           cex = 0.9)
      
      axis(tick = T,side = 1,labels = F)

      mtext("Era",
            side = 1,
            line = 3,
            cex  = 1.25)
      
      par(op) #calling the previous par, not the new=T one 
      add.fig.label(fig.num = 1,extra = 0)

      legend(x = 1,
             y = 0.98,
             legend = c("50 x 50","100 x 100","250 x 250"),
             title  = expression(paste("Site Size (km)")),
             col    = cols, 
             pt.cex = 2,
             cex  = 1.4,
             lty  = linetypes,
             lwd  = 2,
             bty  = "n",
             ncol = 1)
    }
  }

  

##### ~~~~~~~~~~~~~ PANNEL 2: Occupancy vs Temperature  ~~~~~~~~~~~~~~~~ ##### 
  for (ii in 1:length(resolutions)){
    sims.mat <- chains.list[[sprintf("covariate%s",resolutions[ii])]]

    ## figure out x-values for our figure
    x.vals <- seq(from = min(my.data$meanmaxt),
                  to   = max(my.data$meanmaxt),
                  length.out = 101)

    ## function to extract the chains of y-values for a given species
    get.for.sp <- function(ss) {
      get.y.val <- function(xx) {
        range <- my.data$species.ranges[ss,]
        chains <- sims.mat[,'psi.0'] +
          sims.mat[,sprintf('psi.sp[%s]', ss)] +
          sims.mat[,sprintf('psi.meanmaxt[%s]', ss)] * xx +
          sims.mat[,        'psi.meanmaxt.sq']       * xx * xx +
          sims.mat[,sprintf('psi.meanprec[%s]', ss)] * mean(my.data$meanprec[range,]) +
          sims.mat[,sprintf('psi.floral[%s]', ss)]   * mean(my.data$floral[range,]) +
          sims.mat[,        'psi.sitearea']          * mean(my.data$sitearea[range])
      } 
      sapply(x.vals, get.y.val)
    }
    
    message("Extracting chains for temperature, site resolution: ", resolutions[ii])
    ## extract chains by species
    chains.by.sp <- lapply(1:my.data$nspecies, get.for.sp)
    ## and then stitch together into a multidimensional array
    chains.arr <- abind(chains.by.sp, along=3)
    ## now we can use this to make our plot
    if(ii==1){
      plot(NA,
           xlim = range(x.vals),
           ylim = c(0,1),
           xlab = 'Temperature (°C)',
           ylab = NA,
           las  = 1,
           xaxt = "n",
           cex.lab  = cex.lab,
           cex.axis = cex.axis)
    }
    
    add.line <- function(xx, bci=FALSE, cc, lwd=1, ...) {
      to.plot <- apply(xx, 2, function(ii)
        expit(c(mean = mean(ii), quantile(ii, probs = c(0.025,0.975)))))
      ## add mean line
      lines(x = x.vals, y = to.plot["mean",], type = "l", col = cc, lwd = lwd,lty = linetypes[ii])
      ## add BCI, if desired
      if(bci) {
        polygon(x   = c(x.vals, rev(x.vals)),
                y   = c(to.plot["2.5%",],rev(to.plot["97.5%",])),
                col = makeTransparent(sprintf("%s",cc), alpha = 0.25),
                border = NA)
      }
    }
    ## add the mean line across species
    xx <- apply(chains.arr, 1:2, mean)
    add.line(xx, bci = TRUE, cc = cols[ii], lwd = 2, lty = 2)

    if(ii == length(resolutions)){
      ## Final plot arguments 
      op <- par(new = TRUE) #does this command, but returns the previous par
      plot(NA,
           xlim = range(JAGS.era.unscaled[,,'meanmaxt']),
           ylim = c(0,1),
           xlab = '',
           ylab = '',
           las = 1)
      par(op) #calling the previous par, not the new=T one 
      add.fig.label(fig.num = 2,extra = 0)
    }
  }

  

##### ~~~~~~~~~~~~ PANNEL 3: Occupancy  vs Precipitation ~~~~~~~~~~~~~~~ ##### 
  for (ii in 1:length(resolutions)){
    sims.mat <- chains.list[[sprintf("covariate%s",resolutions[ii])]]
    
    ## figure out x-values for our figure
    x.vals <- seq(from = min(my.data$meanprec),
                  to   = max(my.data$meanprec),
                  length.out = 1001)

    ## function to extract the chains of y-values for a given species
    get.for.sp <- function(ss) {
      get.y.val <- function(xx) {
        range <- my.data$species.ranges[ss,]
        chains <- sims.mat[,'psi.0'] +
          sims.mat[,sprintf('psi.sp[%s]', ss)] +
          sims.mat[,sprintf('psi.meanmaxt[%s]', ss)]    * mean(my.data$meanmaxt[range,]) +
          sims.mat[,        'psi.meanmaxt.sq']          * mean(my.data$meanmaxt[range,]^2) +
          sims.mat[,sprintf('psi.meanprec[%s]', ss)]    * xx +
          sims.mat[,sprintf('psi.floral[%s]', ss)]      * mean(my.data$floral[range,]) +
          sims.mat[,        'psi.sitearea']             * mean(my.data$sitearea[range])
      }
      sapply(x.vals, get.y.val)
    }
    message("Extracting chains for precipitation, site resolution: ", resolutions[ii])
    ## extract chains by species
    chains.by.sp <- lapply(1:my.data$nspecies, get.for.sp)
    ## and then stitch together into a multidimensional array
    chains.arr <- abind(chains.by.sp, along=3)
    ## now we can use this to make our plot

    if(ii == 1){
      plot(NA,
           xlim = range(x.vals),
           ylim = c(0,1),
           xlab = 'Precipitation (kg/m²)',
           ylab = 'Occupancy Probability',
           las  = 1,
           xaxt = "n",
           cex.lab = cex.lab,
           cex.axis= cex.axis)
    }
    add.line <- function(xx, bci=FALSE, cc, lwd=1, ...) {
      to.plot <- apply(xx, 2, function(ii)
        expit(c(mean = mean(ii), quantile(ii, probs = c(0.025,0.975)))))
      ## add mean line
      lines(x = x.vals, y = to.plot["mean",], type = "l", col = cc, lwd = lwd,lty = linetypes[ii])
      ## add BCI, if desired
      if(bci) {
        polygon(x   = c(x.vals, rev(x.vals)),
                y   = c(to.plot["2.5%",],rev(to.plot["97.5%",])),
                col = makeTransparent(sprintf("%s",cc), alpha=0.25),
                border = NA)
      }
    }
    ## add the mean line across species
    xx <- apply(chains.arr, 1:2, mean)
    add.line(xx, bci = TRUE, cc = cols[ii], lwd = 2, lty = 2)

    if(ii == length(resolutions)){
      ## Final plot arguments 
      op <-par(new = TRUE) #does this command, but returns the previous par
      plot(NA,
           xlim = range(JAGS.era.unscaled[,,'meanprec']), 
           ylim = c(0,1),
           xlab = '', 
           ylab = '',
           las  = 1)
      par(op) #calling the previous par, not the new=T one 
      
      add.fig.label(fig.num = 3, extra = 0)
    }
  }

  
##### ~~~~~~~~~~~~~~~~~ PANNEL 4: Occupancy vs Floral Resources  ~~~~~~~~~~~~~~~~~~~~~~ ##### 
  for (ii in 1:length(resolutions)){
    sims.mat <- chains.list[[sprintf("covariate%s",resolutions[ii])]]
    
    ## figure out x-values for our figure
    x.vals <- seq(from = min(my.data$floral),
                  to   = max(my.data$floral),
                  length.out = 1001)

    ## function to extract the chains of y-values for a given species
    get.for.sp <- function(ss) {
      get.y.val <- function(xx) {
        range <- my.data$species.ranges[ss,]
        chains <- sims.mat[,'psi.0'] +
          sims.mat[,sprintf('psi.sp[%s]', ss)] +
          sims.mat[,sprintf('psi.meanmaxt[%s]', ss)]    * mean(my.data$meanmaxt[range,]) +
          sims.mat[,        'psi.meanmaxt.sq']          * mean(my.data$meanmaxt[range,]^2) +
          sims.mat[,sprintf('psi.meanprec[%s]', ss)]    * mean(my.data$meanprec[range,])  +
          sims.mat[,sprintf('psi.floral[%s]', ss)]      * xx +
          sims.mat[,        'psi.sitearea']             * mean(my.data$sitearea[range])
      } 
      sapply(x.vals, get.y.val)
    }
    message("Extracting chains for floral resources, site resolution: ", resolutions[ii])
    ## extract chains by species
    chains.by.sp <- lapply(1:my.data$nspecies, get.for.sp)
    ## and then stitch together into a multidimensional array
    chains.arr <- abind(chains.by.sp, along=3)
    ## now we can use this to make our plot
    if(ii==1){
      plot(NA,
           xlim = range(x.vals),
           ylim = c(0,1),
           xlab = 'Floral Resources',
           ylab = NA,
           las  = 1,
           xaxt = "n",
           cex.lab  = cex.lab,
           cex.axis = cex.axis)
    }
    add.line <- function(xx, bci=FALSE, cc, lwd=1, ...) {
      to.plot <- apply(xx, 2, function(ii)
        expit(c(mean = mean(ii), quantile(ii, probs = c(0.025,0.975)))))
      ## add mean line
      lines(x = x.vals, y = to.plot["mean",], type = "l", col = cc,lwd = lwd, lty = linetypes[ii])
      ## add BCI, if desired
      if(bci) {
        polygon(x   = c(x.vals, rev(x.vals)),
                y   = c(to.plot["2.5%",],rev(to.plot["97.5%",])),
                col = makeTransparent(sprintf("%s",cc), alpha=0.25),
                border = NA)
      }
    }
    ## add the mean line across species
    xx <- apply(chains.arr, 1:2, mean)
    add.line(xx, bci = TRUE, cc = cols[ii], lwd = 2, lty = 2)
    
    if(ii==length(resolutions)){
      ## Final plot arguments 
      op <- par(new = TRUE) #does this command, but returns the previous par
      plot(NA,
           xlim = range(JAGS.era.unscaled[,,"floral.all"]), 
           ylim = c(0,1),
           xlab = '', 
           ylab = '',
           las  = 1)
      par(op) #calling the previous par, not the new=T one 
      
      add.fig.label(fig.num = 4,extra = 0.065)
    }
  }
  message("Calling dev.off()")
  dev.off() 
}

## So now we've made that whole function that contains all the
## plotting code, now here is where we're going to actually call it!
## Make sure at the top you've specified your location for the file! 

make.resolutions.figure(full     = TRUE, 
                        cex.lab  = 1.75, 
                        cex.axis = 1.25,
                        cex.main = 1)
