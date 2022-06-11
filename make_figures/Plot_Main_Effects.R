## Written by Hanna M. Jackson - hmj2@sfu.ca

## Put your working directory here: 
workingdirectory <- 'PUT YOUR WORKING DIRECTORY HERE'
setwd(workingdirectory)
message("working directory: ",getwd())


make.main.figure <- function(full=FALSE, site.resolution, cex.lab,
                             cex.axis, cex.main){
  pdf(file   = sprintf("MainEffects_%s.pdf", site.resolution),
      width  = 9,
      height = 9)
  
  ## Setup the layout of the multipannel plots
  m <- matrix(c(1,2,3,4), ncol = 2, nrow = 2, byrow = TRUE)
  layout(m)

  ## Plotting parameters
  par(oma = c(0.2, 0.2, 0, 0.2), 
      mar = c(4.5, 4, 0.1, 0.1),
      mgp = c(2,0.2,0),
      tcl = 0,
      cex.axis = cex.axis,
      cex.main = cex.main,
      pty = 's')

  ## ~~~~~~~~ Initiation  - loading and sourcing files ~~~~~~~~~

  ## Some required functions 
  expit <- function(x) 1/(1+exp(-x))
  logit <- function(x) log(x/(1-x))  


  ## Set the names and locations of the model runs here: 
  if(site.resolution == 100){
    files <- list(covariate = "output/modeloutput-env-100.RData", 
                  era =       'output/modeloutput-era-100.RData') 
  }

  if(site.resolution == 50){
    files <- list(covariate = 'output/modeloutput-env-50.RData',  
                  era =       'output/modeloutput-era-50.RData')  
  }

  if(site.resolution == 250){
    files <- list(covariate = 'output/modeloutput-env-250.RData', 
                  era =       'output/modeloutput-era-250.RData') 
  }
  load(sprintf("data_prep/saved/for_analysis/JAGS.era-%s.RData",site.resolution),verbose=T)
  JAGS.era.unscaled <- JAGS.era
  
  library("abind")
  fig.labs <- c("a","b","c","d")
  add.fig.label <- function(fig.num,extra){
    text(x = par('usr')[2]-par('usr')[2]*(0.08)+extra,
         y = par('usr')[4]*0.95, 
         fig.labs[fig.num],
         cex = 1.5, 
         pos = 4
         )
  }

  ## Function that makes a colour transparent
  makeTransparent <- function(..., alpha=0.5) {
    if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
    alpha = floor(255*alpha)  
    newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
    .makeTransparent = function(col, alpha) {
      rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
    }
    apply(newColor, 2, .makeTransparent, alpha=alpha)
  }
  
  
  
  ##~~~~~~~~~~~ PANNEL 1: Occupancy vs Era  ~~~~~~~~~~~~~

  message("Loading file: ", files[["era"]])
  load(files[["era"]], verbose=TRUE)
  summ     <- bugs$BUGSoutput$summary
  sims.arr <- bugs$BUGSoutput$sims.array
  sims.mat <- bugs$BUGSoutput$sims.matrix  
  
  ## figure out x-values for our figure
  x.vals <- seq(from = 1, 
                to   = 6,
                length = 6)

  ## function to extract the chains of y-values for a given species
  get.for.sp <- function(ss) {
    get.y.val <- function(xx) {
      sp.range <- my.data$species.ranges[ss,]
      chains <- sims.mat[,'psi.0']                   +
        sims.mat[,sprintf('psi.sp[%d]', ss)]         +       
        sims.mat[,sprintf('psi.era[%d]', ss)]   * xx +
        sims.mat[,        'psi.sitearea']       * mean(my.data$sitearea)         
    }
    sapply(x.vals, get.y.val)
  }
  message("Extracting chains for era")
  ## extract chains by species
  chains.by.sp <- lapply(1:my.data$nspecies, get.for.sp)
  ## and then stitch together into a multidimensional array
  chains.arr <- abind(chains.by.sp, along=3)
  
  ## now we can use this to make our plot
  plot(NA,
       xlim = range(x.vals),
       ylim = c(0,1),
       xlab = '',
       ylab = 'Occupancy Probability',
       las  = 1,
       xaxt = "n",
       cex.lab  = cex.lab,
       cex.axis = cex.axis)

  ## add a species-specific line
  add.line <- function(xx, bci=FALSE, cc, lwd=1, ...) {
    to.plot <- apply(xx, 2, function(ii)
      expit(c(mean = mean(ii), quantile(ii, probs = c(0.025,0.975)))))
    ## add mean line
    lines(x = x.vals, y = to.plot["mean",], type = "l", col = cc, lwd = lwd)
    ## add BCI, if desired
    if(bci) {
      polygon(x   = c(x.vals, rev(x.vals)),
              y   = c(to.plot["2.5%",],rev(to.plot["97.5%",])),
              col = makeTransparent(sprintf("%s",cc), alpha=0.25),
              border = NA)
    }
  }
  apply(chains.arr, 3, add.line, cc='gray', lwd=2)
  ## add the mean line across species
  xx <- apply(chains.arr, 1:2, mean)
  add.line(xx, bci=TRUE, cc='red', lwd=2, lty=2)
  
  ## Final plot arguments 
  op <- par(new=TRUE) #does this command, but returns the previous par
  labs <- c("1900-1920",
            "1920-1940",
            "1940-1960",
            "1960-1980",
            "1980-2000",
            "2000-2020")
  text(x = seq_along(labs)+0.1, 
       y = par('usr')[3]-0.01,
       labels = labs,
       srt = 20,
       adj = c(1.1,1.1),
       xpd = NA,
       cex = 0.9)
  axis(tick=T,side=1,labels=F)

  mtext("Era",
        side = 1,
        line = 3,
        cex = 1.25) 
  
  par(op) #calling the previous par, not the new=T one 
  add.fig.label(fig.num=1,extra=0)

  
  
  ## ~~~~~~~~~~~~~~~~ PANNELS 2-4 Occupancy vs Covariates ~~~~~~~~~~~~~~~~~~
  message("Loading file: ", files[["covariate"]])
  load(files[["covariate"]], verbose=TRUE)
  summ     <- bugs$BUGSoutput$summary
  sims.arr <- bugs$BUGSoutput$sims.array
  sims.mat <- bugs$BUGSoutput$sims.matrix 
  
  
  ## ~~~~~~~~~~~~~~~~ PANNEL 2: Occupancy vs Temperature ~~~~~~~~~~~~~~~~~~

  ## figure out x-values for our figure
  x.vals <- seq(from = min(my.data$meanmaxt),
                to   = max(my.data$meanmaxt),
                length.out = 101)

  ## function to extract the chains of y-values for a given species
  get.for.sp <- function(ss) {
    get.y.val <- function(xx) {
      sp.range <- my.data$species.ranges[ss,]
      chains <- sims.mat[,'psi.0'] +
        sims.mat[,sprintf('psi.sp[%s]', ss)] +
        sims.mat[,sprintf('psi.meanmaxt[%s]', ss)] * xx +
        sims.mat[,        'psi.meanmaxt.sq']       * xx * xx +
        sims.mat[,sprintf('psi.meanprec[%s]', ss)] * mean(my.data$meanprec[sp.range,]) +
        sims.mat[,sprintf('psi.floral[%s]', ss)]   * mean(my.data$floral[sp.range,]) +
        sims.mat[,        'psi.sitearea']          * mean(my.data$sitearea[sp.range])
    } 
    sapply(x.vals, get.y.val)
  }
  message("Extracting chains for temperature")
  ## extract chains by species
  chains.by.sp <- lapply(1:my.data$nspecies, get.for.sp)
  ## and then stitch together into a multidimensional array
  chains.arr <- abind(chains.by.sp, along=3)
  ## now we can use this to make our plot

  plot(NA,
       xlim = range(x.vals),
       ylim = c(0,1),
       xlab = 'Temperature (°C)',
       ylab = NA,
       las  = 1,
       xaxt = "n",
       cex.lab  = cex.lab,
       cex.axis = cex.axis)

  ## add a species-specific line
  add.line <- function(xx, sp.range=NULL, avg.trend=FALSE, cc, lwd=1, ...) {
    to.plot <- apply(xx, 2, function(ii)
      expit(c(mean = mean(ii), quantile(ii, probs = c(0.025,0.975)))))

    ## add mean line
    if(!avg.trend) {
      sp.range.vals <- range(my.data$meanmaxt[sp.range,])
      out.range.low  <- x.vals <= min(sp.range.vals)
      out.range.high <- x.vals >= max(sp.range.vals)
      in.range <- x.vals > min(sp.range.vals) & x.vals < max(sp.range.vals)
      lines(x = x.vals[in.range],
            y = to.plot["mean",in.range],
            type = "l",
            col  = 'gray',
            lwd  = lwd)
    }  
    ## add community trend
    if(avg.trend) {
      lines(x=x.vals, y=to.plot["mean",], type="l", col=cc, lwd=lwd)
      polygon(x   = c(x.vals, rev(x.vals)),
              y   = c(to.plot["2.5%",],rev(to.plot["97.5%",])),
              col = makeTransparent(sprintf("%s",cc), alpha=0.25),
              border = NA)
    }
  }
  add.line.for.sp <- function(sp) {
    sp.range <- my.data$species.ranges[sp,]
    add.line(chains.arr[,,sp], sp.range = sp.range, cc = 'gray', lwd = 2)
  }

  sapply(1:dim(my.data$species.ranges)[1], add.line.for.sp)
  
  ## add the mean line across species
  xx <- apply(chains.arr, 1:2, mean)
  add.line(xx, avg.trend = TRUE, cc = 'black', lwd = 2.5, lty = 2)
  
  ## Final plot arguments 
  op <- par(new = TRUE) #does this command, but returns the previous par
  plot(NA,
       xlim = range(JAGS.era.unscaled[,,'meanmaxt']),
       ylim = c(0,1),
       xlab = '',
       ylab = '',
       las  = 1)
  par(op) #calling the previous par, not the new=T one 
  add.fig.label(fig.num = 2, extra = 0)



  
  ## ~~~~~~~~~~~~ PANNEL 3: Occupancy vs Precipitation ~~~~~~~~~~~~~        

  ## figure out x-values for our figure
  x.vals <- seq(from = min(my.data$meanprec),
                to   = max(my.data$meanprec),
                length.out = 1001)

  ## function to extract the chains of y-values for a given species
  get.for.sp <- function(ss) {
    get.y.val <- function(xx) {
      sp.range <- my.data$species.ranges[ss,]
      chains <- sims.mat[,'psi.0'] +
        sims.mat[,sprintf('psi.sp[%s]', ss)] +
        sims.mat[,sprintf('psi.meanmaxt[%s]', ss)]    * mean(my.data$meanmaxt[sp.range,]) +
        sims.mat[,        'psi.meanmaxt.sq']          * mean(my.data$meanmaxt[sp.range,]^2) +
        sims.mat[,sprintf('psi.meanprec[%s]', ss)]    * xx +
        sims.mat[,sprintf('psi.floral[%s]', ss)]      * mean(my.data$floral[sp.range,]) +
        sims.mat[,        'psi.sitearea']             * mean(my.data$sitearea[sp.range])
    }
    sapply(x.vals, get.y.val)
  }
  message("Extracting chains for precipitation")
  ## extract chains by species
  chains.by.sp <- lapply(1:my.data$nspecies, get.for.sp)
  ## and then stitch together into a multidimensional array
  chains.arr <- abind(chains.by.sp, along=3)
  ## now we can use this to make our plot

  plot(NA,
       xlim = range(x.vals),
       ylim = c(0,1),
       xlab = 'Precipitation (kg/m²)',
       ylab = 'Occupancy Probability',
       las  = 1,
       xaxt = "n",
       cex.lab = cex.lab,
       cex.axis = cex.axis)

  ## add a species-specific line
  add.line <- function(xx, sp.range = NULL, avg.trend = FALSE, cc, lwd = 1, ...) {
    to.plot <- apply(xx, 2, function(ii)
      expit(c(mean = mean(ii), quantile(ii, probs = c(0.025,0.975)))))

    ## add mean line
    if(!avg.trend) {
      sp.range.vals <- range(my.data$meanprec[sp.range,])
      out.range.low  <- x.vals <= min(sp.range.vals)
      out.range.high <- x.vals >= max(sp.range.vals)
      in.range <- x.vals > min(sp.range.vals) & x.vals < max(sp.range.vals)
      lines(x = x.vals[in.range],
            y = to.plot["mean",in.range],
            type = "l",
            col  = 'gray',
            lwd  = lwd)
    }
    
    ## add community trend
    if(avg.trend) {
      lines(x = x.vals,
            y = to.plot["mean",],
            type = "l",
            col  = cc,
            lwd  = lwd)
      polygon(x   = c(x.vals, rev(x.vals)),
              y   = c(to.plot["2.5%",],rev(to.plot["97.5%",])),
              col = makeTransparent(sprintf("%s",cc), alpha=0.25),
              border = NA)
    }
    
  }
  add.line.for.sp <- function(sp) {
    sp.range <- my.data$species.ranges[sp,]
    add.line(chains.arr[,,sp], sp.range = sp.range, cc = 'gray', lwd = 2)
  }

  sapply(1:dim(my.data$species.ranges)[1], add.line.for.sp)

  ## add the mean line across species
  xx <- apply(chains.arr, 1:2, mean)
  add.line(xx, avg.trend = TRUE, cc = 'black', lwd = 2.5, lty = 2)
  
  ## Final plot arguments 
  op <-par(new = TRUE) #does this command, but returns the previous par
  plot(NA, xlim = range(JAGS.era.unscaled[,,'meanprec']), 
       ylim = c(0,1),
       xlab = '', 
       ylab = '',
       las  = 1)
  par(op) # calling the previous par, not the new=T one 
  
  add.fig.label(fig.num = 3,extra = 0)
  


  
  ## ~~~~~~~~~~~~~ PANNEL 4: Occupancy vs Floral Resources ~~~~~~~~~~~~~~~~~~       

  ## figure out x-values for our figure
  x.vals <- seq(from = min(my.data$floral),
                to   = max(my.data$floral),
                length.out = 1001)

  ## function to extract the chains of y-values for a given species
  get.for.sp <- function(ss) {
    get.y.val <- function(xx) {
      sp.range <- my.data$species.ranges[ss,]
      chains <- sims.mat[,'psi.0'] +
        sims.mat[,sprintf('psi.sp[%s]', ss)] +
        sims.mat[,sprintf('psi.meanmaxt[%s]', ss)]    * mean(my.data$meanmaxt[sp.range,]) +
        sims.mat[,        'psi.meanmaxt.sq']          * mean(my.data$meanmaxt[sp.range,]^2) +
        sims.mat[,sprintf('psi.meanprec[%s]', ss)]    * mean(my.data$meanprec[sp.range,])  +
        sims.mat[,sprintf('psi.floral[%s]', ss)]      * xx +
        sims.mat[,        'psi.sitearea']             * mean(my.data$sitearea[sp.range])
    } 
    sapply(x.vals, get.y.val)
  }
  message("Extracting chains for floral resources")
  ## extract chains by species
  chains.by.sp <- lapply(1:my.data$nspecies, get.for.sp)
  ## and then stitch together into a multidimensional array
  chains.arr <- abind(chains.by.sp, along=3)
  ## now we can use this to make our plot

  plot(NA,
       xlim = range(x.vals),
       ylim = c(0,1),
       xlab = 'Floral Resources',
       ylab = NA,
       las  = 1,
       xaxt = "n",
       cex.lab  = cex.lab,
       cex.axis = cex.axis)

  ## add a species-specific line
  add.line <- function(xx, sp.range = NULL, avg.trend = FALSE, cc, lwd = 1, ...) {
    to.plot <- apply(xx, 2, function(ii)
      expit(c(mean = mean(ii), quantile(ii, probs = c(0.025,0.975)))))

    ## add mean line
    if(!avg.trend) {
      sp.range.vals <- range(my.data$floral[sp.range,])
      out.range.low  <- x.vals <= min(sp.range.vals)
      out.range.high <- x.vals >= max(sp.range.vals)
      in.range <- x.vals > min(sp.range.vals) & x.vals < max(sp.range.vals)
      lines(x = x.vals[in.range],
            y = to.plot["mean",in.range],
            type = "l",
            col  = 'gray',
            lwd  = lwd)
    }
    
    ## add community trend
    if(avg.trend) {
      lines(x = x.vals,
            y = to.plot["mean",],
            type = "l",
            col  = cc,
            lwd  = lwd)
      polygon(x = c(x.vals, rev(x.vals)),
              y = c(to.plot["2.5%",],rev(to.plot["97.5%",])),
              col = makeTransparent(sprintf("%s",cc), alpha=0.25),
              border = NA)
    }
  }

  add.line.for.sp <- function(sp) {
    sp.range <- my.data$species.ranges[sp,]
    add.line(chains.arr[,,sp], sp.range = sp.range, cc = 'gray', lwd = 2)
  }

  sapply(1:dim(my.data$species.ranges)[1], add.line.for.sp)
  
  ## add the mean line across species
  xx <- apply(chains.arr, 1:2, mean)
  add.line(xx, avg.trend = TRUE, cc = 'black', lwd = 2.5, lty = 2)

  ## Final plot arguments 
  op <- par(new = TRUE) # does this command, but returns the previous par
  plot(NA,
       xlim = range(JAGS.era.unscaled[,,"floral.all"]), 
       ylim = c(0,1),
       xlab = '', 
       ylab = '',
       las = 1)
  par(op) # calling the previous par, not the new=T one 
  add.fig.label(fig.num = 4, extra = 0.065) 
  dev.off() 
}


## Now call that HUGE funtion we just made for each site resolution we
## want to run it for
for (i in c(50,250)){
  message("Starting resolution ", i)
  make.main.figure(full = TRUE, 
                   site.resolution = i,
                   cex.lab  = 1.5,
                   cex.axis = 1.25, 
                   cex.main = 1)
}
