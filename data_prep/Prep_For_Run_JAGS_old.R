# For running on lab mac: 
#setwd('/Users/hannajackson/Dropbox/Bumble_Bee_Climate/climate-bumblebee')

# remember to SAVE ALL and source files before running anything - or else your changes won't be run!

source('init.R')
load(sprintf("JAGS.arr-%s.RData",site.resolution),   verbose=TRUE) # occurence data
load(sprintf('JAGS.visit-%s.RData', site.resolution), verbose=TRUE)
load(sprintf("JAGS.site-%s.RData", site.resolution),  verbose=TRUE)
load(sprintf("JAGS.era-%s.RData", site.resolution), verbose=TRUE)

#load("species.ranges.RData", verbose=TRUE)
load(sprintf("species.ranges-%s.RData", site.resolution), verbose=TRUE)
load(sprintf("master.index-%s.RData", site.resolution), verbose=TRUE)
load(sprintf("visit.arr-%s.RData", site.resolution), verbose=TRUE) 

#JAGS
my.inits<-function(){
   tmp <- array(as.matrix(species.ranges*1), dim = c(nspecies, nsite, nera)) # takes species x site 1s and 0s (1= yes in species range according to 9.Make_Range.R) and then just duplicates it nera times
   tmp2 <- aperm(tmp, c(2,3,1))  # makes it so the dimensions are the same as in JAGS.arr 
   list(Z=tmp2)    
}
my.inits2<-function(){
   tmp <- array(as.matrix(species.ranges*1), dim = c(nspecies, nsite, nera)) # takes species x site 1s and 0s (1= yes in species range according to 9.Make_Range.R) and then just duplicates it nera times
   tmp2 <- aperm(tmp, c(2,3,1))  # makes it so the dimensions are the same as in JAGS.arr 
   list(Z=tmp2,sd=1)    
}

my.inits3<-function(){
   tmp <- array(as.matrix(species.ranges*1), dim = c(nspecies, nsite, nera)) # takes species x site 1s and 0s (1= yes in species range according to 9.Make_Range.R) and then just duplicates it nera times
   tmp2 <- aperm(tmp, c(2,3,1))  # makes it so the dimensions are the same as in JAGS.arr 
   list(Z=tmp2,sd=1,mu.psi.meanmaxt.sq=-2,sigma.psi.meanmaxt.sq=0.01) 
}

if(scale.method=="poly"){
   
   ## ------------------------------------------------------------
   ## add scaled variables (linear AND squared) via poly function

   library(abind)
   JAGS.era <- abind(JAGS.era,
                     meanmaxt.poly.li=array(NA,dim=dim(JAGS.era)[1:2]),
                     meanmaxt.poly.sq=array(NA,dim=dim(JAGS.era)[1:2]),
                     meanprec.poly.li=array(NA,dim=dim(JAGS.era)[1:2]),
                     meanprec.poly.sq=array(NA,dim=dim(JAGS.era)[1:2]),
                     floral.poly.li=array(NA,dim=dim(JAGS.era)[1:2]),
                     floral.poly.sq=array(NA,dim=dim(JAGS.era)[1:2]),
                     along=3)
   
   meanmaxt.poly <- poly(as.vector(JAGS.era[,,'meanmaxt']),2)
   JAGS.era[,,'meanmaxt.poly.li'] <- meanmaxt.poly[,1] #first column of this poly(temp)
   JAGS.era[,,'meanmaxt.poly.sq'] <- meanmaxt.poly[,2] #second column of this poly(temp)
   
   meanprec.poly<-poly(as.vector(JAGS.era[,,'meanprec']),2)
   JAGS.era[,,'meanprec.poly.li'] <- meanprec.poly[,1]
   JAGS.era[,,'meanprec.poly.sq'] <- meanprec.poly[,2]
   
   floral.poly<-poly(as.vector(JAGS.era[,,'floral.all']),2)
   JAGS.era[,,'floral.poly.li'] <- floral.poly[,1]
   JAGS.era[,,'floral.poly.sq'] <- floral.poly[,2]
   
   ## save this too
   out.poly <- list(temp=meanmaxt.poly,
                    prec=meanprec.poly,
                    flor=floral.poly)
   save(out.poly, file=sprintf("out.poly-%s.RData",site.resolution))
   
   my.data <- list(X                 = JAGS.arr[master.index], 
                   nsite             = dim(JAGS.arr)[1],
                   nvisit            = dim(JAGS.arr)[2],
                   nera              = dim(JAGS.arr)[3],
                   nspecies          = dim(JAGS.arr)[4],
                   nind              = nrow(master.index),
                   master.index      = master.index,
                   
                   site              = site,
                   sitevec           = site.vec,
                   visit             = visit,
                   visitvec          = visit.vec,
                   era               = era,
                   eravec            = era.vec,
                   species           = species,
                   
                   sitearea          = JAGS.site[,"sitearea"],
                   
                   meanmaxt.raw     = JAGS.era[,,"meanmaxt"],
                   meanmaxt.li      = JAGS.era[,,'meanmaxt.poly.li'],  
                   meanmaxt.sq      = JAGS.era[,,"meanmaxt.poly.sq"],  
                   
                   meanprec.raw     = JAGS.era[,,"meanprec"],
                   meanprec.li      = JAGS.era[,,"meanprec.poly.li"],  
                   meanprec.sq      = JAGS.era[,,"meanprec.poly.sq"],  
                   
                   floral.raw       = JAGS.era[,,"floral.all"],
                   floral.li        = JAGS.era[,,"floral.poly.li"],  
                   floral.sq        = JAGS.era[,,"floral.poly.sq"],  
    
                   species.ranges    = species.ranges,               
                   visit.arr         = visit.arr
   )
}




if(scale.method=="scale"){
   my.scale<-function(x){
      return((x-mean(x,na.rm=T))/sd(x,na.rm=T))
   }
   
   #Modification to the data - scaling variables and ensuring class = numeric 
   JAGS.site[,"sitearea"]    <- my.scale(JAGS.site[,"sitearea"])
   #JAGS.visit[,,,'nsources']  <-as.numeric(scale(JAGS.visit[,,,'nsources']))
   
   JAGS.era[,,"meanmaxt"]     <- my.scale(JAGS.era[,,"meanmaxt"])
   JAGS.era[,,"meanprec"]     <- my.scale(JAGS.era[,,"meanprec"])
   JAGS.era[,,"landuse"]      <- my.scale(JAGS.era[,,"landuse"])
   JAGS.era[,,"floral.all"]   <- my.scale(JAGS.era[,,"floral.all"])
   # Do not include a scale for the poly data 
   #scale(poly(c(1,2,3),2)[,1])
   #scale(scale(c(1,2,3)))
   #poly(c(1,2,3))
   
   my.data <- list(X                 = JAGS.arr[master.index], #HOW/WHY DOES THIS WORK???
                   nsite             = dim(JAGS.arr)[1],
                   nvisit            = dim(JAGS.arr)[2],
                   nera              = dim(JAGS.arr)[3],
                   nspecies          = dim(JAGS.arr)[4],
                   nind              = nrow(master.index),
                   master.index      = master.index,
                   
                   site              = site,
                   sitevec           = site.vec,
                   visit             = visit,
                   visitvec          = visit.vec,
                   era               = era,
                   eravec            = era.vec,
                   species           = species,
                   
                   sitearea          = JAGS.site[,"sitearea"],
                   
                   meanmaxt          = JAGS.era[,,"meanmaxt"],
                   #meanmint         = JAGS.era[,,'meanmint'],
                   meanprec          = JAGS.era[,,"meanprec"],
                   #nsources         = JAGS.visit[,,'nsources'],
                   landuse           = JAGS.era[,,'landuse'],
                   floral            = JAGS.era[,,"floral.all"],     
                   
                   species.ranges    = species.ranges,               
                   visit.arr         = visit.arr
                   )



   
}



save(my.inits, file="my.inits.RData")
save(JAGS.site, file="JAGS.site.scaled.RData")
save(JAGS.visit, file="JAGS.visit.scaled.RData")
save(JAGS.era, file=sprintf("JAGS.era.scaled-with.poly-%s.RData",
                            site.resolution))


