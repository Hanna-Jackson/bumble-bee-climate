## Written by Hanna M. Jackson - hmj2@sfu.ca

## Welcome! You have found Main.R
## This is the code that runs (almost) everything!

## 1.First we'll set your working directory ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## 2. Then we'll make an object called `args` ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## that we'll use to store all the various parameters that we'll use
## in both our data prep and the model run

## 3. Then we will prep the data! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## This section sources 10 files that go through each step of the data
## processing. When you source a file all it does is load in the
## function that that file has in it. Then all we have to do is call
## that function with the arguments we want and all the code runs. At
## the end of each file the important objects get saved to your
## computer and then when a different function needs that file it just
## loads it! 
##   IMPORTANT NOTE: 
## The data prep requires downloading all the climate data from the
## CHELSA website. If you don't want to do that, you can skip this
## section and move right on to model running and instead load the
## pre-processed data that I provide in the correct directory.

## 4. Then we will run the model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## This requires that you have jags downloaded onto your computer. You
## can find the download link at https://mcmc-jags.sourceforge.io 

## 5. Plotting!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## After you have your model output saved you can then go to the
## `make_figures` folder and make your plots! Make sure that before
## you run the figure code you tell the pdf() funtion the file path so
## it knows where to save

## 6. Thanks for taking the time to look through my work ~~~~~~~~~~~~~~~~~
## Any feedback on the code and analysis is appreciated!
## If you have issues, feel free to email me at hmj2@sfu.ca and I can
## try to help! Good luck!





## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~ 1. Set your working directory ~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

workingdirectory <- '~/Dropbox/Bumble_Bee_Climate/climate-bumblebee/for_github'
setwd(workingdirectory)
message("working directory set to: ",getwd())



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~ 2. Choose Parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Here we set the parameters of the run 
## See the code that these options are referenced in to know what they
## do specifically
## The parameters that are set right now are the ones we used in the
## main text

args <- list(
  ## Parameters used in multiple functions:
  site.resolution = 250, # The site resolution of our grid over North America, here 100 x 100 km
  n.cores.to.use  = 8,   # The number of cores your computer will use for multi-threaded tasks
  start.year      = 1901,# The starting year of the analysis - records before are filtered out
  
  ## format.data 
  visit.bin.time = 5,    # The number of years in a 'visit' to an 'era'. Here 5 years per each 20 year era
  era.bin.time   = 20,   # Number years in each era, here 20 
  
  ## clean.data
  is.subset       = FALSE, # subset the analysis down to a set number of records by setting this to TRUE
  subset.cap      = 0,     # if is.subset is TRUE, will subset to this number of records
  end.year        = 2020,  # THe ending year of the analysis - records after are filtered out
  split.occi      = TRUE,  # Split B. occidentalis into occidentalis and mckayi
  remove.problem.datasets = FALSE, # remove a few datasets that seemed to be suspicious
  
  ## make.bees.spatial
  site.era.min = 2,                                                    ## A site must have this many eras observed otherwise it is removed from the data
  remove.sites.without.certain.number.visits.in.any.era = TRUE,        ## Should sites be filtered out based on number of repeat visits in any era?
  min.site.visits.in.any.era = 2,                                      ## will use this value if the above is set to TRUE- Must be between 0-3 becaue there are only 4 visits/era
  remove.sites.without.certain.number.observations.in.all.eras = FALSE,## Should sites without a certain number of observations in all sites be removed?
  min.site.obs.in.all.eras = 5,                                        ## Will use this threshold if above is set to TRUE
  
  ## run.model
  model   = 'env_model', # The name of a file in the `models` folder that you want to run
  niter   = 4,     # Number of itterations of the MCMC to run
  nburnin = 1,     # Number of those niter that you want to be 'thrown away' at the beginning
  nthin   = 1,        # Thin the number of itterations saved by this many
  nchains = 3,         # Number of chains running at once
  notes   = " "        # Add notes for your own convenience here - doesnt DO anything other than this info will be stored
)



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~ 3. Data Prep ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## This first section of code preps all the data for the analysis
## See the code each of these in the `source` folder


## IMPORTANT NOTE:
## Code that involves the climate data will not be able to be run
## unless you download all the CHELSA climate data from their website.
## I have uploaded my version of what running this code would give
## you so you can run the 'run model' section without having to have
## your computer actually prep all the data 

source('data_prep/1.Make_Sites.R')         
make.sites  (site.resolution = args$site.resolution)

source('data_prep/2.Formatting.R')
format.data (visit.bin.time = args$visit.bin.time,
             era.bin.time   = args$era.bin.time,
             start.year     = args$start.year) 

source('data_prep/3.Cleaning.R') 
clean.data  (is.subset      = args$is.subset, 
             subset.cap     = args$subset.cap,
             end.year       = args$end.year,
             start.year     = args$start.year,
             split.occi     = args$split.occi,
             remove.problem.datasets = args$remove.problem.datasets) 

source('data_prep/4.Making_bees_spatial.R')
make.bees.spatial   (site.era.min    = args$site.era.min,
                     remove.sites.without.certain.number.visits.in.any.era=args$remove.sites.without.certain.number.visits.in.any.era,
                     min.site.visits.in.any.era = args$min.site.visits.in.any.era,
                     remove.sites.without.certain.number.observations.in.all.eras = args$remove.sites.without.certain.number.observations.in.all.eras,
                     min.site.obs.in.all.eras = args$min.site.obs.in.all.eras,
                     site.resolution = args$site.resolution,
                     n.cores.to.use  = args$n.cores.to.use) 

source('data_prep/5.Prep_Climate.R')
prep.climate(site.resolution = args$site.resolution)

source('data_prep/6.Filter_Bees.R')
filter.bees(site.resolution = args$site.resolution) 

source('data_prep/7.Making_Occ_Data.R')
make.occupancy.data(site.resolution = args$site.resolution) 

source('data_prep/8.Make_Range.R')
make.species.ranges(site.resolution = args$site.resolution) 

source('data_prep/9.Master.Index.R')
what.to.iterate.over(site.resolution = args$site.resolution)
                                       
source('data_prep/10.Making_Occ_Covariates.R')
make.covariates(site.resolution = args$site.resolution,
                n.cores.to.use  = args$n.cores.to.use)



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~ 4.Run Model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Source the file that loads and gets everything ready
site.resolution <- args$site.resolution
source('data_prep/11.Prep_For_Run_JAGS.R')

## Source the model file 
source(sprintf("models/%s.R", args$model))

## Visually check that these make sense and are indeed the options you want to run
dim(JAGS.arr)
message(args$site.resolution)
message(args$model)

## Run the model (Finally!)
run.start.time <- Sys.time()
run.date       <- Sys.Date()
message('\nRunning model\n ', args$model, ' with ', args$niter, ' iterations at size', args$site.resolution, '\n')
bugs <- jags(data       = my.data,
             inits      = my.inits,
             parameters.to.save = my.params, 
             model.file = my.model,
             n.iter     = args$niter, 
             n.burnin   = args$nburnin,
             n.thin     = args$nthin,
             n.chains   = args$nchains,
             working.directory = NULL)   

run.time <- as_hms(Sys.time() - run.start.time)

## Save everything
save(modelruns, file="modelruns.RData")
save(my.data, bugs, args, run.time, my.params, file='output/modeloutput.RData')
message('\nSaved: ', args$model, " at ",args$site.resolution, " site resolution with ", bugs$n.iter, " iterations")



