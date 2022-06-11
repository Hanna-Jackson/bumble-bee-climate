## Written by: Hanna M. Jackson - hmj2@sfu.ca 

make.sites <- function(site.resolution){

  ## Dependencies 
  library(maptools)
  library(parallel)
  library(raster)
  library(rgdal)
  library(rgeos)
  library(sp)
  library(spatstat)
  
  ## define initial and target map projections
  prj1 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  prj2 <- "+proj=aea +lat_1=20 +lat_2=60   +lat_0=40 +lon_0=-96  +x_0=0  +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"# Albers Conical Equal Area (Florida Geographic Data Library)
  
  ## load political borders for countries
  dd.countries <- readOGR(dsn='data_prep/data/Countries_WGS84')

  ## re-project 
  dd.countries <- spTransform(dd.countries, CRS(prj1)) 

  ## Crop the world to just Canada, US, and Mexico
  ## First make north america 
  dd.NA <- dd.countries[c(which(dd.countries@data$CNTRY_NAME=="Canada"), 
                          which(dd.countries@data$CNTRY_NAME=="United States"), 
                          which(dd.countries@data$CNTRY_NAME=="Mexico"))
                      , ]
  ## Make a box that will crop to precisely where we want 
  dd.box <- bbox2SP(n = 87,
                    e = -40,
                    s = 12,
                    w = -190,
                    proj4string=CRS(prj1))
  ## Do the actual cropping
  dd.AM <- raster::intersect(dd.NA, dd.box)
  
  ## Merging the countries - we dont want hard site boarders between them
  dd.AM@data$continent <- c("North America", "North America", "North America")
  dd.AM <- unionSpatialPolygons(dd.AM, dd.AM@data$continent)
  
  ## Project to the one we're going to use 
  dd.AM <- spTransform(dd.AM, CRS(prj2))

  message("Create raster cells")
  ## create raster cells
  ## units for this resolution is km squared 
  rr <- raster(extent(dd.AM),
               resolution = c(site.resolution*1000, site.resolution*1000),
               crs = prj2) #given us 10 cols and 10 rows
  rr[] <- 1:ncell(rr)
  dd.grid <- rasterToPolygons(rr, dissolve=TRUE)
  
  ## re-project grid so that it matches map
  dd.grid <- spTransform(dd.grid, crs(proj4string(dd.AM)))

  ## create sites by intersecting grid with map
  dd.sites <- raster::intersect(dd.grid, dd.AM)
  
  ## add site names
  dd.sites@data <- data.frame(site = paste('site.',
                                           formatC(seq(1:length(dd.sites)),# Format so that all sites have same number of digits 001,002 etc
                                                   width = nchar(length(dd.sites@polygons)),
                                                   format = "d",
                                                   flag = "0"),
                                           sep='')
                              ) 
  
  ## get site areas to account for differences in land areas 
  site.areas <- sapply(dd.sites@polygons, function(x) x@area)
  site.areas <- (as.matrix(site.areas))
  site.areas <- data.frame(area=site.areas, sitename=dd.sites@data)
  
  message("Saving Objects")
  save(dd.sites, site.areas, file = sprintf("data_prep/sites/%s/dd.sites.RData", site.resolution))
  save(dd.AM, file="data_prep/saved/dd.AM.RData")
}
