## ###########################################################
## Author:  Roy Burstein (rburstein@idmod.org)
## Date:    June 2019
## Purpose: Aggregate pixels to districts
## Input:   global dist shapefile, DPT3 coverage raster, 
##          worldpop raster
## Output:  Mean coverage estimates for each district by year
## Note:    
## ###########################################################



# load libraries
libs <- c('rgdal', 'raster', 'data.table', 'sf', 'sp','fasterize')
for(l in libs) library(l, character.only = TRUE)


# set paths 
user     <- Sys.info()['login']
root     <- sprintf('C:/Users/%s/Dropbox (IDM)/AMUG/roy/scratch/district_agg', user) # assuming you are linked to AMUG dropbox
datpath  <- sprintf('%s/data/',root) 
datpath <- 'C:/Users/rburstein/Dropbox (IDM)/AMUG/vxdel_equity/pathways/data/rasters/'
outpath  <- sprintf('%s/output/',root) 
shp_path <- "Q:/Data/Polis/polis_geodatabases/downloaded_20190211/WHO_POLIO_GLOBAL_GEODATABASE.gdb"


# load shapefile
shp <- st_read(dsn=shp_path, layer = 'GLOBAL_ADM2')
shp <- subset(shp, as.Date(ENDDATE) >= Sys.Date())

# unique dataset for shp and get an info shapefile
shp$id  <- 1:nrow(shp)
shpinfo <- as.data.table(shp)[, c('id',paste0('ADM',0:2,'_VIZ_NAME'))]


## load vaccine coverage rasters
dpt3 <- vdpt3 <- list()
for(f in list.files(path = sprintf('%s/vaccine/dpt3_coverage/mean',datpath))){
  n <- gsub(x=gsub(x=f,pattern='IHME_AFRICA_DPT_2000_2016_DPT3_COVERAGE_PREV_MEAN',replacement=''),
            pattern='_Y2019M04D01.TIF',replacement='')
  message(sprintf('Now loading: %s',n))
  dpt3[[n]]  <- raster(sprintf('%s/vaccine/dpt3_coverage/mean/%s',datpath,f))
  vdpt3[[n]] <- as.vector(dpt3[[n]])
}


## load population rasters
pop <- vpop <- list()
for(f in list.files(path = sprintf('%s/population/under5',datpath))){
  n <- gsub(x=gsub(x=f,pattern='worldpop_raked_a0004t_1y',replacement=''),
            pattern='_00_00.tif',replacement='')
  message(sprintf('Now loading: %s',n))
  pop[[n]]  <- raster(sprintf('%s/population/under5/%s',datpath,f))
  pop[[n]]  <- raster::crop(pop[[n]], extent(dpt3[[1]]))
  vpop[[n]] <- as.vector(pop[[n]])
}


# rasterize the gadm file
rshp <- fasterize(sf = shp, raster = dpt3[[1]], field = 'id')


# QA check dimensions all make sense
stopifnot(length(rshp) == length(dpt3[[1]]) & length(rshp) == length(pop[[1]]))


# assign admin ids, make long table, drop NAs
d <- data.table()
for(y in 2000:2016){ # 2016 is last year with all data.. end there for now
  message(sprintf('Binding year %s', y))
  tmp <- data.table(id   = as.vector(rshp), 
                    year = y, 
                    dpt3 = vdpt3[[paste0('_',y)]],
                    pop  = vpop[[paste0('_',y)]])
  tmp[,pxid:=1:.N]
  d <- rbind(d,tmp)
}

d <- na.omit(d)


# collapse to admin year and merge adm info
d_summary <- d[, .(dpt3_coverage = weighted.mean(dpt3,pop),
                   under5_pop    = sum(pop)), by = .(id,year)]
d_summary <- merge(d_summary, shpinfo, by = 'id', all.x = TRUE)


# save
write.csv(d_summary, sprintf('%s/dpt3_africa_district_summaries.csv', outpath))


