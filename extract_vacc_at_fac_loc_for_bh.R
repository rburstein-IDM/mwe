

library(data.table)
library(raster)
library(ggplot2)

setwd('C:/Users/rburstein/OneDrive - IDMOD/code/mwe')

datpath  <- sprintf('C:/Users/rburstein/Dropbox (IDM)/AMUG/vxdel_equity/pathways/data/rasters/vaccine/')


# load fac data
d <- fread('epic_facility_locations.csv')


# load dpt3 data
m3 <- raster(sprintf('%s/dpt3_coverage/mean/IHME_AFRICA_DPT_2000_2016_DPT3_COVERAGE_PREV_MEAN_2013_Y2019M04D01.TIF',datpath))
u3 <- raster(sprintf('%s/dpt3_coverage/upper/IHME_AFRICA_DPT_2000_2016_DPT3_COVERAGE_PREV_UPPER_2013_Y2019M04D01.TIF',datpath))
l3 <- raster(sprintf('%s/dpt3_coverage/lower/IHME_AFRICA_DPT_2000_2016_DPT3_COVERAGE_PREV_LOWER_2013_Y2019M04D01.TIF',datpath))
m1 <- raster(sprintf('%s/dpt1_coverage/mean/IHME_AFRICA_DPT_2000_2016_DPT1_COVERAGE_PREV_MEAN_2013_Y2019M04D01.TIF',datpath))
u1 <- raster(sprintf('%s/dpt1_coverage/upper/IHME_AFRICA_DPT_2000_2016_DPT1_COVERAGE_PREV_UPPER_2013_Y2019M04D01.TIF',datpath))
l1 <- raster(sprintf('%s/dpt1_coverage/lower/IHME_AFRICA_DPT_2000_2016_DPT1_COVERAGE_PREV_LOWER_2013_Y2019M04D01.TIF',datpath))


# clean up lat and long
d[, latitude  := gsub('\\(','-',latitude)]
d[, latitude  := gsub('\\)','',latitude)]
d[, longitude := gsub('\\(','-',longitude)]
d[, longitude := gsub('\\)','',longitude)]
d[, latitude   := as.numeric(latitude)]
d[, longitude  := as.numeric(longitude)]


# extract values
d[, dpt3_mean  := extract(m3, cbind(longitude,latitude))]
d[, dpt3_lower := extract(l3, cbind(longitude,latitude))]
d[, dpt3_upper := extract(u3, cbind(longitude,latitude))]

d[, dpt1_mean  := extract(m1, cbind(longitude,latitude))]
d[, dpt1_lower := extract(l1, cbind(longitude,latitude))]
d[, dpt1_upper := extract(u1, cbind(longitude,latitude))]

# output
write.csv(d,'epic_facility_locations_with_dpt.csv')

