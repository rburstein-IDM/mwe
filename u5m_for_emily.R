# What percent of districts/provinces make up 


library(data.table)
library(sf)
library(ggplot2)
library(maps) 

setwd('C:/Users/rburstein/Dropbox (IDM)/IHME/data')

# load and trim data
d <- fread('U5MR_LMICS_admin2_data_full.csv')
d <- d[year == 2017]


# order and get cumulative deaths
d <- d[order(-u5_deaths_mean)]
d <- d[, cum_d  := cumsum(u5_deaths_mean)]
d <- d[, prop_d := cum_d/sum(u5_deaths_mean)]
d <- d[, cum_p  := cumsum(population_under5)]
d <- d[, prop_p := cum_p/sum(population_under5)]


# keep where under 50%
d[, first50 := prop_d<.500001]
d50 <- d[first50==TRUE]

# proportion of total districts
nrow(d50)/nrow(d)*100

# what countries are they in
table(d50$ADM0_NAME)

# load shapefile 
s <- read_sf('./shp/lbd_standard_admin_2.shp')
s <- merge(s, d, by = 'ADM2_CODE')

w <- sf::st_as_sf(map('world', plot = FALSE, fill = TRUE))

png('bottom50.png', height = 2000, width = 2400)
ggplot() + 
  geom_sf(data = s, aes(fill=first50), lwd = 0, colour = NA) +  
  geom_sf(data = w, fill = NA, lwd = 0.5) +
  coord_sf(xlim = c(-80, 140), ylim = c(-40, 50),
           expand = FALSE) + 
  theme_bw()
dev.off()

