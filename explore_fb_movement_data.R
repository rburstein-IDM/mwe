######################################################################
## Author:  Roy Burstein
## Date:    July 2019
## Purpose: Explore Facebook movement data for Seattle to get a better
##          handle on whats in there
######################################################################


rm(list=ls())

# load libs
library(data.table)
library(raster)
library(wellknown)
library(leaflet)
library(geosphere)
library(leaflet.minicharts)
library(ggplot2)
library(plotly)
library(sf)
library(grid)
library(gridExtra)

# get wd (in the AMUG Dropbox)
user <- Sys.info()[['user']]
setwd(sprintf('C:/Users/%s/Dropbox (IDM)/AMUG/roy/scratch/FB-Seattle-Data-June-28-2019', user))


# load in the movement data - they are csv's by 8 hour period
# unfortunately, it looks like we captured a weekend here (june 27th (thrs) through june 30th (sun))
d <- data.table()
for(f in list.files('./seattle_flu_map_june_28_2019_movement_between_tiles/', pattern='.csv')){
  tmp      <- fread(sprintf('./seattle_flu_map_june_28_2019_movement_between_tiles/%s', f))
  tmp[, date := strsplit(gsub('.csv','',f), '_')[[1]][1]]
  tmp[, time := strsplit(gsub('.csv','',f), '_')[[1]][2]]
  d        <- rbind(d, tmp)
}

# to keep things simple, lets just work with June 27th for now since its a thursday
d <- d[date == 20190627, ]

# The geographic line data are stored in Well-Known Text format as linestrings
#  to make them easier to work with for now, lets convert this into four points
#  representing the latitude and longitude of starting point (A) and ending point (B)
#  lonA, latA, lonB, latB. Use a helper function to extract these values out of WKT format
extr.coords <- function(x){
  wkt <- d$Geometry[x]
  out <- wkt2geojson(wkt)$geometry$coordinates
  return( data.table(lonA = out[1,1],
                     latA = out[1,2],
                     lonB = out[2,1],
                     latB = out[2,2]) )
}
coords <- lapply(1:nrow(d), FUN=extr.coords)
coords <- do.call('rbind', coords)

d[, Geometry := NULL]
d <- cbind(d, coords)


# Question: Do stayers (i.e. no tile movement) appear in the dataset?
# Answer:   No, there are no such rows. I am asking Alex at FB about this now.
d[ round(latA,3) == round(latB,3) & round(lonA,2) == round(lonB,3)]


# These data are generated as part of a crisis map, so there is a baseline mean movement and a 'crisis'
#  movement column. THe baseline represents the mean movement between tiles for the the 45 days before
#  the crisis map was kicked off (for that day and time). We will just use baseline movement for now
# Note that there is a masking at 10 people. If there is a movement vector < 10 people, 
#  they are excluded from this dataset
summary(d[['Baseline: People Moving']])
summary(d[['Crisis: People Moving']])

hist(d[['Baseline: People Moving']], breaks = 500)

# rename some key variables so they are easier to work with
setnames(d, c('Baseline: People Moving','Length(km)'),
            c('movers','length'))


# give IDs to the A and B tiles (the ones provided in this dataset dont appear to be unique)
coords <- unique(data.table(lat = c(d$latA,d$latB), lon = c(d$lonA,d$lonB)))
coords[, id := 1:.N]
d <- merge(d, coords, by.x = c('lonB','latB'), by.y = c('lon','lat'), all.x = TRUE)
setnames(d,'id','idB')
d <- merge(d, coords, by.x = c('lonA','latA'), by.y = c('lon','lat'), all.x = TRUE)
d[, idA := id]

### Lets try to visualize the data a bit
# mostly the same tos and fros around the puget sound region
plot(d$lonA,d$latA)
points(d$lonB,d$latB,col='red')

# add a basemap of seattle
dd=d[time == '0000']
m <- leaflet(data = dd) %>% # just morning commuters
  setView(lng = -122.5, lat = 47.6, zoom = 9) %>%
  addTiles() %>%
  addMarkers(~lonA, ~latA, popup=dd$idA) 
m


#  intensity of inflow or outflow

dd <- d[time=='1600',] #subset  (I think 0800 is the middle of the night)
#  try to visualize flows
flows <- gcIntermediate(dd[,c('lonA','latA')], dd[,c('lonB','latB')], sp = TRUE, addStartEnd = TRUE)
flows$time    <- dd$time
flows$start   <- dd[['Starting Region Name']]
flows$movers  <- dd$movers
flows$idA     <- dd$idA
flows$idB     <- dd$idB
hover     <- paste0(dd[['Starting Region Name']], ' (',
                    dd$idA,') to ', dd[['Ending Region Name']],
                       ' (',dd$idB,') ', ' | Movers: ',dd$movers)
  
m <- leaflet() %>% # just morning commuters
  setView(lng = -122.25, lat = 47.6, zoom = 11) %>%
  addTiles() %>%
  addPolylines(data = flows, 
               weight = ~movers/20, label = hover, 
               group = ~idA, color = start,opacity=.2) #%>%
m

# another way of looking at it
m <- leaflet() %>% # just morning commuters
  setView(lng = -122.25, lat = 47.6, zoom = 11) %>%
  addTiles() %>%
  addFlows(
    lat0 = dd$latA + rnorm(nrow(dd), sd=0.001),  
    lng1 = dd$lonB + rnorm(nrow(dd), sd=0.001),  
    lng0 = dd$lonA + rnorm(nrow(dd), sd=0.001), 
    lat1 = dd$latB + rnorm(nrow(dd), sd=0.001), 
    flow = dd$movers, 
    opacity      = 0.3,
    maxThickness = 3)
m



# in reality, these are centroids of pixels. Can try to visualize them as such, showing


# make a heatmap matrix of comings and goings. 
shp <- st_read('./seattleshp/City_Clerk_Neighborhoods.shp')
dd <- highlight_key(d, ~id)
g1 <-
ggplot(data = dd, aes(x=idA, y=idB, size=movers, color=time, group=id)) + 
  geom_point(alpha=.5) + xlab('Coming') + ylab('Going') + theme_bw()
g2 <-
  ggplot(dd) + #geom_sf(data=shp) + 
  geom_point(aes(lonA,latA,group=id)) + theme_bw() 


s <- subplot(g1,g2,nrows=2,shareY=FALSE,shareX=FALSE)
highlight(s, on = "plotly_selected", selectize = TRUE)








# get basemap
# register_google(key='AIzaSyCcGMG6iDYb7K1ITyW2J6KgWHitbZ-Y2tU')
# base <- ggmap(
#   get_stamenmap(c(left   = min(coords$lon),
#                 right  = max(coords$lon),
#                 bottom = min(coords$lat),
#                 top    = max(coords$lat)), 
#               zoom = 9, maptype = "toner-lite"))
# g2= base + geom_point(data=coords, aes(lon,lat,color=id))


