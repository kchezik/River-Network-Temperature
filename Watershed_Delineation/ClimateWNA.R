library(raster);library(sp);library(rgdal);library(GISTools);library(maptools) #GIS tools.
library(tidyverse)

setwd("~/sfu/River_Network_Temperature_Output/shp_climate_grid/")

#Coordinate reference systems.
WGS84 = crs("+init=epsg:4326"); BCAlbers = crs("+init=epsg:3005") 

#Clip grid by fraser basin.
thompson = readOGR(dsn = "./", layer = "pt_grid", p4s = c("+init=epsg:4326"))

#Add columns for Climate WNA
names(thompson)[1] = "ID1"
thompson$ID2 = thompson$ID1
thompson$lat = thompson@coords[,2]
thompson$long = thompson@coords[,1]

#Extract elevation data from the thompson grid dataset.
dem = raster::raster(x = "~/Documents/BC_DEM25/bcdem.tif")
thompson$el = raster::extract(x = dem, y = thompson, method = 'simple')

#Convert to dataframe and save as CSV.
thompson = data.frame(thompson) %>% select(-6,-7,-8)
write_csv(thompson, "~/sfuvault/River-Network-Temperature/Watershed_Delineation/climateBC.csv")


#Split up ClimateWNA results.
annual = read_csv("~/sfu/River_Network_Temperature_Output/predictors/climateBC_2014-2017YT.csv")
season = read_csv("~/sfu/River_Network_Temperature_Output/predictors/climateBC_2014-2017ST.csv")

year_save = function(x,y){
	for(i in unique(x$Year)){
		sub = x %>% filter(Year == i)
		write_csv(x = sub, path = paste("~/sfu/River_Network_Temperature_Output/predictors/",y,as.character(i),".csv", sep = ""))
	}
}

year_save(annual,"YT")
year_save(season,"ST")
