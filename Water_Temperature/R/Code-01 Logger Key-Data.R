#HOBO-Temperature Logger Data

#Key (Manual Transfer from Field Notes)
library(gtools)
load("HoboDat.RData")

site  = 86
data[["key"]] = dplyr::bind_rows(data$key,data.frame(site,logger.id = NA,deployment = NA,date.deployed = NA,date.retrieved = NA, latitude = NA, longitude = NA,retrieval.status = NA,moved = NA))

data$key$deployment[dim(data$key)[1]] = 1
data$key$logger.id[dim(data$key)[1]] = 10504065
data$key$date.deployed[dim(data$key)[1]] = "2014-08-08"
data$key$date.retrieved[dim(data$key)[1]] = NA
data$key$retrieval.status[dim(data$key)[1]] = 1
# 1 = present, 2 = present: out-of-water, 3 = present: buried, 4 = missing: river
# 5 = missing: human, 6 = missing: unknown
#data$key$latitude[dim(data$key)[1]] = 50.54352
#data$key$longitude[dim(data$key)[1]] = -119.0122
tail(data$key)

data[["key"]] = dplyr::bind_rows(data$key,data.frame(site,logger.id = NA,deployment = NA,date.deployed = NA,date.retrieved = NA, latitude = NA, longitude = NA,retrieval.status = NA,moved = NA))

data$key$deployment[dim(data$key)[1]] = 2
data$key$logger.id[dim(data$key)[1]] = 10602134
data$key$date.deployed[dim(data$key)[1]] = "2015-09-06"
data$key$date.retrieved[dim(data$key)[1]] = "2015-10-24"
data$key$retrieval.status[dim(data$key)[1]] = 2
data$key$moved[dim(data$key)[1]] = 1
# 1 = yes: no new GPS, 2 = yes: new GPS, 3 = no
#This refers to whether the logger was placed in a different...
#...location from the previous logger.
tail(data$key)

data[["key"]] = dplyr::bind_rows(data$key,data.frame(site,logger.id = NA,deployment = NA,date.deployed = NA,date.retrieved = NA, latitude = NA, longitude = NA,retrieval.status = NA,moved = NA))

data$key$deployment[dim(data$key)[1]] = 3
data$key$logger.id[dim(data$key)[1]] = 10602091
data$key$date.deployed[dim(data$key)[1]] = "2015-10-24"
data$key$moved[dim(data$key)[1]] = 1
tail(data$key)

#save(data, file = "HoboDat.RData")

#Add Lat./Long. data from GPS files and manually from PDF maps.
library(rgdal)
setwd("~/Documents/Simon_Fraser_University/PhD_Research/Projects/Data/Original Data/GIS Data/Hobo Loggers/")
Logger.GPS <- readOGR(dsn="Hobo Loggers", layer="WGS84_Hobo_Loggers") 
Logger.GPS = select(data.frame(Logger.GPS),ID_S2014,ID_F2014,time,ele,latitude,longitude)

for(i in 1:nrow(Moved)){
	row = which(data$key$logger.id==Moved$name[i] & data$key$deployment==1)
	data$key$latitude[row]=Moved[i,"Latitude"][[1]]
	data$key$longitude[row]=Moved[i,"Longitude"][[1]]
}

#GPS.Inclusive = filter(data$key,is.na(logger.id)==F)
#setwd("~/Documents/Simon_Fraser_University/PhD_Research/Projects/Data/Original Data/GIS Data/Hobo Loggers")
#write.csv(GPS.Inclusive,"tidyData.csv")

#Add Temperature Data
setwd("~/Documents/Simon_Fraser_University/PhD_Research/Projects/Data/Original Data/HOBOware CSV")

add.data = function(FileName){
	print(FileName)
	browser()
	deploy = NULL
	TData =	read.csv(FileName,header=T, skip=1, colClasses=c("NULL","character","numeric",rep("NULL",count.fields(FileName,sep =",")[2]-3)))
	names(TData) = c("date","temperature")
	TData$date = strptime(TData$date, "%m/%d/%y %I:%M:%S %p")
	
	id = as.numeric(substr(FileName,1,8))
	year = substr(FileName,nchar(FileName)-5,nchar(FileName)-4)
	month = substr(FileName,9,nchar(FileName)-8)
	
	if(month == "oct" & year == "14") deploy = 1
	if((month == "sept"|month == "oct") & year == "15") deploy = 2
	
	site = dplyr::filter(data$key, logger.id == id, deployment == deploy)$site
	
	data[[paste("s",site,sep="")]][[paste("d",deploy,sep="")]] <<- TData
}

for(i in "10602134oct2415.csv"){
	add.data(i)
}

setwd("~/Documents/Simon_Fraser_University/PhD_Research/Projects/River-Network-Temperature")
#save(data, file = "HoboDat.RData")

mixedsort(names(data))
Site = data$s42
plot.dat = NULL
for(i in names(Site)){
	plot.dat = rbind(plyr::mutate(Site[[i]],deployment = i),plot.dat)
}
#plot.dat = subset(plot.dat, date > "2015-09-05")
ggplot(plot.dat, aes(x = date, y = temperature))+
	geom_line(aes(colour = as.factor(deployment))) + 
	theme_classic()

dplyr::filter(data$key, retrieval.status == 4|5|6, deployment == 2)

no.year = strptime(format(plot.dat$date, "2015-%m-%d %H:%M:%S"),format = "%Y-%m-%d %H:%M:%S")
dat = data.frame(plot.dat, no.year)
ggplot(dat, aes(x = no.year, y = temperature))+
	geom_line(aes(colour = as.factor(deployment))) + 
	theme_classic()
