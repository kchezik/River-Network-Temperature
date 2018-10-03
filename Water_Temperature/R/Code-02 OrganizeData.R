#Get necessary libraries.
library(tidyverse); library(rstan); library(lubridate);library(rgdal)
#Initiate cores and make sure unchanged models are not recompiled.
options(mc.cores = parallel::detectCores(), auto_write = TRUE)

#Functions
n_init = function(time){
	library(lubridate)
	if(length(time) == 1) cc = 365
	else {
		dur = data.frame(table(as.numeric(as.duration(time[2:length(time)]-time[1:(length(time)-1)])))) %>% 
			arrange(desc(Freq)) %>% as.tbl(.)
		dur = as.numeric(as.character(dur[[1]]))[1]
		cc = as.numeric(dyears(1)/dur[1])
	}
	if(leap_year(year(time[1]))) cc = cc + 1
	else cc = cc
	cc
}

#Read in cleaned data
df = read_rds("~/Documents/rTDataScrub/Data/thompson_daily_water_manual_clean.rds") 
#Read in water parameter estimates
est = read_csv("~/Documents/rTDataScrub/Data/water_param_est.csv")
#Read in climate data
df_inits = read_rds("~/Documents/rTDataScrub/Data/thompson_air_inits.rds") %>% 
	group_by(site) %>% summarise(air_mean = mean(air_mean), air_A = mean(air_A))

# Test Single Clean Data
clean = df %>% filter(error_lab == "water", temperature>0) %>% 
	mutate(year = year(date),
				 doy = yday(date),
				 year_no = year - (min(year)-1)) %>%
	group_by(site, year) %>% mutate(n = n_init(date))

ref = data.frame(site = unique(clean$site), newID = c(1:length(unique(clean$site))))
clean = left_join(clean,ref)

#Compile Model
compile = stan_model("./Water_Temperature/Stan/Temperature_Model.stan")
#Run Model
mod = sampling(compile, data= list(N = nrow(clean),
																	 Y = length(unique(clean$year)),
																	 R = length(unique(clean$site)),
																	 site = clean$newID,
																	 year = clean$year_no,
																	 y = round(clean$temperature,2),
																	 d = clean$doy,
																	 n = clean$n), 
							 iter = 500, chains = 4,
							 control = list(adapt_delta = .8, max_treedepth = 10))

mod = read_rds("./Water_Temperature/Data/result_250.rds")

#Create full date range.
dates_full = data.frame(day = seq.Date(from = ymd("2014-01-01"), to = ymd("2017-12-31"), by = "day"))

#Incorporate full range of dates into each site for fitting and plotting.
full = clean %>% group_by(site) %>% mutate(day = as.Date(date)) %>% 
	do({
		temp = left_join(dates_full,.,by="day")
		temp$d = yday(temp$day)
		temp$doy = yday(temp$day)
		temp$n = sapply(X = temp$day, FUN = n_init)
		temp$year = year(temp$day)
		temp$year_no = temp$year - (min(temp$year)-1)
		temp$site = unique(.$site)
		temp$newID = unique(.$newID)
		temp
	})

#Extract the parameters
cf = extract(mod)
#Temperature posteriors
rtrn = full %>% group_by(site, year_no, day) %>% do({
	temp = (cf$alpha[,.$newID,.$year_no] + cf$A[,.$newID,.$year_no]*cos(2*pi*.$doy/.$n + cf$tau[,.$newID]*pi)) + cf$S[,.$newID,.$year_no]*cos(2*pi*.$doy/(.$n/2) + cf$tao[,.$newID]*pi)
	Q = quantile(temp, probs = c(.025,.1,.5,.9,.975))
	data.frame(t(Q))
})
#Correct names
names(rtrn) = c("site","year_no","day","Q2.5","Q1","Q5","Q9","Q97")
rtrn = left_join(full,rtrn)

#Determine which sites and years should be omited from SSN analysis.
rtrn %>% filter(newID %in% c(26,76,84,86,92)) %>% 
ggplot(., aes(day, temperature, group = year_no)) + geom_point(size = 0.5) + 
	#geom_point(aes(day, Q5), color = "red") +
	geom_ribbon(aes(x = day, ymin = Q2.5, ymax = Q97), fill = "red", alpha = .7) + 
	facet_wrap(~newID)

#Manually list sites and years to be omitted.
omit = data.frame(site = c(4,5,5,6,10,11,12,13,15,16,17,17,19,19,20,20,21,21,22,22,23,23,
														25,26,26,26,27,27,28,28,29,29,30,31,31,32,33,33,34,37,38,38,
														39,43,43,44,44,46,47,47,50,50,51,55,55,56,56,57,58,58,59,60,61,62,64,64,
														65,66,66,67,69,70,71,71,73,74,75,75,76,76,76,77,78,78,79,79,80,82,82,
														83,83,83,83,84,84,84,85,86,86,86,89,91,91,92,92,92,92,93,93,93,96,96,
														97,100,101,102,103,103),
									year_no=c(1,1,4,1, 3, 4, 3, 1, 2, 1, 1, 4, 3, 4, 3, 4, 3, 4, 1, 4, 1, 3,
														 3, 1, 3, 4, 1, 4, 1, 2, 1, 4, 1, 1, 2, 1, 1, 3, 2, 4, 1, 4,
														 4, 3, 4, 1, 3, 2, 1, 3, 3, 4, 3, 1, 2, 1, 3, 3, 1, 3, 1, 1, 1, 1 , 1, 3,
														 1, 1, 3, 3, 3, 4, 1, 4, 3, 1, 1, 4, 1, 3, 4, 3, 1, 4, 3, 4, 4, 3 ,4,
														 1, 2, 3, 4, 1, 2, 4, 3, 1, 3, 4, 4, 3, 4, 1, 2, 3, 4, 1, 2, 3, 3, 4,
														 4,  3,  3,  1,  1,  2))

#Create a data-frame of sites and years as observations for the ssn.
sites = rtrn %>% select(site,year_no) %>% distinct() %>% dplyr::setdiff(.,omit)
sites = rtrn %>% select(site,newID,year,year_no) %>% distinct() %>% left_join(sites,.)
#Read in HOBO Logger lat/long information.
hobo = readOGR(dsn = "~/sfuvault/Data/Original_Data/GIS_Data/Hobo_Loggers/", layer = "sites")
gps = data.frame(site = as.numeric(stringr::str_extract(hobo@data$Name,"\\d+")),
														 AlbersE = hobo@coords[,1],
														 AlbersN = hobo@coords[,2])
sites = left_join(sites,gps, by = "site")
sp = SpatialPoints(coords = sites[,c("AlbersE","AlbersN")], proj4string = CRS("+init=epsg:3005"))
sp = SpatialPointsDataFrame(coords = sp@coords, data = sites %>% select(1:4))
#Write the observation sites for the SSN.
writeOGR(obj = sp, dsn = "./Water_Temperature/lsn/sites/", layer = "obs",
				 driver = "ESRI Shapefile", overwrite_layer = T)


#Organize table of cost rca attributes

#Elevation
area = read_csv("~/sfu/River_Network_Temperature_Output/predictors/csv_stats/DEM_c.csv", col_names = F) %>% rename(gridcode = X1, count = X2) %>% mutate(rcaArea_KM2 = count*0.025^2) %>% arrange(gridcode) %>% 
	dplyr::select(-count)
el = read_csv("~/sfu/River_Network_Temperature_Output/predictors/csv_stats/DEM_m.csv", col_names = F) %>% 
	rename(gridcode = X1, mu_el = X2) %>% mutate(el_w = mu_el*area$rcaArea_KM2) %>% arrange(gridcode) %>%
	dplyr::select(-mu_el, -gridcode)

#Climate
dat1 = read_csv("~/sfu/River_Network_Temperature_Output/predictors/csv_stats/A_2014.csv", col_names = F) %>% rename(gridcode = X1, amp = X2) %>% mutate(year = 2014, amp_w = amp*area$rcaArea_KM2) %>% 
	dplyr::select(-amp)
dat2 = read_csv("~/sfu/River_Network_Temperature_Output/predictors/csv_stats/A_2015.csv", col_names = F) %>% rename(gridcode = X1, amp = X2) %>% mutate(year = 2015, amp_w = amp*area$rcaArea_KM2) %>% 
	dplyr::select(-amp)
dat3 = read_csv("~/sfu/River_Network_Temperature_Output/predictors/csv_stats/A_2016.csv", col_names = F) %>% rename(gridcode = X1, amp = X2) %>% mutate(year = 2016, amp_w = amp*area$rcaArea_KM2) %>% 
	dplyr::select(-amp)
dat4 = read_csv("~/sfu/River_Network_Temperature_Output/predictors/csv_stats/A_2017.csv", col_names = F) %>% rename(gridcode = X1, amp = X2) %>% mutate(year = 2017, amp_w = amp*area$rcaArea_KM2) %>% 
	dplyr::select(-amp)

amp = bind_rows(dat1,dat2,dat3,dat4) %>% spread(key = year, value = amp_w) %>% 
	rename(Aw2014=`2014`,Aw2015=`2015`,Aw2016=`2016`,Aw2017=`2017`) %>% 
	arrange(gridcode) %>% dplyr::select(-gridcode)

dat1 = read_csv("~/sfu/River_Network_Temperature_Output/predictors/csv_stats/CYT2014_MAT.csv",
								col_names = F) %>% rename(gridcode = X1, alpha = X2) %>%
	mutate(year = 2014, alpha_w = alpha*area$rcaArea_KM2) %>% dplyr::select(-alpha)
dat2 = read_csv("~/sfu/River_Network_Temperature_Output/predictors/csv_stats/CYT2015_MAT.csv",
								col_names = F) %>% rename(gridcode = X1, alpha = X2) %>% 
	mutate(year = 2015, alpha_w = alpha*area$rcaArea_KM2) %>% dplyr::select(-alpha)
dat3 = read_csv("~/sfu/River_Network_Temperature_Output/predictors/csv_stats/CYT2016_MAT.csv",
								col_names = F) %>% rename(gridcode = X1, alpha = X2) %>% 
	mutate(year = 2016, alpha_w = alpha*area$rcaArea_KM2) %>% dplyr::select(-alpha)
dat4 = read_csv("~/sfu/River_Network_Temperature_Output/predictors/csv_stats/CYT2017_MAT.csv",
								col_names = F) %>% rename(gridcode = X1, alpha = X2) %>% 
	mutate(year = 2017, alpha_w = alpha*area$rcaArea_KM2) %>% dplyr::select(-alpha)

alpha = bind_rows(dat1,dat2,dat3,dat4) %>% spread(key = year, value = alpha_w) %>% 
	rename(MuW2014=`2014`,MuW2015=`2015`,MuW2016=`2016`,MuW2017=`2017`) %>%
	arrange(gridcode) %>% dplyr::select(-gridcode)


#Glacier
dat = read_csv("~/sfu/River_Network_Temperature_Output/predictors/csv_stats/glaciers.csv", col_names = F)
dat = dat %>% rename("gridcode" = "X1", "glacier_c" = "X2")
glacier = dat %>% mutate(glacier_km2 = 0.025^2*glacier_c) %>% dplyr::select(-glacier_c) %>%
	arrange(gridcode) %>% dplyr::select(-gridcode)

#Lakes
dat = read_csv("~/sfu/River_Network_Temperature_Output/predictors/csv_stats/lakes.csv", col_names = F)
dat = dat %>% rename("gridcode" = "X1", "lake_c" = "X2")
lake = dat %>% mutate(lake_km2 = 0.025^2*lake_c) %>% dplyr::select(-lake_c) %>%
	arrange(gridcode) %>% dplyr::select(-gridcode)

#Logging
dat1 = read_csv("~/sfu/River_Network_Temperature_Output/predictors/csv_stats/logging_year_1_zonal.csv")
names(dat1) = tolower(names(dat1)); dat1 = dat1 %>% rename(logg1_c = count)
dat2 = read_csv("~/sfu/River_Network_Temperature_Output/predictors/csv_stats/logging_year_2_zonal.csv") %>% gather(key = year, value = logg2_c, 2:27) %>% rename(gridcode = GRIDCODE)
dat2$year = as.numeric(dat2$year)

dat = left_join(dat2,dat1) %>% mutate(logg_c = logg2_c+logg1_c)

logging = dat %>% mutate(logg_km2 = 0.025^2*logg_c) %>% 
	group_by(gridcode) %>% 
	summarise(logg_km2_5 = sum(logg_km2[year>2005]),
						logg_km2_10 = sum(logg_km2[year>2000]),
						logg_km2_cum = sum(logg_km2)) %>% 
	arrange(gridcode) %>% dplyr::select(-gridcode)
#dat_year = dat %>% select(1,2,6) %>% spread(key = year, value = logg_km2)
#names(dat_year)[2:27] = paste("logg",names(dat_year)[2:27], sep = "_")
#dat = dat %>% select(-c(2:6)) %>% distinct() %>% left_join(.,dat_year)

#Wildfire
dat1 = read_csv("~/sfu/River_Network_Temperature_Output/predictors/csv_stats/wildfire_year_1_zonal.csv") %>% gather(key = year, value = fire1_c, 2:27) %>% rename(gridcode = GRIDCODE)
dat1$year = as.numeric(dat1$year)
dat2 = read_csv("~/sfu/River_Network_Temperature_Output/predictors/csv_stats/wildfire_year_2_zonal.csv") %>% gather(key = year, value = fire2_c, 2:27) %>% rename(gridcode = GRIDCODE)
dat2$year = as.numeric(dat2$year)

dat = left_join(dat2,dat1) %>% mutate(fire_c = fire2_c+fire1_c)

fire = dat %>% mutate(fire_km2 = 0.025^2*fire_c) %>% 
	group_by(gridcode) %>% 
	summarise(fire_km2_5 = sum(fire_km2[year>2005]),
						fire_km2_10 = sum(fire_km2[year>2000]),
						fire_km2_cum = sum(fire_km2)) %>% 
	arrange(gridcode) %>% dplyr::select(-gridcode)

#Compile full predictor dataset
final = bind_cols(area,el,amp,alpha,lake,glacier,logging,fire)
write_csv(final, "~/sfu/River_Network_Temperature_Output/predictors/csv_stats/full_preds.csv")


#Add new site for each year to prediction dataset.
preds = readOGR("./Water_Temperature/lsn/lsn.ssn/", layer = "preds")
df = preds@data %>% mutate(year = 2014)
df2 = df; df2$year = 2015
df3 = df; df3$year = 2016
df4 = df; df4$year = 2017
dat = bind_rows(df,df2,df3,df4)

sp = SpatialPointsDataFrame(coords = dat[,c("NEAR_X","NEAR_Y")], data = dat, 
														proj4string = preds@proj4string)
writeOGR(obj = sp, dsn = "./Water_Temperature/lsn/preds/", layer = "preds", driver = "ESRI Shapefile")
