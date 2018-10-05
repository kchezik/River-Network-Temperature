library(tidyverse); library(SSN); library(rgdal); library(boot); library(lubridate);
library(foreach); library(doParallel) #Activate Libraries
#Read in the .ssn network data.
network = importSSN(filepath = "~/sfuvault/River-Network-Temperature/Water_Temperature/lsn/lsn.ssn/",
										predpts = "preds", o.write = F) 
#Gather dataframes from observations and predictions for manipulation.
obs = getSSNdata.frame(network, Name = "Obs"); names(obs)[4] = "year"
preds = getSSNdata.frame(network, Name = "preds"); names(preds)[15] = "year"
#Manual correction for points within lakes and no associated costRCA. For improvement, we could/should eliminate prediction sites that coincide with lakes. This would account for thermal refuge offered by the lake.
preds[preds$rid == 1645,"rid"] = 1769
preds[preds$rid == 3139,"rid"] = 3284
#Gather stream network dataframe and correct column names for consistancy in subsequent analysis.
netDat = network@data; names(netDat)[c(5,36,23)] = c("rca_h2o","El_km2","rca_km2")

#Create the distance matrices for the model. Only needs to be done once. Look for the 'distance' folder in the lsn.ssn folder.
#createDistMat(network, predpts = "preds", o.write = F, amongpreds = F)

#Organize predictors

#Calculate climate variables for each site
Climate_Year = function(netID,rid,year,c_var){
	nms = names(netDat); srch = paste(c_var,'\\d+',sep=""); 
	var_nms = nms[!is.na(str_match(nms,srch)[,1])]
	col = var_nms[as.numeric(str_match(var_nms,'\\d+')) == unique(year-2000)]
	netDat[(netDat$netID == unique(netID) & netDat$rid == unique(rid)),col]
}
#Calculate mean annual air temperature and range for the contributing watershed for each site and year.
obs_org = obs %>% group_by(netID,rid,year) %>% mutate(air_mean = Climate_Year(netID,rid,year,"Mu"),
																						 air_A = Climate_Year(netID,rid,year,"A"))
preds_org = preds %>% group_by(netID,rid,year) %>% 
	mutate(air_mean = Climate_Year(netID,rid,year,"Mu"),
				 air_A = Climate_Year(netID,rid,year,"A"))


#Function to calculate site watershed attributes
WA = function(netID, rid, ratio, var, perc = T){
	#Isolate
	edge = netDat[netDat$netID == unique(netID) & netDat$rid == unique(rid),]
	nms = names(edge); U_RCA = nms[!is.na(str_match(nms,var))[,1]]
	h2o = U_RCA[grep("h2",U_RCA)]
	rca = U_RCA[grep("km",U_RCA)]
	#Calculate
	area = edge[,"rca_h2o"] + (1-unique(ratio))*edge[,"rca_km2"]
	effect = edge[,h2o] + (1-unique(ratio))*edge[,rca]
	#Return Percentage
	if(perc){
		rtrn = rep((effect/area),length(ratio))
		if(any(rtrn==0)) rtrn = rtrn + 0.001
	}
	else rtrn = rep(effect,length(ratio))
	return(rtrn)
}
#Function to calculate rca sums and weighted means.
RCA = function(netID, rid, var, cost = T, perc = F, WA = F){
	#Isolate
	edge = netDat[netDat$netID == unique(netID) & netDat$rid == unique(rid),]
	nms = names(edge); U_RCA = nms[!is.na(str_match(nms,var))[,1]]
	h2o = U_RCA[grep("h2",U_RCA)]
	rca = U_RCA[grep("km",U_RCA)]
	if(perc == T){
		if(cost == T) effect = edge[,rca]/edge[,"rca_km2"]
		else effect = edge[,h2o]/edge[,"rca_h2o"]
	}
	else{
		if(cost == T) effect = edge[,rca]
		else effect = edge[,h2o]
	}
	if(WA == T){
		effect = (edge[,rca]*edge[,"rca_km2"] + edge[,h2o]*edge[,"rca_h2o"])/
			(edge[,"rca_km2"]+edge[,"rca_h2o"])
	}
	return(effect)
}
#Summarise the watershed attributes contributing to observation sites.
obs_org = obs_org %>% group_by(netID,rid,ratio) %>% 
	mutate(ca_h2o = WA(netID, rid, ratio, var = "rca", perc = F),
				 ele_h2o = RCA(netID, rid, var = "El", cost = F, WA = T),
				 lake_h2o = WA(netID, rid, ratio, var = "lake"),
				 glacier_h2o = WA(netID, rid, ratio, var = "glacier"),
				 logcum_h2o = WA(netID, rid, ratio, var = "log_cu|log\\w{5}cu"),
				 firecum_h2o = WA(netID, rid, ratio, var = "fire_cu|fire\\w{5}c")
				 ) %>% ungroup() %>% 
	mutate(site_fac = as.factor(newID),
				 year_fac = as.factor(year),
				 scl_air_mn = (air_mean-mean(air_mean))/sd(air_mean),
				 scl_air_A = (air_A-mean(air_A))/sd(air_A),
				 scl_ca = (log(ca_h2o)-mean(log(ca_h2o)))/sd(log(ca_h2o)),
				 scl_ele = (ele_h2o-mean(ele_h2o))/sd(ele_h2o),
				 lgt_lake = logit(lake_h2o),
				 lgt_glacier = logit(glacier_h2o),
				 lgt_logcum = logit(logcum_h2o),
				 lgt_firecum = logit(firecum_h2o))
#Summarise the watershed attributes contributing to prediction sites.
preds_org = preds_org %>% group_by(netID,rid,ratio) %>% 
	mutate(ca_h2o = WA(netID, rid, ratio, var = "rca", perc = F),
				 ele_h2o = RCA(netID, rid, var = "El", cost = F, WA = T),
				 lake_h2o = WA(netID, rid, ratio, var = "lake"),
				 glacier_h2o = WA(netID, rid, ratio, var = "glacier"),
				 logcum_h2o = WA(netID, rid, ratio, var = "log_cu|log\\w{5}cu"),
				 firecum_h2o = WA(netID, rid, ratio, var = "fire_cu|fire\\w{5}c")
				 ) %>% ungroup() %>% 
	mutate(site_fac = locID,
				 year_fac = as.factor(year),
				 scl_air_mn = (air_mean-mean(air_mean))/sd(air_mean),
				 scl_air_A = (air_A-mean(air_A))/sd(air_A),
				 scl_ca = (ca_h2o-mean(ca_h2o))/sd(ca_h2o),
				 scl_ele = (ele_h2o-mean(ele_h2o))/sd(ele_h2o),
				 lgt_lake = logit(lake_h2o),
				 lgt_glacier = logit(glacier_h2o),
				 lgt_logcum = logit(logcum_h2o),
				 lgt_firecum = logit(firecum_h2o))
#Put prediction data back into the network.
	#Must make sure the returned dataframe is not a tibble. 
	#Also, make sure the original rownames are present.
preds_org = data.frame(preds_org)
row.names(preds_org) = row.names(preds)
network <- putSSNdata.frame(preds_org, network, Name = "preds")

#Sample observations.

#Read in Model
mod = read_rds("./Water_Temperature/Data/result_250.rds")
#Aggregate Coefficients
cf = rstan::extract(mod)

#Sample coefficients for each site
cf_samp = function(ID,year,var,row,Y=TRUE){
	v = cf[[var]]
	if(Y)	v[row,ID,year]
	else	v[row,ID]
}
	


#-------------------------> START Sequence <-------------------------#
cl <- makeCluster(4)
registerDoParallel(cl)

rtrn = foreach(i=1:4, .packages = c('tidyverse','SSN','boot'), .errorhandling = "remove") %dopar% {
	
	#Sample a unique row and sampe from each coefs. posterior.
	row = sample(x = c(1:1000), size = 1)
	obs_fit = obs_org %>% group_by(newID, year_no) %>% 
		mutate(alpha = cf_samp(newID, year_no, "alpha", row),
					 A = cf_samp(newID, year_no, "A", row),
					 S = cf_samp(newID, year_no, "S", row),
					 tau = cf_samp(newID, year_no, "tau", row, Y=F),
					 tao = cf_samp(newID, year_no, "tao", row, Y=F))
	
	#Must make sure the returned dataframe is not a tibble.
	obs_fit = data.frame(obs_fit)
	network <- putSSNdata.frame(obs_fit, network, Name = "Obs")
	
	#Models
	
	#Alpha
	#Colder places have been logged less between 1985-2010 than warmer places. Easier to log low down in the watershed? Is the negative effect of logging due to the temporal lag? Maybe, logged places during this time period are now quite forested and 'unforested' places accoring to this dataset have since been logged and are therefore warmer.
	
	mod_alpha = glmssn(formula = alpha~scl_air_mn+scl_air_A+
										 	lgt_firecum*scl_ca+lgt_glacier*scl_ca+lgt_logcum+
										 	lgt_lake*scl_ele,
										 ssn.object = network, EstMeth = "REML", family = "Gaussian", 
										 CorModels = c("site_fac","year_fac","Exponential.tailup"), 
										 addfunccol = "afvArea")
	#Amplitude
	mod_A = glmssn(formula = A~scl_air_mn+scl_air_A+lgt_glacier+scl_ca+lgt_lake,
								 ssn.object = network, EstMeth = "REML", family = "Gaussian", 
								 CorModels = c("year_fac","Exponential.tailup"), 
								 addfunccol = "afvArea")
	#Snow
	mod_snow = glmssn(formula = S~scl_air_mn+scl_ele+lgt_glacier,
										ssn.object = network, EstMeth = "REML", family = "Gaussian", 
										CorModels = c("year_fac","Exponential.tailup"), 
										addfunccol = "afvArea")
	#Annual tau
	mod_tau = glmssn(formula = tau~scl_ele,
									 ssn.object = network, EstMeth = "REML", family = "Gaussian", 
									 CorModels = c("year_fac","Exponential.tailup"), 
									 addfunccol = "afvArea")
	#Season tao
	mod_tao = glmssn(formula = tao~ele_h2o,
									 ssn.object = network, EstMeth = "REML", family = "Gaussian", 
									 CorModels = c("year_fac","Exponential.tailup"), 
									 addfunccol = "afvArea")
	
	#Predict alpha and save to matrix.
	pred_alpha <- predict.glmssn(mod_alpha, predpointsID = "preds")
	alpha_mat = pred_alpha$ssn.object@predpoints@SSNPoints[[1]]@point.data$alpha
	alpha_coef = summary(mod_alpha)
	
	#Predict A and save to matrix.
	pred_A <- predict.glmssn(mod_A, predpointsID = "preds")
	A_mat = pred_A$ssn.object@predpoints@SSNPoints[[1]]@point.data$A
	A_coef = summary(mod_A)
	
	#Predict snow and save to matrix.
	pred_snow <- predict.glmssn(mod_snow, predpointsID = "preds")
	snow_mat = pred_snow$ssn.object@predpoints@SSNPoints[[1]]@point.data$S
	snow_coef = summary(mod_snow)
	
	#Predict tau and save to matrix.
	pred_tau <- predict.glmssn(mod_tau, predpointsID = "preds")
	tau_mat = pred_tau$ssn.object@predpoints@SSNPoints[[1]]@point.data$tau
	tau_coef = summary(mod_tau)
	
	#Predict tao and save to matrix.
	pred_tao <- predict.glmssn(mod_tao, predpointsID = "preds")
	tao_mat = pred_tao$ssn.object@predpoints@SSNPoints[[1]]@point.data$tao
	tao_coef = summary(mod_tao)
	
	write_rds(list(estimates = data.frame(alpha = alpha_mat, A = A_mat, snow = snow_mat,
																				tau = tau_mat, tao = tao_mat),
								 alpha = alpha_coef$fixed.effects.estimates, A = A_coef$fixed.effects.estimates,
								 snow = snow_coef$fixed.effects.estimates,
								 tau = tau_coef$fixed.effects.estimates, tao = tao_coef$fixed.effects.estimates),
						path = paste("./Water_Temperature/Data/iterations/",as.character(i),".rds",sep = ""))
	
	list(estimates = data.frame(alpha = alpha_mat, A = A_mat, snow = snow_mat,
															tau = tau_mat, tao = tao_mat),
			 alpha = alpha_coef$fixed.effects.estimates, A = A_coef$fixed.effects.estimates,
			 snow = snow_coef$fixed.effects.estimates,
			 tau = tau_coef$fixed.effects.estimates, tao = tao_coef$fixed.effects.estimates)
}
write_rds(rtrn, "./Water_Temperature/Data/network_4.rds")
stopCluster(cl)
#-------------------------> END Sequence <-------------------------#

#Visual Diagnostics
cv.out = CrossValidationSSN(mod_tao)
par(mfrow = c(1, 2))
plot(mod_tao$sampinfo$z,
		 cv.out[, "cv.pred"], pch = 19,
		 xlab = "Observed Data", ylab = "LOOCV Prediction")
abline(0, 1)
plot( na.omit( getSSNdata.frame(network)[, "tao"]),
			cv.out[, "cv.se"], pch = 19,
			xlab = "Observed Data", ylab = "LOOCV Prediction SE")

plot(network, "A", lwdLineCol = "afvArea", lwdLineEx = 10, lineCol = "black", pch = 19, xlab = "x-coordinate (m)", ylab = "y-coordinate (m)", asp   =   1)

#The notable interaction of lake and elevation may be explained by the fact that higher elevations have steeper slopes and fewer lakes and therefore a negative relationship with temperature. This is supported by Lisi:2015.


# Thermal Stress Probability Accumulation Analysis

#Sample early (June), mid (July) and late (August) migration dates.
samp_dates = function(single=F){
	library(lubridate); library(tidyverse)
	#Dates by early, mid and late summer migration
	dat = data.frame(migrate = c("early","mid","late"),
						 date = c(ymd("2000-06-01"), ymd("2000-07-01"), ymd("2000-08-01")))
	#Minimum and maximum
	if(single == F){
		dat = dat %>% 
			mutate(year = year(date),
						 min = yday(date), max = yday(date)+as.numeric(days_in_month(date)),
						 start = round(runif(3,min,max),0)) %>% 
			select(migrate,start)
	} else {
		dat = data.frame(migrate = 
										 	c(rep(as.character(dat$migrate[1]),as.numeric(days_in_month(dat$date[1]))),
		rep(as.character(dat$migrate[2]),as.numeric(days_in_month(dat$date[2]))),
		rep(as.character(dat$migrate[3]),as.numeric(days_in_month(dat$date[3])))),
		start = 
			c(seq(yday(dat$date[1]),yday(dat$date[1])+(days_in_month(dat$date[1])-1), by = 1),
			seq(yday(dat$date[2]),yday(dat$date[2])+(days_in_month(dat$date[2])-1), by = 1),
			seq(yday(dat$date[3]),yday(dat$date[3])+(days_in_month(dat$date[3])-1), by = 1)))
	}
	return(dat)
}


cl <- makeCluster(7)
registerDoParallel(cl)
#Sample an early, mid and late summer migration date and calculate temperatures.
probs = plyr::dlply(.data = samp_dates(single=F), .variables = "migrate", .fun = function(x){
	
	#----->Data<-----#
	#SSN Model fits
	fits = read_rds("./Water_Temperature/Data/network_250.rds")
	#Network binary and rid IDs.
	net1 = read_delim("./Water_Temperature/lsn/lsn.ssn/netID1.dat", delim = ",",
										col_names = c("rid","binaryID"), skip = 1, col_types = cols("i","c"))
	net1=as.data.frame(net1)
	net2 = read_delim("./Water_Temperature/lsn/lsn.ssn/netID2.dat", delim = ",",
										col_names = c("rid","binaryID"), skip = 1, col_types = cols("i","c"))
	net2=as.data.frame(net2)
	#Prediction Points
	preds = getSSNdata.frame(network, Name = "preds"); names(preds)[15] = "year"
	preds[preds$rid == 1645,"rid"] = 1769
	preds[preds$rid == 3139,"rid"] = 3284
	
	#----->Functions<-----#
	#Determine downstream nodes.
	dwnStrm = function(net, binary){
		rows = vector(mode = "numeric", length = nchar(binary))
		for(i in 1:nchar(binary)){
			rows[i] = which(net[,2]==substr(binary,start = 1, stop = i))
		}
		return(net[rows,1])
	}
	#Calculate estimated d given a migration start date, distance from outlet and ...
	#...standard travel rate. 
	swim_date = function(start_dates, distance){
		to_thomp = 447; up_thomp = distance/1000
		round(start_dates + (to_thomp+up_thomp)/36,0)
	}
	#Caluculate mean daily temperatures based on estimated d
	temperature = function(df){
		(df$alpha + df$A*cos(2*pi*df$d/df$n+pi*df$tau)) + cos(2*pi*df$d/(df$n/2)+pi*df$tao)
	}
	
	#----->Analysis<-----#
	plyr::dlply(.data = x, .variables = "start", .fun = function(y){
		tmp = t(plyr::ldply(fits, function(z){
			n = if_else(leap_year(ymd(paste(as.character(preds$year),"01","01",sep="-"))),366,365)
			est = mutate(z$estimates, d = swim_date(y$start,preds$upDist),n = n)
			temperature(est)
		})) %>% data.frame(.,row.names = NULL)
		preds %>% group_by(locID) %>% do({
			if(.$netID==1) net = net1
			if(.$netID==2) net = net2
			binary = net[net$rid==.$rid,"binaryID"]
			#Determine downstream prediction points.
			dwn_rid = dwnStrm(net, binary)
			#Probability of >18ºC experience
			pts = tmp[which(preds$year==.$year & preds$upDist<=.$upDist & preds$rid %in% dwn_rid),]
			G19_sum = sum(pts>19);G20_sum = sum(pts>20);G21_sum = sum(pts>21);G22_sum = sum(pts>22)
			G23_sum = sum(pts>23);G24_sum = sum(pts>24); denom = length(pts)*nrow(pts)
			data.frame(Sum19 = G19_sum, G19 = G19_sum/denom,
								 Sum20 = G20_sum, G20 = G20_sum/denom,
								 Sum21 = G21_sum, G21 = G21_sum/denom,
								 Sum22 = G22_sum, G22 = G22_sum/denom,
								 Sum23 = G23_sum, G23 = G23_sum/denom,
								 Sum24 = G24_sum, G24 = G24_sum/denom)
		})
	}, .parallel = F, .paropts = list(.packages = c('tidyverse','lubridate'),
																		.export = c("fits","preds","net1","net2",
																								"swim_date","temperature",
																								"samp_dates","dwnStrm")))
})
write_rds(probs, "./Water_Temperature/Data/probs_flatD.rds")
stopCluster(cl)

#The approximate distance from the Albion test fishery (Parken et al. 2008) to the Thompson as Kamloops is 447km. Adding this initial distance to the distance upstream from Kamloops give the total travel distance. We can then divide by an approximate travel rate of 36km/day (Salinger and Anderson 2006).

probs = read_rds("./Water_Temperature/Data/probs_flatD.rds")
#Summarize by migration group and place in shapefile for mapping.
lim_summary = function(probs, cols, prefix){
	plyr::llply(probs, .fun=function(x){
		tmp = plyr::ldply(x, .fun=function(y){
			data.frame(locID = y$locID, cum = y[[cols[1]]], probs = y[[cols[2]]], year = preds$year)
		})
		#Get the mean by year accross days.
		mean_c = tmp %>% group_by(year,locID) %>% 
			summarise(mean_p = mean(probs), mean_c = round(mean(cum),0))
		#Get the day with the greatest average value.
		max_p = tmp %>% group_by(year,start) %>% summarise(mean_p = mean(probs)) %>% 
			filter(mean_p == max(mean_p)) %>% rename(start_p = start) %>% 
			inner_join(.,tmp,by = c("year","start_p"="start")) %>% select(-mean_p, -cum)
		max_c = tmp %>% group_by(year,start) %>% summarise(mean_c = mean(cum)) %>% 
			filter(mean_c == max(mean_c))  %>% rename(start_c = start) %>% 
			inner_join(.,tmp,by = c("year","start_c"="start")) %>% select(-mean_c, -probs)
		max = left_join(max_p, max_c) 
		#return summarized data by migration.
		final = max %>% left_join(., mean_c, by = c("year","locID")) %>% select(3,1,7,8,2,4,5,6)
		names(final)[c(3:8)] = paste(prefix,names(final)[c(3:8)],sep="_")
		return(final)
	})
}
#Summarize by threshold
G19 = lim_summary(probs, c(2,3), "G19")
G20 = lim_summary(probs, c(4,5), "G20")
G21 = lim_summary(probs, c(6,7), "G21")
G22 = lim_summary(probs, c(8,9), "G22")
G23 = lim_summary(probs, c(10,11), "G23")
G24 = lim_summary(probs, c(12,13), "G24")

#Save shape files

shp_save = function(pred_df, period, yr){
	#Bind together the thresholds
	p_riod = pred_df %>% select(-c(1:7)) %>% 
		left_join(.,G19[[period]]) %>% left_join(.,G20[[period]]) %>%
		left_join(.,G21[[period]]) %>% left_join(.,G22[[period]]) %>%
		left_join(.,G23[[period]]) %>% left_join(.,G24[[period]])
	#Filter out by year
	p_riod = data.frame(network@predpoints@SSNPoints[[1]]@point.coords) %>% 
		bind_cols(., p_riod) %>% filter(year == yr)
	#Make spatial object
	p_riod = SpatialPointsDataFrame(coords = p_riod[,c("coords.x1","coords.x2")], 
												 data = p_riod[,c(3:46)], proj4string = network@proj4string)
	#Write shapefile.
	lyr = paste(period, as.character(yr), sep = "_")
	writeOGR(p_riod, dsn = "./Water_Temperature/lsn/risk_probs/", layer = lyr,
					 overwrite_layer = T, driver = "ESRI Shapefile")
}

#Early
shp_save(preds, "early", 2014)
shp_save(preds, "early", 2015)
shp_save(preds, "early", 2016)
shp_save(preds, "early", 2017)
#Mid
shp_save(preds, "mid", 2014)
shp_save(preds, "mid", 2015)
shp_save(preds, "mid", 2016)
shp_save(preds, "mid", 2017)
#Late
shp_save(preds, "late", 2014)
shp_save(preds, "late", 2015)
shp_save(preds, "late", 2016)
shp_save(preds, "late", 2017)