library(tidyverse); library(SSN); library(rgdal); library(boot); library(lubridate);
library(foreach); library(doParallel) #Activate Libraries
#Read in the .ssn network data.
network = importSSN(filepath = "~/sfuvault/River-Network-Temperature/Water_Temperature/lsn/newLSN.ssn/", predpts = "preds", o.write = F) 
#Gather dataframes from observations and predictions for manipulation.
obs = getSSNdata.frame(network, Name = "Obs"); names(obs)[4] = "year"
preds = getSSNdata.frame(network, Name = "preds"); names(preds)[4] = "year"
#For improvement, we could/should eliminate prediction sites that coincide with lakes. This would account for thermal refuge offered by the lake.

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
obs_org = obs %>% group_by(netID,rid,year) %>%
	mutate(air_mean = Climate_Year(netID,rid,year,"Mu"),
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
	if(area == 0) area = 0.001
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
		num = (edge[,rca]*edge[,"rca_km2"] + edge[,h2o]*edge[,"rca_h2o"])
		denom = (edge[,"rca_km2"]+edge[,"rca_h2o"])
		if(denom == 0) denom = 0.001
		effect = num/denom
	}
	return(effect)
}

#Summarise the watershed attributes contributing to observation sites.
obs_org = obs_org %>% group_by(netID,rid,ratio) %>% 
	mutate(ca_h2o = WA(netID, rid, ratio, var = "rca", perc = F),
				 ele_h2o = RCA(netID, rid, var = "El", cost = F, WA = T),
				 lake_h2o = WA(netID, rid, ratio, var = "lake", perc = T),
				 glacier_h2o = WA(netID, rid, ratio, var = "glacier", perc = T),
				 logcum_h2o = WA(netID, rid, ratio, var = "log_cu|log\\w{5}cu", perc = T),
				 firecum_h2o = WA(netID, rid, ratio, var = "fire_cu|fire\\w{5}c", perc = T))
#Summarise the watershed attributes contributing to prediction sites.
preds_org = preds_org %>% group_by(netID,rid,ratio) %>% 
	mutate(ca_h2o = WA(netID, rid, ratio, var = "rca", perc = F),
				 ele_h2o = RCA(netID, rid, var = "El", cost = F, WA = T),
				 lake_h2o = WA(netID, rid, ratio, var = "lake", perc = T),
				 glacier_h2o = WA(netID, rid, ratio, var = "glacier", perc = T),
				 logcum_h2o = WA(netID, rid, ratio, var = "log_cu|log\\w{5}cu", perc = T),
				 firecum_h2o = WA(netID, rid, ratio, var = "fire_cu|fire\\w{5}c", perc = T))

#Get scaling factors across observation and prediction datasets.
cs = data.frame(vars = names(obs_org)[17:24]) %>% group_by(vars) %>% do({
	nums = c(preds_org[,names(preds_org)==.$vars][[1]],obs_org[,names(obs_org)==.$vars][[1]])
	#Ignore zeros when scaling these log transformed values.
	if(.$vars %in% c("ca_h2o","firecum_h2o","glacier_h2o","lake_h2o","logcum_h2o")){
		data.frame(mean = mean(log(nums[nums>0])), sd = sd(log(nums[nums>0])))
	} else data.frame(mean = mean(nums), sd = sd(nums))
})
#Scale and logit transform data.
obs_org = obs_org %>% ungroup() %>% 
	mutate(site_fac = as.factor(newID),
				 year_fac = as.factor(year),
				 scl_air_mn = (air_mean-cs[cs$vars=="air_mean","mean"][[1]])/cs[cs$vars=="air_mean","sd"][[1]],
				 scl_air_A = (air_A-cs[cs$vars=="air_A","mean"][[1]])/cs[cs$vars=="air_A","sd"][[1]],
				 scl_ca = (log(ca_h2o+0.001)-cs[cs$vars=="ca_h2o","mean"][[1]])/cs[cs$vars=="ca_h2o","sd"][[1]],
				 scl_ele = (ele_h2o-cs[cs$vars=="ele_h2o","mean"][[1]])/cs[cs$vars=="ele_h2o","sd"][[1]],
				 lgt_lake = logit(lake_h2o),
				 lgt_glacier = logit(glacier_h2o),
				 lgt_logcum = logit(logcum_h2o),
				 lgt_firecum = logit(firecum_h2o))

preds_org = preds_org %>% ungroup() %>% 
	mutate(site_fac = as.factor(UID),
				 year_fac = as.factor(year),
				 scl_air_mn = (air_mean-cs[cs$vars=="air_mean","mean"][[1]])/cs[cs$vars=="air_mean","sd"][[1]],
				 scl_air_A = (air_A-cs[cs$vars=="air_A","mean"][[1]])/cs[cs$vars=="air_A","sd"][[1]],
				 scl_ca = (log(ca_h2o+0.001)-cs[cs$vars=="ca_h2o","mean"][[1]])/cs[cs$vars=="ca_h2o","sd"][[1]],
				 scl_ele = (ele_h2o-cs[cs$vars=="ele_h2o","mean"][[1]])/cs[cs$vars=="ele_h2o","sd"][[1]],
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

#Get average estimates of temporal coefficients
data.frame(alpha_sd = mean(cf$alpha), A_sd = mean(cf$A), S_sd = mean(cf$S),
					 tau_sd = mean(cf$tau), tao_sd = mean(cf$tao), sigma_sd = mean(cf$sigma))


#-------------------------> START Sequence <-------------------------#
cl <- makeCluster(7)
registerDoParallel(cl)
start_time = Sys.time()
rtrn = foreach(i=1:250, .packages = c('tidyverse','SSN','boot'), .errorhandling = "remove") %dopar% {
	#Sample a unique row and sampe from each coefs. posterior.
	row = sample(x = c(1:1000), size = 1)
	obs_fit = obs_org %>% group_by(newID, year_no) %>% 
		mutate(alpha = cf_samp(newID, year_no, "alpha", row),
					 A = cf_samp(newID, year_no, "A", row),
					 S = cf_samp(newID, year_no, "S", row),
					 tau = cf_samp(newID, year_no, "tau", row, Y=F),
					 tao = cf_samp(newID, year_no, "tao", row, Y=F),
					 sigma = cf_samp(newID, year_no, "sigma", row, Y=F))
	#Standardize response variables
	obs_fit$std_alpha = log(obs_fit$alpha)/sd(log(obs_fit$alpha))
	obs_fit$std_A = log(obs_fit$A)/sd(log(obs_fit$A))
	obs_fit$std_S = obs_fit$S/sd(obs_fit$S)
	obs_fit$std_tau = obs_fit$tau/sd(obs_fit$tau)
	obs_fit$std_tao = obs_fit$tao/sd(obs_fit$tao)
	obs_fit$std_sigma = obs_fit$sigma/sd(obs_fit$sigma)
	
	#Must make sure the returned dataframe is not a tibble.
	obs_fit = data.frame(obs_fit)
	network <- putSSNdata.frame(obs_fit, network, Name = "Obs")
	
	#Models
	
	#Alpha

		#Colder places have been logged less between 1985-2010 than warmer places. Easier to log low down in the watershed? Is the negative effect of logging due to the temporal lag? Maybe, logged places during this time period have since re-grown and are now quite forested leaving 'unforested' places likely to have been more recently logged. After accounting for air and lake effects, does this temporal mismatch in the data result in a negative relationship? If so, maybe this explains the postive effect of glaciers which are typically high up in the watershed where logging may have occurred more recently such that lower logging and higher glaciers are associated with a small warming trend in the water.
	
		#The notable interaction of lake and elevation may be explained by the fact that higher elevations have steeper slopes and fewer lakes and therefore a negative relationship with temperature. This is supported by Lisi:2015.
	mod_alpha = glmssn(formula = std_alpha~scl_air_mn+scl_air_A+
												 	lgt_lake*scl_ele+lgt_firecum,
										 ssn.object = network, EstMeth = "REML", family = "Gaussian", 
										 CorModels = c("site_fac","year_fac","Exponential.tailup"), 
										 addfunccol = "afvArea")
	#Amplitude
	mod_A = glmssn(formula = std_A~scl_air_mn+scl_air_A+lgt_glacier+scl_ca+lgt_lake,
								 ssn.object = network, EstMeth = "REML", family = "Gaussian", 
								 CorModels = c("site_fac","year_fac","Exponential.tailup"), 
								 addfunccol = "afvArea")
	#Snow
	mod_snow = glmssn(formula = std_S~scl_air_mn+scl_ele+lgt_glacier,
										ssn.object = network, EstMeth = "REML", family = "Gaussian", 
										CorModels = c("site_fac","year_fac","Exponential.tailup"), 
										addfunccol = "afvArea")
	#Annual tau
	mod_tau = glmssn(formula = std_tau~scl_ele,
									 ssn.object = network, EstMeth = "REML", family = "Gaussian", 
									 CorModels = c("site_fac","year_fac","Exponential.tailup"), 
									 addfunccol = "afvArea")
	#Season tao
	mod_tao = glmssn(formula = std_tao~scl_ele,
									 ssn.object = network, EstMeth = "REML", family = "Gaussian", 
									 CorModels = c("site_fac","year_fac","Exponential.tailup"), 
									 addfunccol = "afvArea")
	#Variance
	mod_sigma = glmssn(formula = std_sigma~scl_ele,
									 ssn.object = network, EstMeth = "REML", family = "Gaussian", 
									 CorModels = c("site_fac","year_fac","Exponential.tailup"), 
									 addfunccol = "afvArea")
	
	#Predict alpha and save to matrix.
	pred_alpha <- predict.glmssn(mod_alpha, predpointsID = "preds")
	alpha_mat = pred_alpha$ssn.object@predpoints@SSNPoints[[1]]@point.data$std_alpha
	alpha_mat = exp(alpha_mat*sd(log(obs_fit$alpha)))
	alpha_coef = summary(mod_alpha)
	alpha_varcomp = varcomp(mod_alpha)
	
	#Predict A and save to matrix.
	pred_A <- predict.glmssn(mod_A, predpointsID = "preds")
	A_mat = pred_A$ssn.object@predpoints@SSNPoints[[1]]@point.data$std_A
	A_mat = exp(A_mat*sd(log(obs_fit$A)))
	A_coef = summary(mod_A)
	A_varcomp = varcomp(mod_A)
	
	#Predict snow and save to matrix.
	pred_snow <- predict.glmssn(mod_snow, predpointsID = "preds")
	snow_mat = pred_snow$ssn.object@predpoints@SSNPoints[[1]]@point.data$std_S
	snow_mat = snow_mat*sd(obs_fit$S)
	snow_coef = summary(mod_snow)
	snow_varcomp = varcomp(mod_snow)
	
	#Predict tau and save to matrix.
	pred_tau <- predict.glmssn(mod_tau, predpointsID = "preds")
	tau_mat = pred_tau$ssn.object@predpoints@SSNPoints[[1]]@point.data$std_tau
	tau_mat = tau_mat*sd(obs_fit$tau)
	tau_coef = summary(mod_tau)
	tau_varcomp = varcomp(mod_tau)
	
	#Predict tao and save to matrix.
	pred_tao <- predict.glmssn(mod_tao, predpointsID = "preds")
	tao_mat = pred_tao$ssn.object@predpoints@SSNPoints[[1]]@point.data$std_tao
	tao_mat = tao_mat*sd(obs_fit$tao)
	tao_coef = summary(mod_tao)
	tao_varcomp = varcomp(mod_tao)
	
	#Predict sigma and save to matrix.
	pred_sigma <- predict.glmssn(mod_sigma, predpointsID = "preds")
	sigma_mat = pred_sigma$ssn.object@predpoints@SSNPoints[[1]]@point.data$std_sigma
	sigma_mat = sigma_mat*sd(obs_fit$sigma)
	sigma_coef = summary(mod_sigma)
	sigma_varcomp = varcomp(mod_sigma)
	
	write_rds(list(estimates = data.frame(alpha = alpha_mat, A = A_mat, snow = snow_mat,
																				tau = tau_mat, tao = tao_mat, sigma = sigma_mat),
								 alpha = alpha_coef$fixed.effects.estimates, alpha_v = alpha_varcomp,
								 A = A_coef$fixed.effects.estimates, A_v = A_varcomp,
								 snow = snow_coef$fixed.effects.estimates, snow_v = snow_varcomp,
								 tau = tau_coef$fixed.effects.estimates, tau_v = tau_varcomp,
								 tao = tao_coef$fixed.effects.estimates, tao_v = tao_varcomp,
								 sigma = sigma_coef$fixed.effects.estimates, sigma_v = sigma_varcomp),
						path = paste("./Water_Temperature/Data/iterations/",as.character(i),".rds",sep = ""))
	
	list(estimates = data.frame(alpha = alpha_mat, A = A_mat, snow = snow_mat,
															tau = tau_mat, tao = tao_mat, sigma = sigma_mat),
			 alpha = alpha_coef$fixed.effects.estimates, alpha_v = alpha_varcomp,
			 A = A_coef$fixed.effects.estimates, A_v = A_varcomp,
			 snow = snow_coef$fixed.effects.estimates, snow_v = snow_varcomp,
			 tau = tau_coef$fixed.effects.estimates, tau_v = tau_varcomp,
			 tao = tao_coef$fixed.effects.estimates, tao_v = tao_varcomp,
			 sigma = sigma_coef$fixed.effects.estimates, sigma_v = sigma_varcomp)
}
Sys.time()-start_time
write_rds(rtrn, "./Water_Temperature/Data/network_250.rds")
stopCluster(cl)
#-------------------------> END Sequence <-------------------------#


#------------------------> Model Checking <------------------------#

alpha_cv = CrossValidationSSN(mod_alpha)
A_cv = CrossValidationSSN(mod_A)
snow_cv = CrossValidationSSN(mod_snow)
tau_cv = CrossValidationSSN(mod_tau)
tao_cv = CrossValidationSSN(mod_tao)
sigma_cv = CrossValidationSSN(mod_sigma)

cv = data.frame(obs = c(mod_alpha$sampinfo$z,
									 mod_A$sampinfo$z,
									 mod_snow$sampinfo$z,
									 mod_tau$sampinfo$z,
									 mod_tao$sampinfo$z,
									 mod_sigma$sampinfo$z),
					 cv_pred = c(alpha_cv[,"cv.pred"],
					 						A_cv[,"cv.pred"],
					 						snow_cv[,"cv.pred"],
					 						tau_cv[,"cv.pred"],
					 						tao_cv[,"cv.pred"],
					 						sigma_cv[,"cv.pred"]),
					 model = c(rep("alpha",length(alpha_cv[,"cv.pred"])),
					 						rep("A",length(A_cv[,"cv.pred"])),
					 						rep("snow",length(snow_cv[,"cv.pred"])),
					 						rep("tau",length(tau_cv[,"cv.pred"])),
					 						rep("tao",length(tao_cv[,"cv.pred"])),
					 						rep("sigma",length(sigma_cv[,"cv.pred"]))))

cv$labels = factor(x = cv$model, labels = c(
	expression(italic("A")),
	expression(alpha),
	expression(sigma),
	expression(phi),
	expression(tau[2]),
	expression(tau[1])
	))

p=ggplot(cv, aes(obs,cv_pred)) + geom_point(size = .5) +
	geom_abline(intercept = 0, slope = 1) +
	facet_wrap(~labels, scales = "free", labeller = label_parsed) + 
	ggthemes::theme_tufte() +
	labs(x = "Coefficient Observations", y = "Coefficient Predictions")
ggsave(filename = "cross_val.png", plot = p, device = "png", 
			 path = "./Water_Temperature/images/", width = 6.5, height = 4.5, units = "in", dpi = 500)

#----------------------> END Model Checking <----------------------#



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
start = Sys.time()
#Sample an early, mid and late summer migration date and calculate temperatures.
probs = plyr::dlply(.data = samp_dates(single=T), .variables = "migrate", .fun = function(x){
	
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
	preds = getSSNdata.frame(network, Name = "preds"); names(preds)[4] = "year"
	preds = mutate(preds, rowNo = row_number())
	
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
		celcius = (df$alpha + df$A*cos(2*pi*df$d/df$n+pi*df$tau)) + #Annual
			cos(2*pi*df$d/(df$n/2)+pi*df$tao)   						          #Seasonal
		p19 = pnorm(19,celcius,df$sigma,lower.tail = F)							#Prob. over 19ºC
		p22 = pnorm(22,celcius,df$sigma,lower.tail = F)							#Prob. over 22ºC
		data.frame(celcius,p19,p22)
	}
	#Determine initial location of each recuring prediction sequence by mcmc chain.
	prob_row = seq(1,length(fits)*nrow(preds),by=nrow(preds))-1
	
	#----->Analysis<-----#
	plyr::dlply(.data = x, .variables = "start", .fun = function(y){
		tmp = plyr::ldply(fits, function(z){
			n = if_else(leap_year(ymd(paste(as.character(preds$year),"01","01",sep="-"))),366,365)
			est = mutate(z$estimates, d = swim_date(y$start,preds$upDist),n = n)
			temperature(est)
		})
		preds %>% group_by(locID) %>% do({
			if(.$netID==1) net = net1
			if(.$netID==2) net = net2
			binary = net[net$rid==.$rid,"binaryID"]
			#Determine downstream prediction points.
			dwn_rid = dwnStrm(net, binary)
			#Probability of >#ºC experience by ...
			#Site or ...
			pt = tmp[c(prob_row+.$rowNo),]
			#Migration route
			pts = filter(preds, year == .$year, upDist<=.$upDist, rid %in% dwn_rid) %>% 
				group_by(rowNo) %>% do({
					tmp[c(prob_row+.$rowNo),]
				})
			#Sum all mean temperature estimates above threshold,
			#Determine the mean probability of being above threshold during migration,
			#Calculate the median probability of being above threshold at the migration destination.
			data.frame(
				sum19 = sum(pts$celcius>19), sum22 = sum(pts$celcius>22),
				G19 = sum(pts$celcius>19)/nrow(pts), G22 = sum(pts$celcius>22)/nrow(pts),
				ptG19 = median(pt$p19), ptG22 = median(pt$p22)
			)
		})
	}, .parallel = T, .paropts = list(.packages = c('tidyverse','lubridate'),
																		.export = c("fits","preds","net1","net2",
																								"swim_date","temperature",
																								"samp_dates","dwnStrm")))
})
Sys.time()-start 
write_rds(probs, "./Water_Temperature/Data/probs_flatD.rds")
stopCluster(cl)

#The approximate distance from the Albion test fishery (Parken et al. 2008) to the Thompson as Kamloops is 447km. Adding this initial distance to the distance upstream from Kamloops give the total travel distance. We can then divide by an approximate travel rate of 36km/day (Salinger and Anderson 2006).

probs = read_rds("./Water_Temperature/Data/probs_flatD.rds")
#Summarize by migration group and place in shapefile for mapping.
lim_summary = function(probs, cols, prefix){
	plyr::llply(probs, .fun=function(x){
		tmp = plyr::ldply(x, .fun=function(y){
			data.frame(locID = y$locID, cum = y[[cols[1]]], probs = y[[cols[2]]],
								 probs_s = y[[cols[3]]], year = preds$year)
		})
		#Get the mean by year across days.
		mean_c = tmp %>% group_by(year,locID) %>% 
			summarise(mean_p = median(probs), mean_c = median(cum), mean_s = median(probs_s))
		#Get the first day with the greatest average value.
		max_p = tmp %>% group_by(year,start) %>% summarise(mean_p = median(probs)) %>% 
			filter(mean_p == max(mean_p)) %>% filter(start == min(start)) %>% rename(start_p = start) %>% 
			inner_join(.,tmp,by = c("year","start_p"="start")) %>% select(-mean_p, -cum, -probs_s)
		max_c = tmp %>% group_by(year,start) %>% summarise(mean_c = mean(cum)) %>% group_by(year) %>%  
			filter(mean_c == max(mean_c)) %>% filter(start == min(start)) %>% rename(start_c = start) %>% 
			inner_join(.,tmp,by = c("year","start_c"="start")) %>% select(-mean_c, -probs, -probs_s)
		max = left_join(max_p, max_c)
		#return summarized data by migration.
		final = max %>% left_join(., mean_c, by = c("year","locID")) %>% select(3,1,7,8,9,2,4,5,6)
		names(final)[c(3:9)] = paste(prefix,names(final)[c(3:9)],sep="_")
		return(final)
	})
}
#Summarize by threshold
G19 = lim_summary(probs, c(2,4,6), "G19")
G22 = lim_summary(probs, c(3,5,7), "G22")

#Save shape files
shp_save = function(pred_df, period, yr){
	#Bind together the thresholds and Filter out by year.
	p_riod = pred_df %>% select(-c(1,2,5:9)) %>% 
		left_join(.,G19[[period]]) %>% left_join(.,G22[[period]]) %>% filter(year == yr)
	#Make spatial object
	p_riod = SpatialLinesDataFrame(sl = as.SpatialLines(network), data = p_riod, match.ID = "UID")
	#Write shapefile.
	lyr = paste(period, as.character(yr), sep = "_")
	writeOGR(p_riod, dsn = "~/sfu/River_Network_Temperature_Output/risk_probs/", layer = lyr,
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
