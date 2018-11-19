library(tidyverse); library(SSN); library(rgdal); library(boot); library(lubridate);
library(foreach); library(doParallel) #Activate Libraries

fits = read_rds("./Water_Temperature/Data/network_250.rds")

#Alpha
alpha = plyr::ldply(fits, function(x){
	x$alpha[,c("FactorLevel","Estimate")]
})

p1 = alpha %>% filter(FactorLevel %in% c("scl_air_mn","scl_air_A","scl_ele")) %>% 
	ggplot(aes(Estimate, group=FactorLevel, fill=FactorLevel)) + 
	geom_density(color = "white", adjust = 1.5) +
	scale_fill_manual(values = c("#F26419","#F6AE2D","#45062E"),
										labels = c("Air Amplitude","Air Mean","Elevation")) +
	ggthemes::theme_tufte() +
	labs(x = "", fill = "", 
			 title = expression(paste("Mean Water Temperature (",alpha,")"))) + 
	theme(legend.direction = "horizontal", legend.position = "top",
				plot.title = element_text(hjust = 0.5, size = 10),
				legend.text.align = 0,
				axis.text.y = element_blank(),
				axis.ticks.y = element_blank(),
				axis.title.y = element_blank(),
				axis.text.x = element_text(size = 6),
				axis.title.x = element_text(size = 8),
				legend.text = element_text(size = 6),
				legend.key.size = unit(2,'mm'))

p2 = alpha %>% filter(FactorLevel %in% c("lgt_firecum","lgt_lake",
																	 "lgt_lake:scl_ele")) %>% 
	ggplot(aes(Estimate, group=FactorLevel, fill=FactorLevel)) + 
	geom_density(color = "white", adjust = 1.5) +
	scale_fill_manual(values = c("#F6AE2D","#33658A", "#45062E"),
										labels = c("WildFire","Lakes","Lakes:Elevation")) +
	ggthemes::theme_tufte() +
	labs(x = "Parameter Estimate", fill = "") + 
	theme(legend.direction = "horizontal", legend.position = "top",
				plot.title = element_blank(),
				legend.text.align = 0,
				axis.text.y = element_blank(),
				axis.ticks.y = element_blank(),
				axis.title.y = element_blank(),
				axis.text.x = element_text(size = 6),
				axis.title.x = element_text(size = 8),
				legend.text = element_text(size = 6),
				legend.key.size = unit(2,'mm'))

#pf = gridExtra::grid.arrange(p1, p2, layout_matrix = matrix(data = c(1,2), nrow = 2, ncol = 1))
#ggsave(filename = "SSN_alpha.png", path = "~/sfuvault/River-Network-Temperature/Water_Temperature/images/", plot = pf, width = 6.5, height = 3.5, units = "in", device = "png")


#A
A = plyr::ldply(fits, function(x){
	x$A[,c("FactorLevel","Estimate")]
})

p3 = A %>% filter(FactorLevel %in% c("scl_air_mn","scl_air_A","scl_ca")) %>% 
	ggplot(aes(Estimate, group=FactorLevel, fill=FactorLevel)) + 
	geom_density(color = "white", adjust = 1.5) +
	scale_fill_manual(values = c("#F26419","#F6AE2D","#9ECE9A"),
										labels = c("Air Amplitude","Air Mean",
															 expression(paste("Catchment Area (H"[2],"O)")))) +
	ggthemes::theme_tufte() +
	xlim(0,0.75) +
	labs(x = "", fill = "", 
			 title = expression(paste("Water Temperature Amplitude (",italic("A"),")"))) + 
	theme(legend.direction = "horizontal", legend.position = "top",
				plot.title = element_text(hjust = 0.5, size = 10),
				legend.text.align = 0,
				axis.text.y = element_blank(),
				axis.ticks.y = element_blank(),
				axis.title.y = element_blank(),
				axis.text.x = element_text(size = 6),
				axis.title.x = element_text(size = 8),
				legend.text = element_text(size = 6),
				legend.key.size = unit(2,'mm'))

p4 = A %>% filter(FactorLevel %in% c("lgt_glacier","lgt_lake")) %>% 
	ggplot(aes(Estimate, group=FactorLevel, fill=FactorLevel)) + 
	geom_density(color = "white", adjust = 1.5) +
	scale_fill_manual(values = c("#86BBD8","#33658A"),
										labels = c("Glacier","Lakes")) +
	ggthemes::theme_tufte() +
	xlim(-0.25,0.35) +
	labs(x = "Parameter Estimate", fill = "") + 
	theme(legend.direction = "horizontal", legend.position = "top",
				plot.title = element_blank(),
				legend.text.align = 0,
				axis.text.y = element_blank(),
				axis.ticks.y = element_blank(),
				axis.title.y = element_blank(),
				axis.text.x = element_text(size = 6),
				axis.title.x = element_text(size = 8),
				legend.text = element_text(size = 6),
				legend.key.size = unit(2,'mm'))

#Snow
snow = plyr::ldply(fits, function(x){
	x$snow[,c("FactorLevel","Estimate")]
})

p5 = snow %>% filter(FactorLevel %in% c("scl_air_mn","scl_ele")) %>% 
	ggplot(aes(Estimate, group=FactorLevel, fill=FactorLevel)) + 
	geom_density(color = "white", adjust = 1.5) +
	scale_fill_manual(values = c("#F6AE2D","#45062E"),
										labels = c("Air Mean","Elevation")) +
	xlim(-1.25,-0.5) +
	ggthemes::theme_tufte() +
	labs(x = "", fill = "", 
			 title = expression(paste("Seasonal Hysteresis (",phi,")"))) +
	theme(legend.direction = "horizontal", legend.position = "top",
				plot.title = element_text(hjust = 0.5, size = 10),
				legend.text.align = 0,
				axis.text.y = element_blank(),
				axis.ticks.y = element_blank(),
				axis.title.y = element_blank(),
				axis.text.x = element_text(size = 6),
				axis.title.x = element_text(size = 8),
				legend.text = element_text(size = 6),
				legend.key.size = unit(2,'mm'))
	

p6 = snow %>% filter(FactorLevel %in% c("lgt_glacier")) %>% 
	ggplot(aes(Estimate, group=FactorLevel, fill=FactorLevel)) + 
	geom_density(color = "white", adjust = 1.5) +
	scale_fill_manual(values = c("#86BBD8"),
										labels = c("Glacier")) +
	ggthemes::theme_tufte() +
	labs(x = "Parameter Estimate", fill = "") + 
	xlim(-0.35,-0.25)+
	theme(legend.direction = "horizontal", legend.position = "top",
				plot.title = element_blank(),
				legend.text.align = 0,
				axis.text.y = element_blank(),
				axis.ticks.y = element_blank(),
				axis.title.y = element_blank(),
				axis.text.x = element_text(size = 6),
				axis.title.x = element_text(size = 8),
				legend.text = element_text(size = 6),
				legend.key.size = unit(2,'mm'))

#pf = gridExtra::grid.arrange(p1, p2, p3, p4, layout_matrix = matrix(data = c(1,2,3,4),
#																																		nrow = 2, ncol = 2))

pf = gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6,
														 layout_matrix = matrix(data = c(1,2,3,4,5,6),
																																		nrow = 2, ncol = 3))
ggsave(filename = "SSN.png", path = "~/sfuvault/River-Network-Temperature/Water_Temperature/images/", plot = pf, width = 6.5, height = 3.5, units = "in", device = "png")

#Remaining variable estimates
tau = plyr::ldply(fits, function(x){
	x$tau[,c("FactorLevel","Estimate")]
})

tau %>% filter(FactorLevel == "scl_ele") %>% .$Estimate %>% summary(.)

tao = plyr::ldply(fits, function(x){
	x$tao[,c("FactorLevel","Estimate")]
})

tao %>% filter(FactorLevel == "scl_ele") %>% .$Estimate %>% summary(.)

sigma = plyr::ldply(fits, function(x){
	x$sigma[,c("FactorLevel","Estimate")]
})

sigma %>% filter(FactorLevel == "scl_ele") %>% .$Estimate %>% summary(.)

#Variance Components
alpha_vc = plyr::ldply(fits, function(x){
	x$alpha_v
}) %>% group_by(VarComp) %>% summarise(prop_mn = mean(Proportion))

A_vc = plyr::ldply(fits, function(x){
	x$A_v
}) %>% group_by(VarComp) %>% summarise(prop_mn = mean(Proportion))

snow_vc = plyr::ldply(fits, function(x){
	x$snow_v
}) %>% group_by(VarComp) %>% summarise(prop_mn = mean(Proportion))

tau_vc = plyr::ldply(fits, function(x){
	x$tau_v
}) %>% group_by(VarComp) %>% summarise(prop_mn = mean(Proportion))

tao_vc = plyr::ldply(fits, function(x){
	x$tao_v
}) %>% group_by(VarComp) %>% summarise(prop_mn = mean(Proportion))

sigma_vc = plyr::ldply(fits, function(x){
	x$sigma_v
}) %>% group_by(VarComp) %>% summarise(prop_mn = mean(Proportion))


#3D shifting probabilities over the season
library(plot3D)
probs = read_rds("./Water_Temperature/Data/probs_flatD.rds")

png(filename = "./Water_Temperature/images/Global_Migration_Risk.png",
		width = 3.5, height = 3.5, units = "in", res = 500)

par(mar = c(1,1,.5,0))
d = density(probs$early[[1]]$G19)
scatter3D(x = rep(1,length(d$x)), y = d$x, z = d$y, 
					type ="l", xlim = c(1,92), zlim = c(0,4),
					colvar = NULL, col = "#F6AE2D",
					xlab = "Time (days)", ylab = "Risk Probability", zlab = "",
					phi = 20, theta = 240, bty = "u")

for(i in 2:length(probs$early)){
	d = density(probs$early[[i]]$G19)
	lines3D(x = rep(i,length(d$x)), y = d$x, z = d$y, add = T, col = "#F6AE2D", colvar = NULL)
}
for(i in 1:length(probs$mid)){
	d = density(probs$mid[[i]]$G19)
	ii = i+30
	lines3D(x = rep(ii,length(d$x)), y = d$x, z = d$y, add = T, col = "#F26419", colvar = NULL)
}
for(i in 1:length(probs$late)){
	d = density(probs$late[[i]]$G19)
	ii = i+61
	lines3D(x = rep(ii,length(d$x)), y = d$x, z = d$y, add = T, col = "#45062E90", colvar = NULL)
}
dev.off()



#Chinook Pop. Results

chinook = readOGR(dsn = "~/sfu/River_Network_Temperature_Output/chinook_populations/", layer = "pops")
c_dat = chinook@data
c_dat$wtsh = c(rep("N",8),rep("S",8))


p_eds = preds %>% select(UID,year,rid,upDist)
#Combine the prediction and risk data for 19ºC and 22ºC.
#Then combine chinook location data and preds/results.
c_dat = c_dat %>% group_by(rid) %>% do({
	if(.$migr=="Early"){
		rows = which(p_eds$rid==.$rid)
		p19 = G19$early$G19_mean_p[rows]
		c19 = G19$early$G19_mean_c[rows]
		s19 = G19$early$G19_mean_s[rows]
		p22 = G22$early$G22_mean_p[rows]
		c22 = G22$early$G22_mean_c[rows]
		s22 = G22$early$G22_mean_s[rows]
	}
	if(.$migr=="Mid"){
		rows = which(p_eds$rid==.$rid)
		p19 = G19$mid$G19_mean_p[rows]
		c19 = G19$mid$G19_mean_c[rows]
		s19 = G19$mid$G19_mean_s[rows]
		p22 = G22$mid$G22_mean_p[rows]
		c22 = G22$mid$G22_mean_c[rows]
		s22 = G22$mid$G22_mean_s[rows]
	}
	if(.$migr=="Late"){
		rows = which(p_eds$rid==.$rid)
		p19 = G19$late$G19_mean_p[rows]
		c19 = G19$late$G19_mean_c[rows]
		s19 = G19$late$G19_mean_s[rows]
		p22 = G22$late$G22_mean_p[rows]
		c22 = G22$late$G22_mean_c[rows]
		s22 = G22$late$G22_mean_s[rows]
	}
	data.frame(rid = .$rid, year = c(2014:2017),
						 p19 = p19, c19 = c19, s19 = s19,
						 p22 = p22, c22 = c22, s22 = s22)
}) %>% left_join(.,p_eds) %>% left_join(c_dat, .)

#Plot
dat =  c_dat %>% arrange(desc(wtsh),desc(upDist),desc(c19))
x = dat %>% select(Pop) %>% distinct() %>% .$Pop %>% as.numeric()
dat$Pop = factor(dat$Pop,levels = levels(dat$Pop)[c(x)])
dat$migr = factor(dat$migr, levels = c("Early","Mid","Late"))

p1 = ggplot(dat, aes(Pop,c19,color=as.factor(year),shape = migr)) +
	geom_point() +
	scale_color_manual(values = c("#F26419","#F6AE2D","#86BBD8","#33658A"))+
	ggthemes::theme_tufte() +
	labs(y = expression(paste(italic("n")," >19"*degree*"C")),
			 title = "Cumulative Migratory Thermal Exposure",
			 x = "", color = "Year", shape = "Watershed") + 
	theme(legend.direction = "horizontal", legend.position = c(.8,1.15),
				legend.text.align = 0, 
				plot.title = element_text(hjust = 0.5, size = 10,
																	margin = margin(t = 10, b = 30)),
				axis.text = element_text(size = 7),
				axis.title = element_text(size = 9),
				legend.text = element_text(size = 7),
				legend.title = element_text(size = 8),
				legend.key.size = unit(2,'mm'),
				axis.text.x = element_blank(),
				plot.margin=unit(c(.1,.1,.1,.3),"cm"))+
	guides(shape = F,
				 color = guide_legend(title.position="top", title.hjust = 0.5))

p2 = ggplot(dat, aes(Pop,c22,color=as.factor(year), shape = migr)) +
	geom_point() +
	scale_color_manual(values = c("#F26419","#F6AE2D","#86BBD8","#33658A"))+
	ggthemes::theme_tufte() +
	labs(y = expression(paste(italic("n")," >22"*degree*"C")),
			 x = "Chinook Populations") + 
	theme(axis.text = element_text(size = 7),
				axis.title = element_text(size = 9),
				axis.text.x = element_text(angle = 35, vjust = .6),
				plot.margin=unit(c(.1,.1,.3,.3),"cm"))+
	guides(shape = F, color=F)

dat =  c_dat %>% arrange(desc(wtsh),desc(upDist),desc(s19))
x = dat %>% select(Pop) %>% distinct() %>% .$Pop %>% as.numeric()
dat$Pop = factor(dat$Pop,levels = levels(dat$Pop)[c(x)])
dat$migr = factor(dat$migr, levels = c("Early","Mid","Late"))

p3 = ggplot(dat, aes(Pop,s19,color=as.factor(year),shape = migr)) + 
	geom_point() +
	scale_color_manual(values = c("#F26419","#F6AE2D","#86BBD8","#33658A"))+
	ggthemes::theme_tufte() +
	labs(y = expression(paste(italic("P(y")[italic("s,yr")],italic(")")," >19"*degree*"C")), 
			 title = "Spawn Site | Thermal Exposure Probability",
			 x = "", shape = "Migration") + 
	theme(legend.direction = "horizontal", legend.position = c(.2,1.15),
				legend.text.align = 0, 
				plot.title = element_text(hjust = 0.5, size = 10,
																	margin = margin(t = 10, b = 30)),
				axis.text = element_text(size = 7),
				axis.title = element_text(size = 9),
				legend.text = element_text(size = 7),
				legend.title = element_text(size = 8),
				legend.key.size = unit(2,'mm'),
				axis.text.x = element_blank(),
				plot.margin=unit(c(.1,.1,.1,.3),"cm"))+
	guides(shape = guide_legend(title.position="top", title.hjust = 0.5),
				 color = F)

p4 = ggplot(dat, aes(Pop,s22,color=as.factor(year), shape = migr)) + 
	geom_point() +
	scale_color_manual(values = c("#F26419","#F6AE2D","#86BBD8","#33658A"))+
	ggthemes::theme_tufte() +
	labs(y = expression(paste(italic("P(y")[italic("s,yr")],italic(")")," >22"*degree*"C")),
			 x = "Chinook Populations") + 
	theme(axis.text = element_text(size = 7),
				axis.title = element_text(size = 9),
				axis.text.x = element_text(angle = 35, vjust = .6),
				plot.margin=unit(c(.1,.1,.3,.3),"cm"))+
	guides(shape = F, color=F)


pf = gridExtra::grid.arrange(p1, p2, p3, p4, 
														 layout_matrix = matrix(data = c(1,2,3,4), nrow = 2, ncol = 2))
ggsave(filename = "./Water_Temperature/images/chinook.png",
			 plot = pf, device = "png", width = 7.5, height = 5,
			 units = "in", dpi = 500)

	#Reorder factors
dat$migr = factor(dat$migr, levels = c("Early","Mid","Late"))
#Plots
p1 = ggplot(dat, aes(migr, p19, color = as.factor(year))) + geom_jitter(width = .1, height = 0) +
	scale_color_manual(values = c("#F26419","#F6AE2D","#86BBD8","#33658A")) +
	labs(y = expression(paste(italic("P(y")[italic("s,yr")],italic(")")," >19"*degree*"C"))) +
	ggthemes::theme_tufte()+
	theme(legend.direction = "horizontal", legend.position = "top",
				legend.title = element_blank(),
				axis.text = element_text(size = 7),
				axis.text.x = element_blank(),
				axis.title.x = element_blank()) +
	guides(color = guide_legend(title.position="top")) +
	facet_wrap(~wtsh)

p2 = ggplot(dat, aes(migr, p22, color = as.factor(year))) + geom_jitter(width = .1, height = 0) +
	scale_color_manual(values = c("#F26419","#F6AE2D","#86BBD8","#33658A")) +
	labs(y = expression(paste(italic("P(y")[italic("s,yr")],italic(")")," >22"*degree*"C")), 
			 x = "Migration Period") +
	ggthemes::theme_tufte()+
	theme(legend.position = "none",
				legend.title = element_blank(),
				axis.text = element_text(size = 7),
				axis.title = element_text(size = 9),
				strip.text = element_blank()) +
	facet_wrap(~wtsh)

pf = gridExtra::grid.arrange(p1, p2,
														 layout_matrix = matrix(data = c(1,2), nrow = 2, ncol = 1))
ggsave(filename = "./Water_Temperature/images/chinook_migr.png",
			 plot = pf, device = "png", width = 3.75, height = 5,
			 units = "in", dpi = 500)
