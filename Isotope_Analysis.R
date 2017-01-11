setwd("~/Documents/Simon_Fraser_University/PhD_Research/Projects/River-Network-Temperature")
library(ggplot2);library(dplyr)
dat = read.csv(file = "Isotope Sources of Variation.csv", header = T, stringsAsFactors = F) %>% select(1:9)
dat = dat[-c(1,2),] %>% tbl_df()
names(dat)[1:7] = c("Source", "D_VarL", "D_Mean", "D_VarH", "O_VarL", "O_Mean", "O_VarH")

dat_p =dat %>% na.omit(D_Mean) %>% arrange(desc(D_Mean)) 
#Levels = as.character(dat_p$Source); dat_p$Source = factor(dat_p$Source, levels = Levels)
Levels = dat_p$Location; dat_p$Location = factor(dat_p$Location, levels = Levels)
ggplot(dat_p, aes(Location, D_Mean, colour = Water_State))+
	geom_point() +
	geom_linerange(aes(ymax = D_VarH, ymin = D_VarL)) +
	coord_flip() +
	theme_bw() +
	geom_hline(yintercept = 0) +
	labs(y = "Deuterium Mean and ~Variation", x = "~Location")
ggsave(filename = "P1_Deuterium.pdf", width = 21, height = 14, units = "in")

dat_p =dat %>% na.omit(O_Mean) %>% arrange(desc(O_Mean))
#Levels = as.character(dat_p$Source); dat_p$Source = factor(dat_p$Source, levels = Levels)
Levels = dat_p$Location; dat_p$Location = factor(dat_p$Location, levels = Levels)
ggplot(dat_p, aes(Location, O_Mean, colour = Water_State))+
	geom_point() +
	geom_linerange(aes(ymax = O_VarH, ymin = O_VarL)) +
	coord_flip() +
	theme_bw() +
	geom_hline(yintercept = 0) +
	labs(y = "O18 Mean and ~Variation", x = "~Location")
ggsave(filename = "P2_Oxygen18.pdf", width = 21, height = 14, units = "in")


ggplot(dat, aes(O_Mean, D_Mean, colour = Water_State, label = Location)) +
	geom_point() +
	geom_errorbar(aes(ymax = D_VarH, ymin = D_VarL)) +
	geom_errorbarh(aes(xmax = O_VarH, xmin = O_VarL)) +
	geom_text(check_overlap = F, nudge_y = 50) + #, nudge_x = -4) +
	theme_bw() +
	geom_hline(yintercept = 0) +
	geom_vline(xintercept = 0) +
	labs(y = "Deuterium Mean and ~Variation", x = "O18 Mean and ~Variation")
ggsave(filename = "P3_D_vs_18O.pdf", width = 21, height = 14, units = "in")





path = "~/Documents/Simon_Fraser_University/PhD_Research/Projects/Data/Original Data/IAEA WISER Isotope Data"
Precip = read.csv(paste(path,"GNIP.csv", sep = "/"), header = T)
Stream = read.csv(paste(path,"GNIR.csv", sep = "/"), header = T)

Precip = Precip %>% select(Name.of.site, Country, Date, O18, H2, Precipitation, Air.Temperature) %>% tbl_df
date = as.Date(paste(as.character(Precip$Date),"-01", sep = ""))
month = as.factor(format(date, "%b"))
year = format(date, "%Y")
Precip = Precip %>% mutate(month,year)

snow.rain = apply(Precip, 1, function(x){
	if(x["month"] == "Jan"| x["month"] == "Feb" | x["month"] == "Mar") "Snow"
	else "Rain"
})
Precip = Precip %>% mutate(snow.rain)

Stream = Stream %>% select(Name.of.site, Country, Date, O18, H2) %>% tbl_df()

ggplot(Precip, aes(O18,H2,color = snow.rain))+
	stat_smooth(method = "lm") +
	geom_point() +
	theme_bw() +
	facet_wrap(~Name.of.site, scales = "free_x")


dat = Precip %>% dplyr::filter(month=="Jan" | month =="Feb" | month =="Mar" | month =="Jun" | month =="Jul" | month =="Aug" & Name.of.site=="OLIVER")
ggplot(dat, aes(O18,H2, group = year))+
	stat_smooth(method = "lm",se = F, color = "grey") +
	geom_point(aes(color = snow.rain)) +
	theme_bw() +
	facet_wrap(~year)


###################################
# When to collect in the Thompson #
###################################

#http://bcrfc.env.gov.bc.ca/data/asp/realtime/
#Mount Cook:1E02P

data1 = read.csv(file = "http://bcrfc.env.gov.bc.ca/data/asp/realtime/data/1E02P.csv", header = T)
#Create Site ID and Date without time variables.
Date = unlist(lapply(strsplit(as.character(data1$Date),split = " "),function(x) x[1]))
#Rename columns, summarize by date and organize dataframe.
data1 = data.frame(data1,Date = as.Date(Date))
data1 = data1 %>% group_by(Date) %>% summarize(Temp_Max = max(Temp..C.,na.rm = T), Temp_Min = min(Temp..C.,na.rm = T), Precipitation = max(RGauge..mm.,na.rm = T), Accumulated_Precip = max(RGauge..mm.+SWE..mm.,na.rm = T), Snow_Water_Equivalent = max(SWE..mm.,na.rm = T), Snow_Depth = max(SD..cm.,na.rm = T))
data1 = na.omit(data1)
#Tidy the datafame.
temp = reshape2::melt(data1[,c(2:7)],na.rm = T)
Pillow_ID = rep("1E02P",nrow(temp))
data2 = data.frame(Pillow_ID=Pillow_ID, Date = rep(data1$Date,6), temp)

ggplot(data2, aes(Date,value)) +
	geom_point() +
	theme_bw() +
	facet_wrap(~variable, scales = "free_y")

year = as.numeric(format(data2$Date, "%Y")); month = as.factor(format(data2$Date, "%b"))
data2 = bind_cols(data2[,c(1,2)],data.frame(year,month),data2[,c(3,4)])

data = read.csv(file = "~/Documents/Simon_Fraser_University/PhD_Research/Projects/Data/Original Data/GIS Data/BC Snow Survey/Locations Automated/Daily_Snow_Depth_Auto_71-15.csv", header = T)
data$Date = as.Date(data$Date)
year = as.numeric(format(data$Date, "%Y")); month = as.factor(format(data$Date, "%b"))
data = bind_cols(data[,c(1,2)],data.frame(year,month),data[,c(3,4)])
data = data %>% dplyr::filter(year == 2010 | year == 2011) %>% dplyr::filter(Pillow_ID == "1E02P")
data = data[-which(data$variable=="Snow_Depth" & data$value > 400),]
data %>% group_by(variable,year) %>% dplyr::filter(value == max(value,na.rm=T)) %>%  dplyr::filter(Date == max(Date))

data3 = bind_rows(data,data2)
ggplot(data3, aes(Date,value, colour = as.factor(year))) +
	geom_point() +
	theme_bw() +
	theme(legend.position = "none") +
	facet_wrap(~variable, scales = "free_y") +
	geom_vline(xintercept = as.numeric(as.Date(c("2011-06-24", "2010-06-24", "2016-06-24")))) +
	geom_vline(xintercept = as.numeric(as.Date(c("2011-04-16", "2010-04-08", "2016-04-16", "2016-02-03"))), colour = "blue") +
	geom_vline(xintercept = as.numeric(as.Date(c("2011-05-20", "2010-05-13", "2016-05-20", "2016-05-13"))), colour = "red")
ggsave(filename = "Thompson_Precipitation.pdf",width = 11, height = 4.5)




