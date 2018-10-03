library(tidyverse); library(stringr)
setwd("/Users/kylechezik/sfu/River_Network_Temperature_Output/predictors/csv_stats")

files = dir(); full = list()
for(i in files){
	dat = read_csv(i, col_names = F); names(dat) = c("GRIDCODE","COUNT")
	dat = mutate(dat, YEAR = as.numeric(str_extract(i, "\\d+")) + 1900)
	full[i] = list(dat)
}

df = data.frame(do.call("rbind",full), row.names = NULL)
df = spread(df,key = YEAR, value = COUNT)

write_csv(df, path = "../wildfire_year_2_zonal.csv")
