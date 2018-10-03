#This script filters the logging_year.tif file by year, which are then summarized by watershed within the watersheds.shp file using the 'zonalstats' tool within the 'rasterstats' module in python. Turning the watersheds.shp file into a .geojson file is necessary for it to be read by the 'zonalstats' tool. This entire process is run in parallel and outputs .csv files using a 'json2csv' tool. Raster .tif files will be cleaned up as they are created and summarized but will require at least 20 GB of space while the program is running. This program uses 'parallel GUI', 'Fiona', 'Rasterio', 'gdal' and 'json2csv'.

seq 85 110 > years.txt

cat years.txt | parallel -I @@ \
	'gdal_calc.py -A /Users/kylechezik/sfu/River_Network_Temperature_Output/predictors/wildfires_year_1.tiff --A_band=1 --outfile=/Users/kylechezik/sfu/River_Network_Temperature_Output/predictors/logging_year_1_@@.tiff --calc="where(A==@@, @@, -32768)" --NoDataValue=-32768 --type="Int16"; fio cat /Users/kylechezik/sfu/River_Network_Temperature_Output/shp_costRCAs/costrcas.shp | rio zonalstats -r /Users/kylechezik/sfu/River_Network_Temperature_Output/predictors/logging_year_1_@@.tiff --sequence --prefix "year_" --stats "count" | json2csv -k properties.GRIDCODE,properties.year_count >> /Users/kylechezik/sfu/River_Network_Temperature_Output/predictors/csv_stats/stats_@@.csv; rm /Users/kylechezik/sfu/River_Network_Temperature_Output/predictors/logging_year_1_@@.tiff'

rm years.txt
