all: clip_region rasterize_streams

# Clip the BC DEM to the region of concern.
clip_region:
	gdalwarp -overwrite -s_srs '+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs' -t_srs EPSG:3005 -r cubicspline -tr 25 -25 -srcnodata -999 -dstnodata -32768 -cutline ./shp_clip_region/Clip.shp -crop_to_cutline /Users/kylechezik/Documents/BC_DEM25/bcdem.tif Thompson_UL.tif

# Rasterize the streams we would like burned into the DEM.
rasterize_streams:
	gdal_rasterize -l brn_strms -at -burn 5 -init 0 -a_srs EPSG:3005 -te 1291224.0525986417196691 578027.7679750561947003 1598624.0525986417196691 898402.7679750561947003 -tr 25 -25 ./shp_burn_streams/brn_strms.shp brn_strms_rast.tif 

# Burn streams into DEM.
burn_streams:
	gdal_calc.py -A Thompson_UL.tif -B brn_strms_rast.tif --outfile=burn.tif --calc="A-B" 
# Cleanup files.
	rm Thompson_UL.tif

# Preprocess DEM
sink_route:
	saga_cmd ta_preprocessor 1 -ELEVATION "burn.tif" -SINKROUTE "sink_route"

# Rasterize pour point.
rasterize_pour_points:
	gdal_rasterize -l pour_2pts -a id -a_srs EPSG:3005 -a_nodata -32768 -te 1291224.0525986417196691 578027.7679750561947003 1598624.0525986417196691 898402.7679750561947003 -tr 25 -25 ./shp_pour_points/pour_2pts.shp pour_2pts.tif 

# Delineate Thompson Watershed
watershed:
	saga_cmd ta_channels 1 -ELEVATION "burn.tif" -CHANNELS "pour_2pts.tif" -SINKROUTE "sink_route.sdat" -BASINS "watershed"
	rm pour_2pts.tif sink_route*

vectorize_watershed:
	saga_cmd shapes_grid 6 -GRID "watershed.sdat" -POLYGONS "./watersheds/watershed.shp"
	rm watershed.*
#Used the Delete holes tool in QGIS to clean up small pockets in the watershed and the smooth geometry tool to remove self-intersection along the outer edge.

# Clip DEM to Thompson watershed.
clip_DEM_watershed:
	gdalwarp -overwrite -s_srs '+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs' -t_srs EPSG:3005 -r cubicspline -tr 25 -25 -srcnodata -999 -dstnodata -32768 -cutline ./watersheds/smooth.shp -crop_to_cutline /Users/kylechezik/Documents/BC_DEM25/bcdem.tif DEM.tif
	rm burn.tif

# Clip forest change year/type to region of concern.
clip_forest_year:
	unzip ../../Forest_Change/C2C_Change_Type.zip -d ../../Forest_Change/
	unzip ../../Forest_Change/C2C_Change_Year.zip -d ../../Forest_Change/
	gdalwarp -overwrite -s_srs '+proj=lcc +lat_1=49 +lat_2=77 +lat_0=49 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs' -t_srs EPSG:3005 -r near -tr 25 -25 -srcnodata 0 -dstnodata -32768 -te 1302474.053 581677.768 1554024.053 890452.768 ../../Forest_Change/C2C_change_year.tif forest_year_extent.tif
	gdalwarp -overwrite -cutline ./watersheds/smooth.shp forest_year_extent.tif forest_year.tif
	rm forest_year_extent*

# Legend: 1 = Wildfire, 2 = Harvest, 5 = Lower Condifence Wildfire, 6 = Lower Confidence Harvest
clip_forest_type:
	gdalwarp -overwrite -s_srs '+proj=lcc +lat_1=49 +lat_2=77 +lat_0=49 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs' -t_srs EPSG:3005 -r near -tr 25 -25 -srcnodata 0 -dstnodata -32768 -te 1302474.053 581677.768 1554024.053 890452.768 ../../Forest_Change/C2C_change_type.tif forest_type_extent.tif
	gdalwarp -overwrite -cutline ./watersheds/smooth.shp forest_type_extent.tif forest_type.tif
	rm forest_type_extent*

# Separate wildfire from logging.
wildfire:
	gdal_calc.py -A forest_type.tif --calc "A*(A==1)|A*(A==5)" --NoDataValue 0 --overwrite --outfile wildfire.tif
logging:
	gdal_calc.py -A forest_type.tif --calc "A*(A==2)|A*(A==6)" --NoDataValue 0 --overwrite --outfile logging.tif

# Extract wildfire and logging by year.
wildfire_year:
	gdal_calc.py -A wildfire.tif -B forest_year.tif --calc "B*(A==1)" --NoDataValue 0 --overwrite --outfile wildfire_year_1.tiff
	gdal_calc.py -A wildfire.tif -B forest_year.tif --calc "B*(A==5)" --NoDataValue 0 --overwrite --outfile wildfire_year_2.tiff

logging_year:
	gdal_calc.py -A logging.tif -B forest_year.tif --calc "B*(A==2)" --NoDataValue 0 --overwrite --outfile logging_year_1.tiff
	gdal_calc.py -A logging.tif -B forest_year.tif --calc "B*(A==6)" --NoDataValue 0 --overwrite --outfile logging_year_2.tiff

cleanup_type:
	rm logging.ti* wildfire.ti* forest*

#setopt extended_glob may need to set this up to use extended glob in zsh.
clean_zip:
	rm -- ../../Forest_Change/^*.zip

# Clip waterbodies and rasterize
clip_rasterize_waterbodies:
	saga_cmd shapes_polygons 11 -CLIP "/Users/kylechezik/sfu/River_Network_Temperature_Output/shp_watersheds/watershed.shp" -S_INPUT "/Users/kylechezik/sfuvault/Data/Original_Data/GIS_Data/BC_Lakes/FWLKSPL/FWLKSPL_polygon.shp" -S_OUTPUT "/Users/kylechezik/sfu/River_Network_Temperature_Output/lakes/vec_lakes.shp" -MULTIPLE "0"
#	gdal_rasterize -l vec_lakes -at -a wbid -a_nodata -32768 -te 1302474.053 581677.768 1554024.053 890452.768 -a_srs EPSG:3005 -tr 25 -25 -ot Int16 /Users/kylechezik/sfu/River_Network_Temperature_Output/lakes/vec_lakes.shp /Users/kylechezik/sfu/River_Network_Temperature_Output/lakes/lakes.tif
#	rm /Users/kylechezik/sfu/River_Network_Temperature_Output/lakes/vec_l*

# Clip glacier and rasterize
clip_rasterize_glaciers:
	ogr2ogr -f "ESRI Shapefile" -s_srs EPSG:4326 -t_srs EPSG:3005 -overwrite vec_glacier1.shp ../../Glacier_Data_Clark_etal_2015/02_rgi60_WesternCanadaUS/02_rgi60_WesternCanadaUS.shp
	saga_cmd shapes_polygons 11 -CLIP "./watersheds/watershed.shp" -S_INPUT "vec_glacier1.shp" -S_OUTPUT "vec_glaciers.shp" -MULTIPLE "0" -DISSOLVE "0"
	gdal_rasterize -l vec_glaciers -at -burn 625 -a_nodata -32768 -te 1302474.053 581677.768 1554024.053 890452.768 -tr 25 -25 vec_glaciers.shp glaciers.tif
	rm vec_g*

# Rasterize climate data. (example: `make rasterize_climate VAR=MAT FILE=climateBC_2014-2017`)
rasterize_climate:
# Create vrt file
	echo "<OGRVRTDataSource>\n<OGRVRTLayer name=\"$(FILE)\">\n<SrcDataSource>/Users/kylechezik/sfu/River_Network_Temperature_Output/predictors/climate_csv/"$(FILE)".csv</SrcDataSource>\n<GeometryType>wkbPoint</GeometryType>\n<GeometryField encoding=\"PointFromColumns\" x=\"Longitude\" y=\"Latitude\" z=\"$(VAR)\"/>\n</OGRVRTLayer>\n</OGRVRTDataSource>" > clim.vrt
# Rasterize
	gdal_grid -a nearest:radius1=0.01:radius2=0.01:angle=0:nodata=-32768 -zfield "$(VAR)" -txe -121.756 -117.763 -tye 49.9807 52.9306 -outsize 252 309 -of GTiff -a_srs EPSG:4326 -ot FLoat64 -l $(FILE) clim.vrt ~/sfu/River_Network_Temperature_Output/predictors/$(FILE)_$(VAR).tiff
# Downscale and reproject
	gdalwarp -s_srs EPSG:4326 -t_srs EPSG:3005 -tr 25 25 -r cubicspline -overwrite -srcnodata -32768 -dstnodata -32768 /Users/kylechezik/sfu/River_Network_Temperature_Output/predictors/$(FILE)_$(VAR).tiff /Users/kylechezik/sfu/River_Network_Temperature_Output/predictors/$(FILE)_$(VAR)_25.tiff
# Clip climate raster by shapefile
	gdalwarp -q -cutline ~/sfu/River_Network_Temperature_Output/shp_watersheds/smooth.shp -tr 25 25 /Users/kylechezik/sfu/River_Network_Temperature_Output/predictors/$(FILE)_$(VAR)_25.tiff /Users/kylechezik/sfu/River_Network_Temperature_Output/predictors/C$(FILE)_$(VAR).tiff
# Remove old TIFF and CSV files
	rm /Users/kylechezik/sfu/River_Network_Temperature_Output/predictors/$(FILE)_$(VAR)*

# Difference seasonal min/max climate data for annual amplitude
amp_diff:
	gdal_calc.py -A $(MIN) -B $(MAX) --outfile=$(YEAR) --calc="B-A" 
