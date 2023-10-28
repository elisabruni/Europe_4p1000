library(raster)
library(rgdal)

#############
#IMPORT DATA
#############
#Import raster data
str_name<-'1out_SOC_er_gaec_Nd_lvs_4.5_ch4_v18_e22.tif' 
file_name <- paste0('/Users/ebruni/Desktop/DOTTORATO/DATI/LUCASSOC/NPP/',str_name)
imported_raster=raster(file_name,band=10) #Cinput(tC/ha/y), average 2015-2019

#imported_raster=raster(file_name,band=9) #NPP (tC/ha/y), average 2015-2019

#Import point data you want to extract
lonlat_file_name <-paste0('lonlat_crop_sites_LUCAS.txt')
lonlat_file <-read.table(lonlat_file_name, header=TRUE)

#################
#VISUALIZE RASTER
#################
#See info on raster
print(imported_raster)
val_raster <- getValues(imported_raster)
min(val_raster,na.rm=TRUE)
max(val_raster,na.rm=TRUE)
sum(val_raster,na.rm=TRUE)

# val_raster_old <- getValues(imported_raster_old)
# min(val_raster_old,na.rm=TRUE)
# max(val_raster_old,na.rm=TRUE)
# sum(val_raster_old,na.rm=TRUE)

#Plot raster
plot(imported_raster)
#Plot histogram raster
hist(imported_raster,maxpixels=22000000)

################
#EXTRACT POINTS
################
#Translate extent (meters) into lon-lat coordinates
proj_raster <- projectRaster(imported_raster, crs='+proj=longlat +datum=WGS84')

val_raster_proj <- getValues(proj_raster)
min(val_raster_proj,na.rm=TRUE)
max(val_raster_proj,na.rm=TRUE)
sum(val_raster_proj,na.rm=TRUE)

#Transform lon-lat file into a spatial object
coordinates(lonlat_file)=~LONG+LAT
mypoints = SpatialPoints(lonlat_file,proj4string = CRS('+proj=longlat +datum=WGS84'))
#Extract values from raster at mypoints
myvalues = extract(proj_raster, mypoints)
#Add data to lon-lat file
#lonlat_file$Cin <- myvalues
lonlat_file$NPP <- myvalues
#Plot raster
plot(proj_raster)
#Plot location of the lon-lat points
points(lonlat_file$LONG,lonlat_file$LAT, pch=0, cex = 0.1,col="black")

write.csv(lonlat_file,'/Users/ebruni/Desktop/DOTTORATO/DATI/LUCASSOC/LUGATO/Cinput_Daycent_extracted.csv')
#write.csv(lonlat_file,'/Users/ebruni/Desktop/DOTTORATO/DATI/LUCASSOC/LUGATO/NPP_Daycent_extracted.csv')

#SOURCES
#https://www.neonscience.org/resources/learning-hub/tutorials/extract-values-rasters-r
#https://proj.org/operations/projections/utm.html
#https://www.rdocumentation.org/packages/raster/versions/3.5-15/topics/projectRaster
#https://gisday.wordpress.com/2014/03/24/extract-raster-values-from-points-using-r/comment-page-1/
#https://datacarpentry.org/r-raster-vector-geospatial/03-raster-reproject-in-r/
