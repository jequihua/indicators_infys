# load packages
library("readxl")
library("sp")
library("raster")
library("rgdal")
library("dplyr")
library("purrr")
library("sf")

# loud source functions
source("./indices.R")
source("./data_manipulation_functions.R")

# display excel sheets
excel_sheets("./data/Arbolado_D_2015_entregable.xlsx")

# read in excel sheet
data_table = read_excel("./data/Arbolado_D_2015_entregable.xlsx")

# select only relevant variables
data_table =  data_table %>% select(UPMID,Sitio,latitud,longitud,`Fasesucesional/tipovegetación`,
                                    azimut,distancia,familia,genero,especie,nombre_comun,diametro_normal,
                                    altura_total,diametro_copa_este_oeste,diametro_copa_norte_sur)

# add site counts variable
data_table = data_table %>% group_by(UPMID) %>% mutate(Sitio_counts = n_distinct(Sitio))

# filter trees in conglomerates with less than 4 
data_table = data_table %>% filter(Sitio_counts == 4)

# add latitud and longitud variables in decimal degrees
data_table = data_table %>% mutate(lat=toDecimal_degrees_v(latitud)) 
data_table = data_table %>% mutate(lon=toDecimal_degrees_v(longitud)) 

# clean missing coordinates
data_table = data_table[complete.cases(data_table[,c("lon","lat")]),]

# to a spatial object
coordinates(data_table)=~lon+lat
proj4string(data_table)= "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

# project site coordinates
data_table_reproj = spTransform(data_table,
                          CRS("+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +datum=WGS84 +units=m +no_defs"))

# obtain coordinates per tree (using site center coords, azimuth and distance to each tree)
trees_coords = azimuth_coords_v(data_table$distancia,data_table$azimut,coordinates(data_table_reproj)[,1],coordinates(data_table_reproj)[,2])

# output tree data with coordinates
out_data = data.frame(data_table_reproj@data,arbol_lon=trees_coords$x,arbol_lat=trees_coords$y,stringsAsFactors = FALSE)

coordinates(out_data)=~arbol_lon+arbol_lat
proj4string(out_data) = CRS("+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +datum=WGS84 +units=m +no_defs")

# to lat lon
out_data = spTransform(out_data,
                       CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))

# write out a new shapefile (including .prj component)
writeOGR(out_data, "./data/arboles_infys_2015.shp", "arboles_infys_2015", driver="ESRI Shapefile")
