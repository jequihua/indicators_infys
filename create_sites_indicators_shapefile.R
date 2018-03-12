
# load packages
library("dismo")
library("rgeos")

# initialize sites (circles) list
sites = list()

for (i in 1:nrow(output_matrix))
{

   site_centre = output_matrix[i,c("lat","lon")]

   coordinates(site_centre)=~lon+lat
   
   projection(site_centre)=CRS("+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +datum=WGS84 +units=m +no_defs")
   

   # 12 m
   circle = polygons(circles(site_centre,12, lonlat=FALSE))
   
   sites[[i]]=circle
}

# join all sites (circles)
joined = do.call(bind,sites) 


# associate indicators data.frame
joined = SpatialPolygonsDataFrame(Sr=joined,
                                  data=output_matrix,
                                  FALSE)

coordinates(joined)

projection(joined)=CRS("+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +datum=WGS84 +units=m +no_defs")

joined2 = spTransform(joined,CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))

# escribimos lo anterior a disco
writeOGR(joined,
         "./data/sites_indicators3.shp",
         "sites_indicators3",
         driver="ESRI Shapefile",
         overwrite_layer=TRUE)
