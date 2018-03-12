
# load packages
library("readxl")
library("sp")
library("raster")
library("rgdal")
library("dplyr")
library("purrr")
library("sf")
library("readr")

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

# unique conglomerates
conglomerates = unique(data_table$UPMID)

# build species list
species_list = paste0(data_table$genero," ",data_table$especie)
species_list = unique(species_list[species_list!="zz-sin genero ZZ-sin especie"])

# initialize output matrix
output_matrix = data.frame(matrix(0,length(conglomerates)*4,17))
colnames(output_matrix) = c("UPMID","Sitio","lon","lat","eberhardt","clarkevans",
                            "simpson","shannon","meanRenyi","richness","eveness","fAlpha",
                            "gMingle","gDomi","SCI","ESCI_1","ESCI_2")

# intialize counter
contador = 0

for (i in 1:length(conglomerates))
{
  print(i)

  for (j in 1:4)
  {
  
  contador = contador+1
  
  # conglomerado id
  output_matrix[contador,1]=conglomerates[i]
  
  # sitio
  output_matrix[contador,2]=j
  
  # choose one site in one conglomerate
  site = data_table %>% filter(UPMID==conglomerates[i],Sitio==j)
  
  if(nrow(site)!=0)
  {
  
  # to a spatial object
  coordinates(site)=~lon+lat
  proj4string(site)= "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
  
  # project site coordinates
  site_reproj = spTransform(site,
                            CRS("+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +datum=WGS84 +units=m +no_defs"))
  
  # obtain coordinates per tree (using site center coords, azimuth and distance to each tree)
  trees_coords = azimuth_coords_v(site$distancia,site$azimut,coordinates(site_reproj)[,1],coordinates(site_reproj)[,2])
  
  # new tree's data.frame
  trees_data_table = data.frame(treeid=1:nrow(trees_coords),trees_coords,as.data.frame(site_reproj))
  
  # lon
  output_matrix[contador,3]=coordinates(site_reproj)[1,1]
  
  # lat
  output_matrix[contador,4]=coordinates(site_reproj)[1,2]
  
  tree_species = paste0(trees_data_table$genero," ",trees_data_table$especie)
  
  # remove unknown species
  tree_species_clean = tree_species[tree_species!="zz-sin genero ZZ-sin especie"]
  
  species_matrix = data.frame(species_list,abundance=rep(0,length(species_list)),stringsAsFactors = FALSE)
  
  abundances = table(tree_species_clean)
  
  species_matrix[match(names(abundances),species_list),2]= abundances
  
  output_matrix[contador,7]=calc_simpson(species_matrix[,2])
  
  output_matrix[contador,8]=calc_shannon(species_matrix[,2])
  
  output_matrix[contador,9]=mean(calc_renyi(species_matrix[,2]),na.rm=TRUE)
  
  output_matrix[contador,10]=calc_richness(species_matrix[,2])
  
  output_matrix[contador,11]=calc_eveness(species_matrix[,2])
  
  output_matrix[contador,12]=calc_fishersAlpha(species_matrix[,2])
  
  if(nrow(site_reproj)>=6)
  {
  
  ##### populate indicators at site level

  output_matrix[contador,5]=calc_eberhardt(trees_coords)
  
  output_matrix[contador,6]=calc_clarkevans(trees_coords,trees_data_table[1,c("lon","lat")])

  output_matrix[contador,13]=calc_gadowsming(trees_coords,tree_species)

  output_matrix[contador,14]=calc_gadowsdomi(trees_coords,trees_data_table$altura_total)

  esci_list = calc_esci(trees_coords,trees_data_table$diametro_normal)

  if(!is.na(esci_list))
  {
    
    output_matrix[contador,15]=esci_list$SCI
    output_matrix[contador,16]=esci_list$ESCI_1
    output_matrix[contador,17]=esci_list$ESCI_2
  }
  
  }
  
  }
  
  
  }
  
}

for (i in 1:ncol(output_matrix))
{
  output_matrix[is.nan(output_matrix[,i]),i] =0
  output_matrix[is.na(output_matrix[,i]),i] =0
}

write_csv(output_matrix,"./data/sites_indicators.csv")
