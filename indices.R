
# load packages
library("vegan")
library("spatstat")
library("geometry")

########################################################## Indices that depend on tree position and distribution
##########################################################

########################################################## Indices based on (tree) nearest neighbours

# Eberhardt index
calc_eberhardt = function(tree_coords)
{
  # tree coords expected as (x,y) columns
  distances_to_nearest = k_nearest(tree_coords)
  
  eberhardt = mean(distances_to_nearest^2,na.rm=TRUE)/mean(distances_to_nearest,na.rm=TRUE)^2
  
  return(eberhardt)
}

# Clark-Evans index
calc_clarkevans =  function(tree_coords,site_centre_coords,site_radius=13)
{
  owin_disc = disc(radius=site_radius,centre=as.numeric(site_centre_coords))
  
  # tree coords expected as (x,y) columns
  point_pattern = ppp(tree_coords[,1],tree_coords[,2],window=owin_disc,check=TRUE)
  
  clarkevans = clarkevans(point_pattern,correction="none")
  
  return(clarkevans)
}

# Ripley's K-function (is not an index)
calc_ripleysk = function(tree_coords,site_centre_coords,site_radius=11.28)
{
  owin_disc = disc(radius=site_radius,centre=as.numeric(site_centre_coords))
  
  # tree coords expected as (x,y) columns
  point_pattern = ppp(tree_coords[,1],tree_coords[,2],window=owin_disc,check=TRUE)
  
  ripleysk = Kest(X=point_pattern, correction="isotropic")
  
  return(ripleysk)
}

# test for Complete Spatial Randomness (CRS) (is not an index)
# returns p.value of test
calc_CRStest = function(tree_coords,site_centre_coords,site_radius=11.28,verbose=FALSE)
{
  owin_disc = disc(radius=site_radius,centre=as.numeric(site_centre_coords))
  
  # tree coords expected as (x,y) columns
  point_pattern = ppp(tree_coords[,1],tree_coords[,2],window=owin_disc,check=TRUE)
  
  CRStest = dclf.test(point_pattern, Lest, nsim=99,verbose=verbose)

  return(CRStest$p.value)
}

########################################################## Species diversity and mix indices
##########################################################

########################################################## Species diversity

# Simpsons index
calc_simpson = function(tree_sp_vector,base = exp(1))
{
  simpson = diversity(tree_sp_vector, "simpson",base = base)
  
  return(simpson)
}

# Shannons index
calc_shannon = function(tree_sp_vector)
{
 shannon = diversity(tree_sp_vector, "shannon")
 
 return(shannon)
}

# Renyi diversity index
calc_renyi = function(tree_sp_vector, scales=c(0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, Inf))
{
  renyi = renyi(tree_sp_vector,scales=scales)
  
  return(renyi)
}

# Species richness
calc_richness = function(tree_sp_vector)
{
  richness = specnumber(tree_sp_vector)
  
  return(richness)
}

# Species eveness
calc_eveness = function(tree_sp_vector)
{
  eveness = diversity(tree_sp_vector)/log(specnumber(tree_sp_vector))
  
  return(eveness)
}

# Fishers alpha 
calc_fishersAlpha = function(tree_sp_vector)
{
  fishersAlpha = fisher.alpha(tree_sp_vector)
  
  return(fishersAlpha)
}

########################################################## Species and attributes mix

# Gadows mingling index
calc_gadowsming = function(tree_coords,tree_species,number_of_neighbours=4)
{
  if (nrow(tree_coords)!=length(tree_species))
  {
    print("number of coordinates and tree species must match!")
  }
  else
  {
  k_nearest_each = k_nearest(tree_coords,k_nearest=number_of_neighbours,distance=FALSE)

  percentage_of_matches = list()
  
  for (i in 1:nrow(k_nearest_each))
  {
    
    neighbour_trees = tree_species[k_nearest_each[i,]]
    
    percentage_of_matches[[i]] = sum(neighbour_trees == tree_species[i])/number_of_neighbours
    
  }
  
  gadowsming =  mean(unlist(percentage_of_matches))
  
  return(gadowsming)
  
  }
}

# Gadows dominance index
calc_gadowsdomi = function(tree_coords,tree_sizes,number_of_neighbours=4)
{
  if (nrow(tree_coords)!=length(tree_sizes))
  {
    print("number of coordinates and tree sizes must match!")
  }
  else
  {
    k_nearest_each = k_nearest(tree_coords,k_nearest=number_of_neighbours,distance=FALSE)
    
    percentage_of_matches = list()
    
    for (i in 1:nrow(k_nearest_each))
    {
      
      neighbour_trees = tree_sizes[k_nearest_each[i,]]
      
      percentage_of_matches[[i]] = sum(neighbour_trees == tree_sizes[i])/number_of_neighbours
      
    }
    
    gadowsdomi =  mean(unlist(percentage_of_matches))
    
    return(gadowsdomi)
    
  }
}


# Gadows uniform angle index TODO
#calc_gadowsunif = function(tree_coords,tree_sizes,number_of_neighbours=4)
#{

#}

########################################################## Structural complexity indices
##########################################################

# ESCI indices
calc_esci = function(tree_coords,tree_sizes,tree_size_threshold=7.5,circular_plot_radius=12)
{
  tree_coords[is.na(tree_coords[,1]),1]=mean(tree_coords[,1],na.rm=TRUE)
  tree_coords[is.na(tree_coords[,2]),2]=mean(tree_coords[,2],na.rm=TRUE)
  tree_sizes[is.na(tree_sizes)]=mean(tree_sizes,na.rm=TRUE)

  if (nrow(tree_coords)!=length(tree_sizes))
  {
    print("number of coordinates and tree sizes must match!")
  }
  else
  {
  tree_coords = tree_coords[tree_sizes>=tree_size_threshold,]
  tree_sizes = tree_sizes[tree_sizes>=tree_size_threshold]
  
  new_tree_data = as.matrix(data.frame(x=tree_coords[,1],
                                       y=tree_coords[,2],
                                       tree_sizes=tree_sizes))

  # Delaunay tesselation
  if(nrow(tree_coords)>=3)
  {
    tesselation = delaunayn(new_tree_data[,1:2],options="QbB")

    triangles = triangles_area(tesselation,new_tree_data)
    
    EVector = triangles$evector
    area = triangles$area
    
    EVector[is.na(EVector)] = 0 
    area = subset(area, area > 0)
    
    
    # calculate Vector Ruggedness Measure (vrm) 
    vrm = ((sum(EVector[,1])^2 + sum(EVector[,2])^2 + sum(EVector[,3])^2)^0.5)/length(area) 
    vrm1 = 2-vrm
    
    area_length = length(area) # Number of triangles in TIN
    tin_area =  sum(area) # TIN area
    proj_area = convhulln(new_tree_data[,1:2], option="FA")$vol # projected TIN area
    vector_rugg = vrm1 # vector ruggedness measure
    SCI = sum(area)/proj_area # Structural Complexity Index (SCI): Ratio of 'real' area of triangles and projected area of triangles.
    ESCI_1 = SCI*vrm1   # ESCI' = SCI * vrm
    ESCI_2 = SCI*vrm1*(1+area_length*10/(pi*circular_plot_radius^2)) # ESCI = SCI * vrm * (1 + trees per 10 ma^2)

    esci = list(area_length=area_length,
                tin_area=tin_area,
                proj_area=proj_area,
                vector_rugg=vector_rugg,
                SCI=SCI,
                ESCI_1=ESCI_1,
                ESCI_2=ESCI_2)
    
  }
  
  else
  {
    esci=NA
  }
  
  return(esci)
  }
}


