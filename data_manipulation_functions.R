# load packages
library("purrr")
library("kknn")


# degrees-minutes-seconds string to decimal degrees
toDecimal_degrees = function(dms_string,split=" ")
{
  if(is.na(dms_string))
  {
    return(NA)
  }
  else
  {
      numbers = as.numeric(strsplit(dms_string,split=split)[[1]])

      if (numbers[1]<0)
      {
        decimal_degrees = numbers[1]-(numbers[2]/60)-(numbers[3]/3600)
      }
      else
      {
        decimal_degrees = numbers[1]+(numbers[2]/60)+(numbers[3]/3600)
      }
  
      return(decimal_degrees)
  }
}

# vectorized toDecimal_degrees
toDecimal_degrees_v = Vectorize(toDecimal_degrees)

# reference + azimuth and distance to coordinates of new point
azimuth_coords = function(distance,azimuth,ref_x,ref_y)
{
  delta_x = distance * sin(azimuth)
  delta_y = distance * cos(azimuth)
  
  out_x = ref_x + delta_x
  out_y = ref_y + delta_y
  
  return (c(out_x,out_y))
}

# vectorized azimuth_coords
azimuth_coords_v = function(distance_v,azimuth_v,ref_x_v,ref_y_v)
{

  args = list(distance_v,azimuth_v,ref_x_v,ref_y_v)
  
  out_coords = args %>% pmap(azimuth_coords)
  
  out_coords <- data.frame(matrix(unlist(out_coords), nrow=length(distance_v), byrow=T))
  
  names(out_coords)=c("x","y")

  return(out_coords)
}

# find nearest (spatial) neighbours
k_nearest = function(xy_df,k_nearest=1,distances=TRUE)
{
  xy_df = data.frame(id=1:nrow(xy_df),xy_df)
  
  nearest <- kknn(id~., train=xy_df, test=xy_df, distance = 2,
                  kernel = "rectangular",k=k_nearest+1) 

  if (distances)
  {
    # distance to nearest neighbour
    return(nearest$D[,2:(k_nearest+1)])
  }
  else
  {
    # index of nearest neighbour
    return(nearest$C[,2:(k_nearest+1)])
  }
}

# find angle between two 2d-vectors
angle = function(x,y)
{
  x = as.numeric(x)
  y = as.numeric(y)
  
  theta= as.numeric(acos( sum(x*y) / ( sqrt(sum(x * x)) * sqrt(sum(y * y)) ) ))

  return(theta)
}


#calculate area of triangles in a Delauny tesselation matrix and the origin coordinates plus a weight variable (e.g. tree height)
triangles_area = function(d,p2)
{
  area <- vector(,1)
  Kreuz <- matrix(NA, nrow=length(d[,1]),ncol=4)
  EVector<- matrix(NA, nrow=length(d[,1]),ncol=3)
  for (n in 1:length(d[,1])) {  
    k<-p2[d[n,],]                  
    ifelse(k[1,1]==k[2,1]&k[1,1]==k[3,1], k<-NA,{   
      ifelse(k[1,2]==k[2,2]&k[1,2]==k[3,2], k<-NA,{ 
        
        ab<-k[1,]-k[2,]  
        ac<-k[1,]-k[3,]  
        bc<-k[2,]-k[3,]  
        a<-(ab[1]^2+ab[2]^2+ab[3]^2)^0.5 
        b<-(ac[1]^2+ac[2]^2+ac[3]^2)^0.5 
        c<-(bc[1]^2+bc[2]^2+bc[3]^2)^0.5 
        s<-(a+b+c)/2                      
        area[n]<-(s*(s-a)*(s-b)*(s-c))^0.5  
        
        Kreuz[n,1]<-(ab[2]*ac[3]-ab[3]*ac[2])   
        ifelse((ab[1]*ac[2]-ab[2]*ac[1]) < 0, Kreuz[n,2]<--(ab[3]*ac[1]-ab[1]*ac[3]), Kreuz[n,2]<-(ab[3]*ac[1]-ab[1]*ac[3])) 
        ifelse((ab[1]*ac[2]-ab[2]*ac[1]) < 0, Kreuz[n,3]<-abs((ab[1]*ac[2]-ab[2]*ac[1])), Kreuz[n,3]<-(ab[1]*ac[2]-ab[2]*ac[1])) 
        
        Kreuz[n,4]<-((ab[2]*ac[3]-ab[3]*ac[2])^2 + (ab[3]*ac[1]-ab[1]*ac[3])^2 + (ab[1]*ac[2]-ab[2]*ac[1])^2)^0.5 
        EVector[n,]<-Kreuz[n,1:3]/Kreuz[n,4] 
      })})
  }
  output = list(evector=EVector,area=area)
  return(output)
}

# build abundance vector from potential species list and an observed species vector TODO
#abundance_vector = function(species_list,observed_species)
#{
  
#}
