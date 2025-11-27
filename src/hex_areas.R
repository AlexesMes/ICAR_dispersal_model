hex_areas <- function(win, cell_d = 400000){
  require(dplyr)
  require(stringr)
  require(rnaturalearth)
  require(sf)
  
  #Spatial Window for Analyses ----
  sampling_win <- win
  
  #Buffer polygon ----
  buffer_win <- st_buffer(sampling_win, dist=60000) #60km buffer to ensure all sites fall in sample_win
  
  #Generate Spatial Hexagons ----
  cell_diameter <- cell_d
  hex_grid <- st_make_grid(buffer_win, square=FALSE,  cellsize = cell_diameter) #makes an hexagonal grid (default: what = "polygons)
  
  #Clip to sampling window with coastal buffer ----
  hex_grid_clipped <- st_intersection(st_as_sf(hex_grid), buffer_win)
  
  #Drop sliver polgons ----
  hex_grid_clipped <- st_collection_extract(hex_grid_clipped, "POLYGON")
  
  #Assign hex IDs ----
  hex_grid_clipped <- hex_grid_clipped %>%
    rename(geometry = x) %>%
    mutate(area_ID = row_number(),
           area_center = st_centroid(hex_grid_clipped))
  
  # ##Transform back to WGS84 ----
  # hex_grid_clipped <- st_transform(hex_grid_clipped, crs = 4326)
  
  #Return Output ----
  return(hex_grid_clipped)
  
}
