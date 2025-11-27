# Load Libraries and Data ----
library(sf)
library(stringr)
library(dplyr)
library(here)
library(ggplot2)
library(ggthemes)
library(parallel)
library(RColorBrewer)
library(cowplot)
library(ggforce)
library(deldir)
library(units)
library(igraph)
library(rnaturalearth)

source(here('src','hex_areas.R'))

#-------------------------------------------------------------------------------
# Generate Spatial Window for Analyses: Box same size as Europe ----

#-------------------------------------------------------------------------------
##Sampling window: Europe ----
#List of countries in europe (excluding france and norway which need to be added separately if needed) ---
europe_countries <- c("Albania", "Latvia", "Andorra", "Liechtenstein",
                      "Armenia", "Lithuania", "Austria", "Luxembourg",
                      "Azerbaijan", "Malta", "Belarus", "Moldova",
                      "Belgium", "Monaco", "Bosnia and Herzegovina", "Montenegro",
                      "Bulgaria", "Netherlands", "Croatia", 
                      "Cyprus", "Poland", "Czechia", "Portugal",
                      "Denmark", "Romania", "Estonia","Finland", "San Marino",
                      "Macedonia", "Republic of Serbia", "Slovakia",
                      "Georgia", "Slovenia", "Germany", "Spain", "Greece", "Sweden",
                      "Hungary", "Switzerland", "Ireland", "Turkey",
                      "Italy", "Ukraine", "Kosovo", "United Kingdom")
#-------------------------------------------------------------------------------
##Generate Spatial Window for Analyses
#Sampling window: Bounding box roughly the centered on mainland Europe ----

sampling_europe_win <- ne_countries(country=europe_countries, returnclass = "sf")

#Create bounding box in correct CRS
sampling_win_proj <- st_transform(sampling_europe_win, crs = 3035) #project to CRS for LAEA Europe
sampling_win_proj <- st_as_sfc(st_bbox(sampling_win_proj))

hex_area_win_proj <- hex_areas(sampling_win_proj, cell_d = 600000) #hex_areas of 600km in diameter #for mapping: CRS for LAEA Europe (crs=3035)

##Plot sample window and hex areal units
ggplot(data = hex_area_win_proj) +
  #geom_sf(data = sampling_win_proj, aes(color = "grey50"), lwd=1) + #See sampling_win
  geom_sf() + #hex grid
  geom_sf(data = hex_area_win_proj$area_center, color="purple", size=1.5, alpha=0.7) +
  geom_sf_text(data = hex_area_win_proj$area_center, aes(label=hex_area_win_proj$area_ID), size=4, color="black") +
  coord_sf(crs=st_crs(3035)) +
  theme(panel.background = element_rect(fill = "lightblue",
                                        colour = "lightblue",
                                        size = 0.5,
                                        linetype = "solid"),
        legend.position = "none")

#-------------------------------------------------------------------------------
##Delaunay triangulation between hex centers --- 
#Note: should not be done in WGS84, because deldir() assumes Euclidean geometry
hex_area_centers_proj <- st_as_sf(hex_area_win_proj$area_center)
center_coords <- st_coordinates(hex_area_centers_proj) 
del <- deldir(center_coords, id=hex_area_win_proj$area_ID)
tiles <- tile.list(del)

##Add center_coords in constants
constants_trig <- list()
constants_trig$center_coords <- center_coords

##Remove external edges that are outside sampling window ---
#Create edges as linestrings
st_segment <- function(r){st_linestring(t(matrix(unlist(r), 2, 2)))} #function to make LINESTRINGs from Delaunay segments
edges_coords <- del$delsgs[ , 1:4]
edges_coords$geom <- st_sfc(sapply(1:nrow(edges_coords), 
                                   function(i){st_segment(edges_coords[i,])}, simplify=FALSE), crs = 3035) #create sf object of segments (in projected CRS)


#Filter edges to those within the sampling window (projected) 
check_edges_ext <- st_within(edges_coords$geom, st_buffer(sampling_win_proj, dist=60000)) #check if edges are inside the (buffered) sampling window
ind_ext_edges <- which(is.na(as.numeric(check_edges_ext))) # indices of edges which are external -- to be removed if we aren't considering ocean movement
##Check
if (length(ind_ext_edges > 0)){
    internal_edges <- edges_coords[-ind_ext_edges, ]
  } else {
    internal_edges <- edges_coords} #if there are no external edges, retain all the original edges
plot(sampling_win_proj)
plot(internal_edges$geom, add=T) #compare to plot(edges_coords$geom, add=T)


##Plot delaunay triangulation ---
ggplot() +
  geom_sf(data = st_buffer(sampling_win_proj, 40000), fill=NA, color = "grey50", linewidth=1) + #sampling window with coastal buffer
  geom_sf(data = hex_area_win_proj) + #hex grid
  geom_sf_text(data = hex_area_win_proj, aes(label = area_ID), size=4, alpha=0.8) + #hex grid labels 
  geom_sf(data = hex_area_win_proj$area_center, size=2, alpha=1, color = "purple") + #hex-origins
  geom_delaunay_segment(aes(x=center_coords[,1], y=center_coords[,2]),
                        alpha=0.3,
                        colour='purple',
                        size=0.8) +
  labs(x = "Longitude", y = "Latitude") +
  theme(panel.background = element_rect(fill = "lightblue",
                                        colour = "lightblue",
                                        size = 0.5,
                                        linetype = "solid"),
        legend.position = "none")

#-------------------------------------------------------------------------------
# ##Compute Great-Arc Distances in km between area centers ---
# #We take the center of each area k to be a point representing that whole area
# hex_area_win <- st_transform(hex_area_win_proj, crs = 4326) ##returns WGS84 projection -- NB: to be used for calculations
# sampling_win <- st_transform(sampling_win_proj, crs = 4326)
# hex_area_centers <- st_transform(hex_area_centers_proj, crs = 4326)
# hex_dist_mat <- set_units(st_distance(hex_area_centers), 'km') #Inter-area distance matrix in km: each area's distance from every other area.

##Compute Distances in km between area centers ---
hex_dist_mat  <- set_units(st_distance(hex_area_centers_proj), 'km') #in crs=3035

#-------------------------------------------------------------------------------
##Transitions dataframe --
transitions <- del$delsgs %>% 
  dplyr::select(region1_id = ind1, 
                region2_id = ind2,
                region1_x = x1,
                region1_y = y1,
                region2_x = x2,
                region2_y = y2) %>% 
  mutate(region1_id = as.integer(region1_id), 
         region2_id = as.integer(region2_id)) %>% 
  rowwise() %>% 
  mutate(distance = hex_dist_mat[region1_id, region2_id]) %>% #Great-arc distance between transitions in km
  filter(distance < set_units(610, "km")) #remove transitions between hexs which aren't adjacent #note: change depending on chosen hex size 

if (length(ind_ext_edges > 0)){
  transitions <- transitions[-ind_ext_edges,] #remove transitions outside the sampling window
} else {
  transitions <- transitions} #if there are no external edges, retain all the transitions

#Transform transitions into usable format to save in constants
edge_info <- as.data.frame(transitions)
constants_trig$n_trans <- nrow(transitions) #Number of transitions
constants_trig$edge_id1 <- edge_info$region1_id
constants_trig$edge_id2 <- edge_info$region2_id 
constants_trig$edge_dist <- edge_info$distance


#-------------------------------------------------------------------------------
## Save everything on a R image file ----
#Save transitions between hex areas
save(#del, 
     #tiles,
     edge_info,
     constants_trig,
     file=here('data','trig.RData')) 

#Save sampling window separately
constants_sw<- list()
constants_sw$countries <- c(europe_countries, "France", "Norway")
constants_sw$n_areas  <- nrow(hex_area_win_proj)

save(constants_sw, sampling_win_proj, hex_area_win_proj, file=here('data','sample_window.RData'))