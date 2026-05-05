# Load Libraries and spatial data ----
library(here)
library(dplyr)
library(tidyr)
library(stringr)
library(truncnorm)
library(cascsim)
library(corrplot)
library(ggplot2)
library(ggridges)
library(rnaturalearth)
library(nimbleCarbon)
library(rcarbon)
library(sf)
library(viridis)
library(cowplot)
library(wesanderson)
library(latex2exp)
library(gridExtra)
library(grid)
library(gridBase)
library(diagram)
library(coda)
library(graphics)
library(ggthemes)
library(ggpattern)
library(rlist)
library(patchwork)
library(tibble)

rm(list = ls())
`%!in%` <- Negate(`%in%`)

#===============================================================================
#Set-up
barcolours1 = c("skyblue","dodgerblue","darkblue","darkgreen")
barcolours2 = c("plum","orchid","purple","darkgreen")

#===============================================================================
## Load Data ----

source(here('src', 'grad_funcs.R'))
source(here('src', 'hex_areas.R'))
source(here('src', 'accuracy_precision.R'))
source(here('src', 'plotting_funcs.R'))

load(here('data','trig.RData')) #Load nodes and edges between hex area centroids
load(here('data','sample_window.RData')) 

#===============================================================================
##FIGURE 1 -- Map of sampling window, with delaunay triangulation and hex areas 

#Convert edge_info to sf LINESTRING geometry
edge_info_sf <- edge_info %>%
  # Keep only the filtered transitions you already created
  mutate(
    geometry = purrr::pmap(
      list(region1_x, region1_y, region2_x, region2_y),
      ~ st_linestring(matrix(c(..1, ..2, ..3, ..4), ncol = 2, byrow = TRUE))
    )
  ) %>%
  st_as_sf(crs = 3035)  #projected CRS

#Plot
y <- ggplot() +
      geom_sf(data = st_buffer(sampling_win_proj, 40000), fill=NA, color = "grey50", linewidth=1) +
      geom_sf(data = hex_area_win_proj) +
      geom_sf_text(data = hex_area_win_proj, aes(label = area_ID), size=7, alpha=0.8) +
      geom_sf(data = hex_area_win_proj$area_center, size=3, alpha=0.6, color = "purple") +
      geom_sf(data = edge_info_sf, colour='purple', size=0.8, alpha=0.3) +
      labs(x = "Longitude", y = "Latitude") +
      theme(panel.background = element_rect(fill = "lightblue",
                                            colour = "lightblue",
                                            size = 0.5,
                                            linetype = "solid"),
            legend.position = "none",
            axis.text = element_text(size=13),
            text = element_text(size=17))

#Output
pdf(file=here('output','figures','sample_win_hex.pdf'), width=15, height=8)
grid.arrange(y, ncol=1, padding=0)
dev.off()

#===============================================================================
#####PAPER SECTION: Wave of Advance Simulation
#===============================================================================
###SIMULATION 1: TACTICAL WOA SIMULATION (with calendar dates)

##Load data
load(here('data','tactical_sim_woa_noerror.RData')) #simulated data
load(here('output', 'Womblemodel_tactical_woa_noerrors.RData')) #inferred data

#Extract values
sites <- sim_dataset$sites
sites_sf <- sim_dataset$sites_sf
siteInfo <- sim_dataset$siteInfo
sim_df <- sim_dataset$sim_df
constants <- sim_dataset$constants
sampling_win_proj <- sim_dataset$sampling_win_proj
hex_area_win_proj <-sim_dataset$hex_area_win_proj

#------------
#Hex areas with and without out sites
Hex_with_sites <- unique(siteInfo$area_id)
Hex_without_sites <- which(rep(1:81) %!in% Hex_with_sites)

#------------
#Calculate average accuracy and precision
ci_95 = credible_interval(out_womble_model, 0.95)
sim_a <- constants$true_a
sim1_model_accuracy = accuracy(sim_a, ci_95)
sim1_model_precision = precision(sim_a, ci_95)

#Calculate accuracy and precision in each area 
sim1_precision_in_hex <- precision_in_each_area(sim_a, ci_95)
sim1_accuracy_in_hex <- accuracy_in_each_area(sim_a, ci_95)

#-------------------------------------------------------------------------------
##FIGURE 2A -- Map of sites, simulated arrival times and sampling window
#Plot
site_map <- ggplot(data = hex_area_win_proj) +
    geom_sf(data = sampling_win_proj, color = "grey50") +  # sampling window border
    geom_sf(aes(fill = constants$true_a)) +
    scale_fill_viridis_c(name = "Arrival time (in BP)    ", option="F", direction=-1,  limits = c(4500, 6800)) +
    guides(fill = guide_colorbar(direction = "horizontal", barwidth = 13)) + #horizontal legend
    geom_sf(data = sites_sf, size = 3, alpha = 0.8) +  # sites
    theme(
      panel.background = element_rect(fill = "lightblue", colour = "lightblue"),
      plot.margin = margin(t = 0.5 * 10, unit = "mm"), #extra space above panel for legend
      legend.position = c(0.999, 1.05),
      legend.justification = c(1, 1), 
      legend.title.position = "left",
      legend.text = element_text(size=14),
      legend.title = element_text(size=16, face="bold"),
      legend.margin = margin(t = 5, r = 30, b = 5, l = 20),  # padding around legend
      legend.background = element_rect(colour="grey70", size=0.5, linetype="solid"),
      axis.title = element_blank(),
      axis.text = element_text(size=15))

#Output
pdf(file=here('output','figures','sim1_sites.pdf'), width=10, height=8.5)
grid.arrange(site_map, ncol=1, padding=0)
dev.off()

#------------
##FIGURE 2B -- Map of Inferred arrival times
out.comb.tac_icar.model  <- do.call(rbind, out_womble_model)
post.a.model.i  <- out.comb.tac_icar.model[,paste0('a[',1:81,']')]  %>% round() 
med.model.i  <- apply(post.a.model.i, 2, median) #Extract arrival times for tactical icar model

median_hex_dates_mod.i <- hex_area_win_proj %>% 
  filter(area_ID %in% 1:81) %>% 
  mutate(median_date = med.model.i,
         contains_sites = as.factor(case_when(area_ID %in% Hex_with_sites ~ 1, area_ID %in% Hex_without_sites ~ 0))) 

#Plot
modi <- ggplot(data = median_hex_dates_mod.i) +
  geom_sf(data = st_buffer(sampling_win_proj, 40000), color = "grey50") + #sampling window with coastal buffer
  geom_sf(aes(fill = median_date)) + #hex grid #alpha=contains_sites
  scale_fill_viridis_c(name = "Arrival time (in BP)    ", option="F", direction=-1,  limits = c(4500, 6800)) +
  guides(fill = guide_colorbar(direction = "horizontal", barwidth = 13)) + #horizontal legend
  theme(
    panel.background = element_rect(fill = "lightblue", colour = "lightblue"),
    plot.margin = margin(t = 0.5 * 10, unit = "mm"), #extra space above panel for legend
    legend.position = c(0.999, 1.05),
    legend.justification = c(1, 1), 
    legend.title.position = "left",
    legend.text = element_text(size=14),
    legend.title = element_text(size=16, face="bold"),
    legend.margin = margin(t = 5, r = 30, b = 5, l = 20),  # padding around legend
    legend.background = element_rect(colour="grey70", size=0.5, linetype="solid"),
    axis.title = element_blank(),
    axis.text = element_text(size=15))

#Output
pdf(file=here('output','figures','sim1_arrivaltime.pdf'), width=10, height=8.5)
grid.arrange(modi, ncol=1, padding=unit(0,"mm"), clip=FALSE)
dev.off()

#------------
##FIGURE 2C -- Map of difference between simulated and median inferred arrival times

#Extract 2.5% and 97.5% quantiles for the inferred model
ci_lower <- apply(post.a.model.i, 2, quantile, probs = 0.025)
ci_upper <- apply(post.a.model.i, 2, quantile, probs = 0.975)

#Compute difference between simulated (true_a) and inferred median
median_diff <- median_hex_dates_mod.i %>%
  mutate(diffence = constants$true_a - med.model.i,
         outside_CI = constants$true_a < ci_lower | constants$true_a > ci_upper)

#Plot
fig2c <- ggplot(data = median_diff) +
  geom_sf(data = st_buffer(sampling_win_proj, 40000), color = "grey50") + #sampling window with coastal buffer
  geom_sf_pattern(
    aes(fill = diffence, pattern = outside_CI),
    #color = NA,  # no hex borders
    pattern_color = "grey60",
    pattern_angle = 60,
    pattern_density = 0.1,
    pattern_spacing = 0.02) +
  scale_fill_gradient2(name = "Arrival time difference",low = "royalblue1",mid = "white",high = "coral2",midpoint = 0,limits = c(-300,300),oob = scales::squish ) +
  scale_pattern_manual(values = c(`TRUE` = "stripe", `FALSE` = "none")) +
  guides(fill = guide_colorbar(direction = "horizontal", barwidth = 13),pattern = "none") +
  theme(
    panel.background = element_rect(fill = "lightblue", colour = "lightblue"),
    plot.margin = margin(t = 0.5 * 10, unit = "mm"), #extra space above panel for legend
    legend.position = c(0.999, 1.05),
    legend.justification = c(1, 1), 
    legend.title.position = "left",
    legend.text = element_text(size=14),
    legend.title = element_text(size=16, face="bold"),
    legend.margin = margin(t = 5, r = 30, b = 5, l = 20),  # padding around legend
    legend.background = element_rect(colour="grey70", size=0.5, linetype="solid"),
    axis.title = element_blank(),
    axis.text = element_text(size=15))

#Output
pdf(file=here('output','figures','sim1_arrivaltime_difference.pdf'), width=10, height=8.5)
grid.arrange(fig2c, ncol=1, padding=unit(0,"mm"), clip=FALSE)
dev.off()

#===============================================================================
#####PAPER SECTION: Comparing the ICAR model to a Phase model
#===============================================================================
###SIMULATION 2: TACTICAL ICAR SIMULATION (with calibrated radiocarbon dates and one covariate)

##Load data
load(here('data','tactical_sim_icar.RData')) #simulated data
load(here('output', 'Womblemodel_tactical_icar.RData')) #inferred data

#Extract values
sites <- sim_dataset$sites
sites_sf <- sim_dataset$sites_sf
siteInfo <- sim_dataset$siteInfo
sim_df <- sim_dataset$sim_df
constants <- sim_dataset$constants
sampling_win_proj <- sim_dataset$sampling_win_proj
hex_area_win_proj <-sim_dataset$hex_area_win_proj

#------------
#Hex areas with and without out sites
Hex_with_sites <- unique(siteInfo$area_id)
Hex_without_sites <- which(rep(1:81) %!in% Hex_with_sites)

#------------
#Calculate average accuracy and precision
ci_95 = credible_interval(out_womble_model, 0.95)
sim_a <- constants$true_a
sim2_model_accuracy = accuracy(sim_a, ci_95)
sim2_model_precision = precision(sim_a, ci_95)

#-------------------------------------------------------------------------------
##FIGURE 3A -- Map of sites, simulated arrival times and sampling window
#Plot
site_map <- ggplot(data = hex_area_win_proj) +
  geom_sf(data = sampling_win_proj, color = "grey50") +  # sampling window border
  geom_sf(aes(fill = constants$true_a)) +
  scale_fill_viridis_c(name = "Arrival time (in BP)    ", option="F", direction=-1,  limits = c(790, 2400)) +
  guides(fill = guide_colorbar(direction = "horizontal", barwidth = 10)) + #horizontal legend
  geom_sf(data = sites_sf, size = 3, alpha = 0.8) +  # sites
  theme(
    panel.background = element_rect(fill = "lightblue", colour = "lightblue"),
    plot.margin = margin(t = 0.5 * 10, unit = "mm"), #extra space above panel for legend
    legend.position = c(0.999, 1.05),
    legend.justification = c(1, 1), 
    legend.title.position = "left",
    legend.text = element_text(size=14),
    legend.title = element_text(size=16, face="bold"),
    legend.margin = margin(t = 5, r = 30, b = 5, l = 20),  # padding around legend
    legend.background = element_rect(colour="grey70", size=0.5, linetype="solid"),
    axis.title = element_blank(),
    axis.text = element_text(size=15))

#Output
pdf(file=here('output','figures','sim2_sites2.pdf'), width=10, height=8.5)
grid.arrange(site_map, ncol=1, padding=0)
dev.off()

#------------
##FIGURE 3B -- Map of Inferred arrival times

out.comb.tac_icar.model  <- do.call(rbind, out_womble_model)
post.a.model.i  <- out.comb.tac_icar.model[,paste0('a[',1:81,']')]  %>% round() 
med.model.i  <- apply(post.a.model.i, 2, median) #Extract arrival times for tactical icar model

median_hex_dates_mod.i <- hex_area_win_proj %>% 
  filter(area_ID %in% 1:81) %>% 
  mutate(median_date = med.model.i,
         contains_sites = as.factor(case_when(area_ID %in% Hex_with_sites ~ 1, area_ID %in% Hex_without_sites ~ 0))) 

#Plot
modi <- ggplot(data = median_hex_dates_mod.i) +
  geom_sf(data = st_buffer(sampling_win_proj, 40000), color = "grey50") + #sampling window with coastal buffer
  geom_sf(aes(fill = median_date)) + #hex grid #alpha=contains_sites
  scale_fill_viridis_c(name = "Arrival time (in BP)    ", option="F", direction=-1,  limits = c(790, 2400)) +
  guides(fill = guide_colorbar(direction = "horizontal", barwidth = 10)) + #horizontal legend
  theme(
    panel.background = element_rect(fill = "lightblue", colour = "lightblue"),
    plot.margin = margin(t = 0.5 * 10, unit = "mm"), #extra space above panel for legend
    legend.position = c(0.999, 1.05),
    legend.justification = c(1, 1), 
    legend.title.position = "left",
    legend.text = element_text(size=14),
    legend.title = element_text(size=16, face="bold"),
    legend.margin = margin(t = 5, r = 30, b = 5, l = 20),  # padding around legend
    legend.background = element_rect(colour="grey70", size=0.5, linetype="solid"),
    axis.title = element_blank(),
    axis.text = element_text(size=15))

#Output
pdf(file=here('output','figures','sim2_arrivaltime.pdf'), width=10, height=8.5)
grid.arrange(modi, ncol=1, padding=0)
dev.off()


#------------
##FIGURE 3C -- Map of difference between simulated and median inferred arrival times

#Extract 2.5% and 97.5% quantiles for the inferred model
ci_lower <- apply(post.a.model.i, 2, quantile, probs = 0.025)
ci_upper <- apply(post.a.model.i, 2, quantile, probs = 0.975)

#Compute difference between simulated (true_a) and inferred median
median_diff <- median_hex_dates_mod.i %>%
  mutate(diffence = constants$true_a - med.model.i,
         outside_CI = constants$true_a < ci_lower | constants$true_a > ci_upper)

#Plot
fig3c <- ggplot(data = median_diff) +
  geom_sf(data = st_buffer(sampling_win_proj, 40000), color = "grey50") + #sampling window with coastal buffer
  geom_sf_pattern(
    aes(fill = diffence, pattern = outside_CI),
    #color = NA,  # no hex borders
    pattern_color = "grey60",
    pattern_angle = 60,
    pattern_density = 0.1,
    pattern_spacing = 0.02) +
  scale_fill_gradient2(name = "Arrival time difference",low = "royalblue1",mid = "white",high = "coral2",midpoint = 0,limits = c(-300,300),oob = scales::squish ) +
  scale_pattern_manual(values = c(`TRUE` = "stripe", `FALSE` = "none")) +
  guides(fill = guide_colorbar(direction = "horizontal", barwidth = 13),pattern = "none") +
  theme(
    panel.background = element_rect(fill = "lightblue", colour = "lightblue"),
    plot.margin = margin(t = 0.5 * 15, unit = "mm"), #extra space above panel for legend
    legend.position = c(0.999, 1.05),
    legend.justification = c(1, 1), 
    legend.title.position = "left",
    legend.text = element_text(size=14),
    legend.title = element_text(size=16, face="bold"),
    legend.margin = margin(t = 5, r = 30, b = 5, l = 20),  # padding around legend
    legend.background = element_rect(colour="grey70", size=0.5, linetype="solid"),
    axis.title = element_blank(),
    axis.text = element_text(size=15))

#Output
pdf(file=here('output','figures','sim2_arrivaltime_difference.pdf'), width=10, height=8.5)
grid.arrange(fig3c, ncol=1, padding=unit(0,"mm"), clip=FALSE)
dev.off()

#------------
##FIGURE 3D -- Map of Wombling Boundaries (highlighting important boundaries)

#With the tactical simulation data from hierarchical wombling model
post.model.tac_womble_nab  <- out.comb.tac_icar.model[,paste0('nabla[',1:208,']')]  %>% round()

#Extract differences in arrival times for tactical wombling model
med.model.tac_womble_nab  <- apply(post.model.tac_womble_nab, 2, median)

#Extract proportion of MCMC sample differences which are significant over a specified time difference
prop_model_tac_womble_nab  <- data.frame(x = 1:208,
                                         y = sapply(as.data.frame(post.model.tac_womble_nab), 
                                                    prop_gthan_threshold, 
                                                    threshold = 500))
#Add info to edges dataframe
edge_info.i <- edge_info %>%
  mutate(mean_gradient = med.model.tac_womble_nab, #50% quantile
         prob_BLV = prop_model_tac_womble_nab$y, #% of distribution > specified threshold
         boundary = mapply(function(a, b) {intersection <- st_intersection(hex_area_win_proj$geometry[[a]], hex_area_win_proj$geometry[[b]])
         if (st_is_empty(intersection) || st_is(intersection, "MULTILINESTRING")) return(st_linestring()) else return(intersection)},
         edge_info$region1_id,
         edge_info$region2_id, SIMPLIFY = FALSE)) #shared boundary between two subareas

#Create nodes
nodes <- st_coordinates(hex_area_win_proj$area_center)

#Create boundary segments
boundaries <- st_sf(prob_BLV = edge_info.i$prob_BLV,
                    geometry = st_sfc(edge_info.i$boundary))
st_crs(boundaries) <- 3035  # Set CRS for correct Europe projection

#Plot
womble_plot <- ggplot(data = median_hex_dates_mod.i) +
  geom_sf(data = st_buffer(sampling_win_proj, 40000), fill = "grey80", color = "grey40") + #sampling window with coastal buffer
  geom_sf(aes(alpha=0.01), color = "grey70") + #hex grid 
  geom_sf(data = boundaries, lwd=3.5, aes(alpha=prob_BLV), color = "red") +
  scale_alpha_continuous(name = "Boundary Probability", range = c(0, 1), guide = guide_legend(direction = "horizontal")) +  # Use for continuous alpha values
  ggtitle(expression(zeta == 500 ~ years)) +
  geom_sf_label(aes(label = area_ID), size=7) +
  theme(panel.background = element_rect(fill = "lightblue",
                                        colour = "lightblue",
                                        size = 0.5,
                                        linetype = "solid"),
        plot.margin = margin(t = 0.5 * 5, unit = "mm"), #extra space above panel for legend
        legend.position = c(0.999, 1.05),
        legend.justification = c(1, 1), 
        legend.title.position = "left",
        legend.text = element_text(size=14),
        legend.title = element_text(size=16, face="bold"),
        legend.margin = margin(t = 5, r = 30, b = 5, l = 20),  # padding around legend
        legend.background = element_rect(colour="grey70", size=0.5, linetype="solid"),
        axis.text = element_text(size=15),
        axis.title = element_blank(),
        title = element_text(size=16, face="bold"))

#Output
pdf(file=here('output','figures','sim2_womble.pdf'), width=10, height=8.5)
grid.arrange(womble_plot, ncol=1, padding=0)
dev.off()

#-------------------------------------------------------------------------------
##Comparing ICAR simulation (sim 2) with phasemodel simulation (sim 7)

#Load ICAR data
load(here('data','tactical_sim_icar.RData')) #simulated data
load(here('output', 'Womblemodel_tactical_icar.RData')) #inferred data

#Extract values
sites_sf <- sim_dataset$sites_sf
siteInfo <- sim_dataset$siteInfo
sim_df <- sim_dataset$sim_df
constants <- sim_dataset$constants
sampling_win_proj <- sim_dataset$sampling_win_proj
hex_area_win_proj <-sim_dataset$hex_area_win_proj
tmp.a.sim2 = extract(out_womble_model)
out_womble_model2 <- out_womble_model
sim2_a <- constants$true_a
sites_sim2 <- sim_dataset$sites

#Load Phasemodel data
load(here('output', 'Womblemodel_tactical_phasemodel.RData')) #inferred data

tmp.a.sim7 = extract(out_womble_model)
out_womble_model7 <- out_womble_model
sim7_a <- constants$true_a
sites_sim7 <- sites

#------------
#Hex areas with and without out sites
Hex_with_sites <- unique(siteInfo$area_id)
Hex_without_sites <- which(rep(1:81) %!in% Hex_with_sites)

#-------------------------------------------------------------------------------
##FIGURE 4 -- Plot and compare posteriors of arrival times for both models

pdf(file=here('output','figures','sim2_vs_sim7_posteriors.pdf'), width=10, height=6, pointsize=4)
par(mar = c(5, 5, 4, 2))   #pad space around plot
plot(NULL, xlim=c(3250, 500), ylim=c(0.5, 17), xlab=paste('Arrival time (BP),', TeX('$a_k$')), ylab=paste('Area,', TeX('$k$')), cex.lab = 2, axes=F)

iseq.a = seq(1,by=1,length.out=16)
abline(h=seq(1,by=1,length.out=16), col='lightgrey')

counter <- 1 #indexing counter
for (i in c(1,14,11,15,24,31,32,45,50,53,59,61,65,68,79,81))
{
  #Plot bars from sim 5 and sim 6 in area i
  post.bar(tmp.a.sim2[,i], i=iseq.a[counter], h=0.5, a= sim2_a[[i]], barcolours=barcolours1)
  post.bar(tmp.a.sim7[,i], i=iseq.a[counter]+0.3, h=0.5, a= sim7_a[[i]], barcolours=barcolours2)
  counter <- counter + 1
}

axis(2, at=iseq.a, labels = paste0(c(1,14,11,15,24,31,32,45,50,53,59,61,65,68,79,81)), las=2, cex.axis=1.7)
axis(1, at = BCADtoBP(c(-1300, -1100, -900, -700, -500, -300, -100, 100, 300, 500, 700, 900, 1100, 1300)), labels=c('1300BC','1100BC','900BC','700BC','500BC','300BC', '100BC', '100AD', '300AD', '500AD', '700AD', '900AD', '1100AD','1300AD'), tck=-0.01, cex.axis=1.7)
axis(3, at = seq(3200, 500, -200), labels=paste0(seq(3200, 500, -200),'BP'), tck=-0.01, cex.axis=1.7)
axis(1, at = BCADtoBP(c(-1200,-1000,-800, -600, -400, -200, 1, 200, 400, 600, 800, 1000, 1200, 1400)), labels=NA, tck=-0.01) #Minor tick marks
axis(3, at = seq(3200, 500, -50), labels=NA, tck=-0.01) #Minor tick marks
box()
#Legend
rect(xleft=3150, xright=2300, ybottom=3, ytop=6.5, border="darkgrey", col="white", lwd=2.5)
post.bar(c(3100,3000,2800,2700,2600,2400,2300), i=5.5, h=0.9, a=2550)
text(x=2650, y=6, "50% HPDI", cex=2)
text(x=2950, y=6,"95% HPDI", cex=2)
text(x=2700, y=5, "Median Posterior", cex=2)
text(x=2550, y=4.5, "Simulated value", cex=2)
segments(y0=5.3, y1=4.6, x0=2550, col="darkgrey", lwd=2)
text(x=2750, y=4, "ICAR Model", col="black", cex=2)
segments(x0=3050, x1=3000, y0=4, col="dodgerblue", lwd=4)
text(x=2750, y=3.5, "Phasemodel", col="black", cex=2)
segments(x0=3050, x1=3000, y0=3.5, col="orchid", lwd=4)
theme(legend.position = "none")
dev.off()


#-------------------------------------------------------------------------------
##FIGURE 5 -- Sample size per hexagon vs. precision for both icar and phasemodels (sim 2 vs sim 7)

#Number of sites in area
sites2_in_areas_summarise <- sites_sim2 %>% 
  group_by(area_id) %>% 
  summarize(n_sites = n_distinct(site_id), .groups="drop") %>%
  complete(area_id = full_seq(1:max(area_id), 1), fill = list(n_sites = 0))

#CHECK: sites7_in_areas_summarise$n_sites==sites2_in_areas_summarise$n_sites
#  sites7_in_areas_summarise <- sites_sim7 %>% 
#  group_by(area_id) %>% 
#  summarize(n_sites = n_distinct(site_id), .groups="drop") %>%
#  complete(area_id = full_seq(min(area_id):max(area_id), 1), fill = list(n_sites = 0))

#Precision in each area
ci_95_sim2 = credible_interval(out_womble_model2, 0.95)
precision2_in_hex <- precision_in_each_area(out_womble_model2, ci_95_sim2)

ci_95_sim7 = credible_interval(out_womble_model7, 0.95)
precision7_in_hex <- precision_in_each_area(out_womble_model7, ci_95_sim7)

#Precision
df_precision <- data.frame(n_sites = sites2_in_areas_summarise$n_sites,
                           precision2 = precision2_in_hex,
                           precision7 = precision7_in_hex)

df_long <- df_precision %>% 
  rename("ICAR Model" = precision2, "Phasemodel" = precision7) %>% 
  pivot_longer(cols = c("ICAR Model", "Phasemodel"), 
               names_to = "simulation", 
               values_to = "precision") %>% 
  mutate(simulation = factor(simulation))

#Plot
p0 <- ggplot(df_long, aes(x = n_sites, y = precision, color = simulation, shape = simulation)) +
  geom_point(size = 3) +
  labs(x = "Number of Sites",
       y = "Precision (in years)",
       color="Model",
       shape = "Model") +
  scale_shape_manual(values = c("ICAR Model" = 16, "Phasemodel" = 17)) +
  scale_x_continuous(breaks = round(seq(0, 35, by = 5),1)) +
  scale_y_continuous(breaks = round(seq(0, 3000, by = 500),1), limits =c(0,3000)) +
  theme_minimal() +
  theme(axis.text = element_text(size=14),
        text = element_text(size=18),
        legend.position = c(0.96, 0.96),
        legend.justification = c(1, 1), 
        legend.title.position = "top",
        legend.text = element_text(size=14),
        legend.title = element_text(size=16, face="bold"),
        legend.margin = margin(t = 5, r = 30, b = 5, l = 20),  # padding around legend
        legend.background = element_rect(colour="grey70", size=0.5, linetype="solid"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=2))

#Output
pdf(file=here('output','figures','precision_vs_nsites_sim27.pdf'), width=10, height=8)
grid.arrange(p0, ncol=1, padding=0)
dev.off()

#===============================================================================
#####PAPER SECTION: Robustness to variation in sample size and sampling intensity
#===============================================================================
##SIMULATION 8: Examining the effect of number of sites

load(here('output', 'Womblemodel_tactical_icar_sitenumber.RData'))

#Calculate accuracy and precision
accuracy_models <- list()
precision_models <- list()
precision_range <- list()

for (i in seq(1:length(output_varysitenumber_models))){
  model <- output_varysitenumber_models[[i]][[1]]
  
  ci_95 = credible_interval(model, 0.95)
  sim_a <- output_varysitenumber_models[[i]][[4]]$true_a
  model_acc = accuracy(sim_a, ci_95)
  model_precis = precision(sim_a, ci_95)
  model_precis_each_hex_range = sd(precision_in_each_area(sim_a, ci_95))
   
  accuracy_models[[i]] <- model_acc
  precision_models[[i]] <- model_precis
  precision_range[[i]] <- model_precis_each_hex_range
} 

model_metrics <- data.frame(n_sites = seq(200,2000,200),
                            accuracy = unlist(accuracy_models),
                            precision = unlist(precision_models),
                            precison_sd = unlist(precision_range))

#-------------------------------------------------------------------------------
##FIGURE 6A -- Accuracy vs.number of sites
p1 <- ggplot(model_metrics, aes(x = n_sites, y = accuracy)) +
  geom_point(size = 3.5) +
  scale_x_continuous(breaks = round(seq(min(model_metrics$n_sites), max(model_metrics$n_sites), by = 200),1)) +
  scale_y_continuous(breaks = round(seq(0, 1, by = 0.1),1), limits =c(0,1)) +
  labs(x = "Number of Sites",
       y = "Accuracy") +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=2),
        axis.text = element_text(size=15),
        text = element_text(size=19)) 

#Output
pdf(file=here('output','figures','accuracy_vs_nsites.pdf'), width=10, height=8)
grid.arrange(p1, ncol=1, padding=0)
dev.off()

#---------------
##FIGURE 6B -- Precision vs. number of sites
p2 <- ggplot(model_metrics, aes(x = n_sites, y = precision)) +
  geom_errorbar(aes(ymin = pmax(precision - precison_sd,0),
                    ymax = precision + precison_sd),
                width = 50, 
                linewidth = 0.8,
                color = "darkgrey") +
  geom_point(size = 3.5) +
  scale_x_continuous(breaks = round(seq(min(model_metrics$n_sites), max(model_metrics$n_sites), by = 200),1)) +
  scale_y_continuous(breaks = round(seq(0, max(model_metrics$precision)+400, by = 200),1), limits =c(0,2600)) +
  labs(x = "Number of Sites",
       y = "Precision (in years)") +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=2),
        axis.text = element_text(size=15),
        text = element_text(size=19)) 

#Output
pdf(file=here('output','figures','precision_vs_nsites.pdf'), width=10, height=8)
grid.arrange(p2, ncol=1, padding=0)
dev.off()

#===============================================================================
##SIMULATION 9: Examining the effect of sampling intensity

load(here('output', 'Womblemodel_tactical_icar_clustering.RData'))

#Calculate accuracy and precision
accuracy_models <- list()
precision_models <- list()
precision_range <- list()

for (i in seq(1:length(output_varyclusterdeg_models))){
  model <- output_varyclusterdeg_models[[i]][[1]]
  
  ci_95 = credible_interval(model, 0.95)
  sim_a <- output_varyclusterdeg_models[[i]][[4]]$true_a
  model_acc = accuracy(sim_a, ci_95)
  model_precis = precision(sim_a, ci_95)
  model_precis_each_hex_range = sd(precision_in_each_area(sim_a, ci_95))
  
  accuracy_models[[i]] <- model_acc
  precision_models[[i]] <- model_precis
  precision_range[[i]] <- model_precis_each_hex_range
}

model_metrics <- data.frame(cluster_deg = seq(0,1,0.1),
                            accuracy = unlist(accuracy_models),
                            precision = unlist(precision_models),
                            precison_sd = unlist(precision_range))
#-------------------------------------------------------------------------------
##FIGURE 7 -- Examine clustering of sites
source(here('src', 'sim_data.R'))

plots <- list()
for (c in c(0,0.8,1)){
  #Generate simulated data set
  sim_dataset <- sim_data(with_calibration = FALSE, seed=123, k = c, n_sites = 800, n_dates = 2400)
  
#Plot
p <- ggplot(data = hex_area_win_proj) +
  geom_sf(data = sampling_win_proj, color = "grey50") +  # sampling window border
  geom_sf() +
  geom_sf(data = sim_dataset$sites_sf, size = 2, alpha = 0.5) +  # sites  
  labs(title=paste0("c = ", c)) +
  theme(
    panel.background = element_rect(fill = "lightblue", colour = "lightblue"),
    legend.position = "bottom",
    axis.title = element_blank(),
    axis.text = element_text(size=13))

plots <- list.append(plots, p)}

#Output
pdf(file=here('output','figures','clustering_sites.pdf'), width=10, height=3)
wrap_plots(plots, nrow=1, padding=0)
dev.off()

#-------------------------------------------------------------------------------
##FIGURE 8A -- Accuracy vs. sampling intensity
p1 <- ggplot(model_metrics, aes(x = cluster_deg, y = accuracy)) +
  geom_point(size = 3.5) +
  scale_x_continuous(breaks = round(seq(min(model_metrics$cluster_deg), max(model_metrics$cluster_deg), by = 0.1),1)) +
  scale_y_continuous(breaks = round(seq(0, 1, by = 0.1),1), limits =c(0,1)) +
  labs(x = "Sampling Intensity",
       y = "Accuracy") +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=2),
        axis.text = element_text(size=15),
        text = element_text(size=19)) 

#Output
pdf(file=here('output','figures','accuracy_vs_sampintesity.pdf'), width=10, height=8)
grid.arrange(p1, ncol=1, padding=0)
dev.off()

#---------------
##FIGURE 8B -- Precision vs. sampling intensity
p2 <- ggplot(model_metrics, aes(x = cluster_deg, y = precision)) +
  geom_errorbar(aes(ymin = pmax(precision - precison_sd,0),
                    ymax = precision + precison_sd),
                width = 0.03,   # adjust width of error bars
                linewidth = 0.8,
                color = "darkgrey") +
  geom_point(size = 3.5) +
  scale_x_continuous(breaks = round(seq(min(model_metrics$cluster_deg), max(model_metrics$cluster_deg), by = 0.1),1)) +
  scale_y_continuous(breaks = round(seq(0, max(model_metrics$precision)+1000, by = 200),1), limits =c(0,1600)) +
  labs(x = "Sampling Intensity",
       y = "Precision (in years)") +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=2),
        axis.text = element_text(size=15),
        text = element_text(size=19)) 


#Output
pdf(file=here('output','figures','precision_vs_sampintensity.pdf'), width=10, height=8)
grid.arrange(p2, ncol=1, padding=0)
dev.off()