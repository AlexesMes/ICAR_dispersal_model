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

#-------------------------------------------------------------------------------
##SUPPLEMENTARY Diagnostics: Traceplots and metric table 

#Create list to store recorded plots
plots <- vector("list", 81)

#Output
pdf(file = here("output", "supplementary_figures", "traceplots_woa.pdf"), width = 10, height = 15, onefile=TRUE)
par(mfrow = c(11,8), mar = c(2,0,2,0), oma = c(0, 0, 0, 0), mgp = c(1, 0.2, 0))
for (i in 1:81) {
  param_name <- paste0("a[", i, "]")
  plot.new()
  traceplot(out_womble_model[, param_name], main = TeX(paste0("$a[", i, "]$")), smooth = TRUE)
  plots[[i]] <- recordPlot()
}
dev.off()

#Diagnostics table
out.comb.tac_icar.model  <- do.call(rbind, out_womble_model)
post.a.model.i  <- out.comb.tac_icar.model[,paste0('a[',1:81,']')]  %>% round() 
med.model.i  <- apply(post.a.model.i, 2, median)

diagnostic_df <- data.frame(median_posterior = paste(med.model.i, "BP"),
                            HPDI95_low = paste(round(ci_95[1,]), "BP"),
                            HPDI95_high = paste(round(ci_95[2,]), "BP"),
                            rhat = round(rhat_womble_model$psrf[1:81,1],2),
                            ESS = round(ess_womble_model[1:81]))

write.csv(diagnostic_df,file=here('output','tables','diagnostics_woa.csv'), row.names = TRUE)

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

#-------------------------------------------------------------------------------
##SUPPLEMENTARY FIGURE -- Posterior distributions of arrival times
#Select parameters a and b (i.e. start and end date of occupation in the region)
sim_a <- constants$true_a
sim_b <- constants$true_b

#Plot
pdf(file=here('output',"supplementary_figures",'sim1_posteriors.pdf'), width=10, height=15, pointsize=4)
par(mar = c(5, 5, 4, 2))   #pad space around plot
plot(NULL, xlim=c(6800, 4600), ylim=c(3, 79), xlab=paste('Arrival time (BP),', TeX('$a_k$')), ylab=paste('Area,', TeX('$k$')), cex.lab = 2, axes=F)
tmp.a = extract(out_womble_model)
iseq.a = seq(1,by=1,length.out=81)
abline(h=seq(1,by=1,length.out=81), col='lightgrey')

counter <- 1 #indexing counter
for (i in c(1:81)) #all relevant hex areas
{
  #Plot bar in area i
  post.bar(tmp.a[,i], i=iseq.a[counter], h=0.5, a= sim_a[[i]])
  counter <- counter + 1
}

axis(2, at=iseq.a, labels = paste0(c(1:81)), las=2, cex.axis=1.7)
axis(1, at = BCADtoBP(c(-4900, -4700, -4500, -4300, -4100, -3900, -3700, -3500, -3300, -3100, -2900, -2700)), labels=c('4900BC','4700BC','4500BC','4300BC', '4100BC', '3900BC', '3700BC', '3500BC', '3300BC', '3100BC', '2900BC','2700BC'), tck=-0.01, cex.axis=1.7)
axis(3, at = seq(6800, 4600, -200), labels=paste0(seq(6800, 4600, -200),'BP'), tck=-0.01, cex.axis=1.7)
axis(1, at = BCADtoBP(c(-4800, -4600, -4400, -4200, -4000, -3800, -3600, -3400, -3200, -3000, -2800)), labels=NA, tck=-0.01) #Minor tick marks
axis(3, at = seq(6800, 4600, -50), labels=NA, tck=-0.01) #Minor tick marks
box()
#Legend
post.bar(c(6900,6800,6600,6500,6400,6200,6100), i=77, h=0.9, a=6750)
text(x=6550, y=78, "50% HPDI", cex=1.5)
text(x=6300, y=78,"95% HPDI", cex=1.5)
text(x=6450, y=76, "Median Posterior", cex=1.5)
text(x=6700, y=76, "Simulated value", cex=1.5)
rect(xleft=6850, xright=6150, ybottom=75, ytop=79, border="darkgrey", col=NA, lwd=2)
theme(legend.position = "none")
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
##SUPPLEMENTARY Diagnostics: Traceplots and metric table 

#Create list to store recorded plots
plots <- vector("list", 81)

#Output
pdf(file = here("output", "supplementary_figures", "traceplots_icar.pdf"), width = 10, height = 15, onefile=TRUE)
par(mfrow = c(11,8), mar = c(2,0,2,0), oma = c(0, 0, 0, 0), mgp = c(1, 0.2, 0))
for (i in 1:81) {
  param_name <- paste0("a[", i, "]")
  plot.new()
  traceplot(out_womble_model[, param_name], main = TeX(paste0("$a[", i, "]$")), smooth = TRUE)
  plots[[i]] <- recordPlot()
}
dev.off()

#Diagnostics table
out.comb.tac_icar.model  <- do.call(rbind, out_womble_model)
post.a.model.i  <- out.comb.tac_icar.model[,paste0('a[',1:81,']')]  %>% round() 
med.model.i  <- apply(post.a.model.i, 2, median)

diagnostic_df <- data.frame(median_posterior = paste(med.model.i, "BP"),
                            HPDI95_low = paste(round(ci_95[1,]), "BP"),
                            HPDI95_high = paste(round(ci_95[2,]), "BP"),
                            rhat = round(rhat_womble_model$psrf[1:81,1],2),
                            ESS = round(ess_womble_model[1:81]))

write.csv(diagnostic_df,file=here('output', 'tables','diagnostics_icar.csv'), row.names = TRUE)

#-------------------------------------------------------------------------------
##SUPPLEMENTARY FIGURE -- Map of sites, presence of covariate and sampling window
#Plot
site_map <- ggplot(data = hex_area_win_proj) +
  geom_sf(data = sampling_win_proj, color = "grey50") +  # sampling window border
  geom_sf() +
  geom_sf(aes(fill = env_type)) + # color by combined type
  scale_fill_manual(
    values = c("No forest" = "grey90", "Forest" = "forestgreen")) +
  geom_sf(data = sites_sf, size = 3, alpha = 0.5) +  # sites
  geom_sf_label(aes(label = area_ID), size=7) +                     # area labels
  theme(
    panel.background = element_rect(fill = "lightblue", colour = "lightblue"),
    legend.title = element_blank(),
    legend.position = "top",
    legend.text=element_text(size=18),
    axis.title = element_blank(),
    axis.text = element_text(size=14))

#Output
pdf(file=here('output','supplementary_figures','sim2_sites.pdf'), width=10, height=8)
grid.arrange(site_map, ncol=1, padding=0)
dev.off()

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
##FIGURE 3C -- Map of Wombling Boundaries (highlighting important boundaries)

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
  geom_sf(aes(alpha=0.01), color = "grey60") + scale_alpha(range = c(0, 1)) + #hex grid 
  geom_sf(data = boundaries, lwd=3, aes(alpha=prob_BLV), color = "red") +
  geom_sf(data = hex_area_win_proj$area_center, size=2, alpha=1, color = "grey40") + #hex-centers
  scale_alpha_continuous(range = c(0, 1)) +  # Use for continuous alpha values
  ggtitle(paste0('c = 500', ' years')) +
  geom_sf_label(aes(label = area_ID), size=7) +
  theme(panel.background = element_rect(fill = "lightblue",
                                        colour = "lightblue",
                                        size = 0.5,
                                        linetype = "solid"),
        legend.position = "none",
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
post.bar(c(3100,3000,2800,2700,2600,2400,2300), i=5.5, h=0.9, a=2550)
text(x=2650, y=6, "50% HPDI", cex=2)
text(x=2950, y=6,"95% HPDI", cex=2)
text(x=2700, y=5, "Median Posterior", cex=2)
text(x=2550, y=4.5, "Simulated value", cex=2)
segments(y0=5.3, y1=4.6, x0=2550, col="darkgrey", lwd=2)
rect(xleft=3150, xright=2300, ybottom=3, ytop=6.5, border="darkgrey", col=NA, lwd=2.5)
text(x=2750, y=4, "ICAR Model", col="black", cex=2)
segments(x0=3050, x1=3000, y0=4, col="dodgerblue", lwd=4)
text(x=2750, y=3.5, "Phasemodel", col="black", cex=2)
segments(x0=3050, x1=3000, y0=3.5, col="orchid", lwd=4)
theme(legend.position = "none")
dev.off()

#-------------------------------------------------------------------------------
##SUPPLEMENTARY FIGURE -- Posterior distributions of arrival times

pdf(file=here('output','supplementary_figures','sim2_vs_sim7_posteriors_even.pdf'), width=10, height=15, pointsize=4)
par(mar = c(5, 5, 4, 2))   #pad space around plot
plot(NULL, xlim=c(3250, 500), ylim=c(1.5, 39.5), xlab=paste('Arrival time (BP),', TeX('$a_k$')), ylab=paste('Area,', TeX('$k$')), cex.lab = 2, axes=F)

iseq.a = seq(1,by=1,length.out=40)
abline(h=seq(1,by=1,length.out=40), col='lightgrey')

counter <- 1 #indexing counter
for (i in seq(2,80,2)) #all even hex areas #for odd hex areas use: seq(1,81,2)
{
  #Plot bars from sim 2 and sim 7 in area i
  post.bar(tmp.a.sim2[,i], i=iseq.a[counter], h=0.5, a= sim2_a[[i]], barcolours=barcolours1)
  post.bar(tmp.a.sim7[,i], i=iseq.a[counter]+0.3, h=0.5, a= sim7_a[[i]], barcolours=barcolours2)
  counter <- counter + 1
}

axis(2, at=iseq.a, labels = paste0(seq(2,80,2)), las=2, cex.axis=1.7) #for odd hex areas use: seq(1,81,2)
axis(1, at = BCADtoBP(c(-1300, -1100, -900, -700, -500, -300, -100, 100, 300, 500, 700, 900, 1100, 1300)), labels=c('1300BC','1100BC','900BC','700BC','500BC','300BC', '100BC', '100AD', '300AD', '500AD', '700AD', '900AD', '1100AD','1300AD'), tck=-0.01, cex.axis=1.7)
axis(3, at = seq(3200, 500, -200), labels=paste0(seq(3200, 500, -200),'BP'), tck=-0.01, cex.axis=1.7)
axis(1, at = BCADtoBP(c(-1200,-1000,-800, -600, -400, -200, 1, 200, 400, 600, 800, 1000, 1200, 1400)), labels=NA, tck=-0.01) #Minor tick marks
axis(3, at = seq(3200, 500, -50), labels=NA, tck=-0.01) #Minor tick marks
box()
#Legend
post.bar(c(3100,3000,2800,2700,2600,2400,2300), i=10.5, h=0.9, a=2550)
text(x=2600, y=11, "50% HPDI", cex=2)
text(x=2950, y=11,"95% HPDI", cex=2)
text(x=2700, y=10, "Median Posterior", cex=2)
text(x=2550, y=9.5, "Simulated value", cex=2)
segments(y0=10.3, y1=9.6, x0=2550, col="darkgrey", lwd=2)
rect(xleft=3150, xright=2300, ybottom=8, ytop=11.5, border="darkgrey", col=NA, lwd=2.5)
text(x=2750, y=9, "ICAR Model", col="black", cex=2)
segments(x0=3050, x1=3000, y0=9, col="dodgerblue", lwd=4)
text(x=2750, y=8.5, "Phasemodel", col="black", cex=2)
segments(x0=3050, x1=3000, y0=8.5, col="orchid", lwd=4)
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
ci_50_sim2 = credible_interval(out_womble_model2, 0.50)
precision2_in_hex <- precision_in_each_area(out_womble_model2, ci_50_sim2)

ci_50_sim7 = credible_interval(out_womble_model7, 0.50)
precision7_in_hex <- precision_in_each_area(out_womble_model7, ci_50_sim7)

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
  scale_y_continuous(breaks = round(seq(0, 1600, by = 500),1), limits =c(0,1700)) +
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

#-------------------------------------------------------------------------------
##SUPPLEMENTARY FIGURE -- (A) ICAR Precision and accuracy vs. number of sites

#Accuracy (hit/miss) in each area
ci_50_sim2 = credible_interval(out_womble_model2, 0.50)
accuracy2_in_hex <- accuracy_in_each_area(sim2_a, ci_50_sim2)

df_ICAR_precis_acc <- data.frame(n_sites = sites2_in_areas_summarise$n_sites,
                           precision = precision2_in_hex,
                           accuracy = accuracy2_in_hex)

p1 <- ggplot(df_ICAR_precis_acc, aes(x = n_sites, y = precision, color = accuracy)) + 
  geom_point(size = 3) +
  labs(x = "Number of Sites",
    y = "Precision (in years)",
    color = "Accuracy (50% HPDI)") +
  scale_shape_manual(values = 16) +
  scale_color_manual(values = c("FALSE" = "green", "TRUE"= "blue")) +
  scale_x_continuous(breaks = round(seq(0, 35, by = 5),1)) +
  scale_y_continuous(breaks = round(seq(0, 1600, by = 500),1), limits =c(0,1700)) +
  theme_minimal() +
  theme(axis.text = element_text(size=14),
        text = element_text(size=18),
        legend.position = c(0.96, 0.96),
        legend.justification = c(1, 1), 
        legend.title.position = "top",
        legend.text = element_text(size=14),
        legend.title = element_text(size=16, face="bold"),
        legend.margin = margin(t = 5, r = 30, b = 5, l = 20), 
        legend.background = element_rect(colour="grey70", size=0.5, linetype="solid"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=2))

pdf(file=here('output','supplementary_figures','precision_vs_nsites_sim27_icar_acc.pdf'), width=10, height=8)
grid.arrange(p1, ncol=1, padding=0)
dev.off()

#---------------
##SUPPLEMENTARY FIGURE -- (B) Phasemodel Precision and accuracy vs. number of sites
#Accuracy (hit/miss) in each area
ci_50_sim2 = credible_interval(out_womble_model2, 0.50)
accuracy2_in_hex <- accuracy_in_each_area(sim2_a, ci_50_sim2)

ci_50_sim7 = credible_interval(out_womble_model7, 0.50)
accuracy7_in_hex <- accuracy_in_each_area(sim7_a, ci_50_sim7)

df_Phase_precis_acc <- data.frame(n_sites = sites2_in_areas_summarise$n_sites,
                                 precision = precision7_in_hex,
                                 accuracy = accuracy7_in_hex)

p2 <- ggplot(df_Phase_precis_acc, aes(x = n_sites, y = precision, color = accuracy)) + 
  geom_point(size = 3) +
  labs(x = "Number of Sites",
    y = "Precision (in years)",
    color = "Accuracy (50% HPDI)") +
  scale_shape_manual(values = 17) +
  scale_color_manual(values = c("FALSE" = "green", "TRUE"= "blue")) +
  scale_x_continuous(breaks = round(seq(0, 35, by = 5),1)) +
  scale_y_continuous(breaks = round(seq(0, 1600, by = 500),1), limits =c(0,1700)) +
  theme_minimal() +
  theme(axis.text = element_text(size=14),
        text = element_text(size=18),
        legend.position = c(0.96, 0.96),
        legend.justification = c(1, 1), 
        legend.title.position = "top",
        legend.text = element_text(size=14),
        legend.title = element_text(size=16, face="bold"),
        legend.margin = margin(t = 5, r = 30, b = 5, l = 20),
        legend.background = element_rect(colour="grey70", size=0.5, linetype="solid"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=2))

pdf(file=here('output','supplementary_figures','precision_vs_nsites_sim27_phase_acc.pdf'), width=10, height=8)
grid.arrange(p2, ncol=1, padding=0)
dev.off()

#---------------
##SUPPLEMENTARY TABLE -- Average precision per hexagon for icar and phasemodel
df_precision_diff <- df_precision %>% 
  mutate(precision_diff = precision2-precision7) %>% 
  group_by(n_sites) %>% 
  summarise(
    n_areas = n(),
    avg_precision_diff = paste(round(mean(precision_diff, na.rm = TRUE)), "years"),
    sd_precision_diff = paste(round(sd(precision_diff, na.rm = TRUE)),"years"))

write.csv(df_precision_diff,file=here('output','tables','icar_phase_hex_precis.csv'), row.names = FALSE)


#===============================================================================
#####PAPER SECTION: The effect of the 800-400 BC calibration plateau
#===============================================================================
#SIMULATION 3: OUT PLATEAU ICAR SIMULATION (with calibrated radiocarbon dates)

#Load data
load(here('data','tactical_sim_icar_withoutplat.RData')) #simulated data
load(here('output', 'Womblemodel_tactical_withoutplat_errors.RData')) #inferred data

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
sim3_model_accuracy = accuracy(sim_a, ci_95)
sim3_model_precision = precision(sim_a, ci_95)

#-------------------------------------------------------------------------------
##SUPPLEMENTARY Diagnostics: Traceplots and metric table 

#Create list to store recorded plots
plots <- vector("list", 81)

#Output
pdf(file = here("output", "supplementary_figures", "traceplots_outplateau_cal.pdf"), width = 10, height = 15, onefile=TRUE)
par(mfrow = c(11,8), mar = c(2,0,2,0), oma = c(0, 0, 0, 0), mgp = c(1, 0.2, 0))
for (i in 1:81) {
  param_name <- paste0("a[", i, "]")
  plot.new()
  traceplot(out_womble_model[, param_name], main = TeX(paste0("$a[", i, "]$")), smooth = TRUE)
  plots[[i]] <- recordPlot()
}
dev.off()

#Diagnostics table
out.comb.tac_icar.model  <- do.call(rbind, out_womble_model)
post.a.model.i  <- out.comb.tac_icar.model[,paste0('a[',1:81,']')]  %>% round() 
med.model.i  <- apply(post.a.model.i, 2, median)

diagnostic_df <- data.frame(median_posterior = paste(med.model.i, "BP"),
                            HPDI95_low = paste(round(ci_95[1,]), "BP"),
                            HPDI95_high = paste(round(ci_95[2,]), "BP"),
                            rhat = round(rhat_womble_model$psrf[1:81,1],2),
                            ESS = round(ess_womble_model[1:81]))

write.csv(diagnostic_df,file=here('output','tables','diagnostics_outplateau_cal.csv'), row.names = TRUE)

#-------------------------------------------------------------------------------
##SUPPLEMENTARY FIGURE -- Map of sites, presence of covariates and sampling window
#Plot
site_map <- ggplot(data = hex_area_win_proj) +
  geom_sf(data = sampling_win_proj, color = "grey50") +  # sampling window border
  geom_sf() +
  geom_sf(aes(fill = env_type)) + # color by combined type
  scale_fill_manual(
    values = c("Neither" = "grey90", "Forest Only" = "forestgreen", "Water Only" = "skyblue", "Forest & Water" = "hotpink4")) +
  geom_sf(data = sites_sf, size = 3, alpha = 0.5) +  # sites
  geom_sf_label(aes(label = area_ID), size=7) +                     # area labels
  theme(
    panel.background = element_rect(fill = "lightblue", colour = "lightblue"),
    legend.title = element_blank(),
    legend.position = "top",
    legend.text=element_text(size=17),
    axis.title = element_blank(),
    axis.text = element_text(size=13))

#Output
pdf(file=here('output','supplementary_figures','sim3_sites.pdf'), width=10, height=8)
grid.arrange(site_map, ncol=1, padding=0)
dev.off()

#-------------------------------------------------------------------------------
##FIGURE 6A -- Map of sites, simulated arrival times and sampling window
#Plot
site_map <- ggplot(data = hex_area_win_proj) +
  geom_sf(data = sampling_win_proj, color = "grey50") +  # sampling window border
  geom_sf(aes(fill = constants$true_a)) +
  scale_fill_viridis_c(name = "Arrival time (in BP)    ", option="F", direction=-1,  limits = c(790, 2400)) +
  guides(fill = guide_colorbar(direction = "horizontal", barwidth = 10)) + #horizontal legend
  geom_sf(data = sites_sf, size = 3, alpha = 0.5) +  # sites
  theme(
    panel.background = element_rect(fill = "lightblue", colour = "lightblue"),
    plot.margin = margin(t = 0.5 * 10, unit = "mm"), #extra space above panel for legend
    legend.position = c(0.999, 1.05),
    legend.justification = c(1, 1), 
    legend.title.position = "left",
    legend.text = element_text(size=14),
    legend.title = element_text(size=16, face="bold"),
    legend.margin = margin(t = 5, r = 30, b = 5, l = 20),
    legend.background = element_rect(colour="grey70", size=0.5, linetype="solid"),
    axis.title = element_blank(),
    axis.text = element_text(size=15))

#Output
pdf(file=here('output','figures','sim3_sites2.pdf'), width=10, height=8.5)
grid.arrange(site_map, ncol=1, padding=0)
dev.off()

#------------
##FIGURE 6B -- Map of Inferred arrival times
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
pdf(file=here('output','figures','sim3_arrivaltime.pdf'), width=10, height=8.5)
grid.arrange(modi, ncol=1, padding=0)
dev.off()

#------------
##FIGURE 7A -- Map of Wombling Boundaries (highlighting important boundaries)

#With the tactical simulation data from hierarchical wombling model
post.model.tac_womble_nab  <- out.comb.tac_icar.model[,paste0('nabla[',1:208,']')]  %>% round()

#Extract differences in arrival times for tactical wombling model
med.model.tac_womble_nab  <- apply(post.model.tac_womble_nab, 2, median)

#Extract proportion of MCMC sample differences which are significant over a specified time difference
prop_model_tac_womble_nab  <- data.frame(x = 1:208,
                                         y = sapply(as.data.frame(post.model.tac_womble_nab),
                                                    prop_gthan_threshold,
                                                    threshold = 400))

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
                    geometry = st_sfc(edge_info.i$boundary)) #lapply(edge_info.i$boundary[[a]], st_coordinates(a))
st_crs(boundaries) <- 3035  # Set CRS for correct Europe projection

#Plot
pdf(file=here('output','figures','sim3_womble.pdf'), width=10, height=8)
ggplot(data = median_hex_dates_mod.i) +
  geom_sf(data = st_buffer(sampling_win_proj, 40000), fill = "grey80", color = "grey40") + #sampling window with coastal buffer
  geom_sf(aes(alpha=0.01), color = "grey60") + scale_alpha(range = c(0, 1)) + #hex grid
  geom_sf(data = boundaries, lwd=3, aes(alpha=prob_BLV), color = "red") +
  geom_sf(data = hex_area_win_proj$area_center, size=2, alpha=1, color = "grey40") + #hex-centers
  scale_alpha_continuous(range = c(0, 1)) +  # Use for continuous alpha values
  ggtitle(paste0('c = 400', ' years')) +
  geom_sf_label(aes(label = area_ID), size=7) +
  theme(panel.background = element_rect(fill = "lightblue",
                                        colour = "lightblue",
                                        size = 0.5,
                                        linetype = "solid"),
        legend.position = "none",
        axis.text = element_text(size=15),
        axis.title = element_blank(),
        title = element_text(size=16, face="bold"))
dev.off()

#-------------------------------------------------------------------------------
##SUPPLEMENTARY FIGURE -- Posterior distributions of arrival times
#Select parameters a and b (i.e. start and end date of occupation in the region)
sim_a <- constants$true_a
sim_b <- constants$true_b

#Plot
pdf(file=here('output',"supplementary_figures",'sim3_posteriors.pdf'), width=10, height=15, pointsize=4)
par(mar = c(5, 5, 4, 2))   #pad space around plot
plot(NULL, xlim=c(3200, 600), ylim=c(3, 79), xlab=paste('Arrival time (BP),', TeX('$a_k$')), ylab=paste('Area,', TeX('$k$')), cex.lab = 2, axes=F)
tmp.a = extract(out_womble_model)
iseq.a = seq(1,by=1,length.out=81)
abline(h=seq(1,by=1,length.out=81), col='lightgrey')

counter <- 1 #indexing counter
for (i in c(1:81)) #all relevant hex areas
{
  #Plot bar in area i
  post.bar(tmp.a[,i], i=iseq.a[counter], h=0.5, a= sim_a[[i]])
  counter <- counter + 1
}

axis(2, at=iseq.a, labels = paste0(c(1:81)), las=2, cex.axis=1.7)
axis(1, at = BCADtoBP(c(-1300, -1100, -900, -700, -500, -300, -100, 100, 300, 500, 700, 900, 1100, 1300)), labels=c('1300BC','1100BC','900BC','700BC','500BC','300BC', '100BC', '100AD', '300AD', '500AD', '700AD', '900AD', '1100AD','1300AD'), tck=-0.01, cex.axis=1.7)
axis(3, at = seq(3200, 600, -200), labels=paste0(seq(3200, 600, -200),'BP'), tck=-0.01, cex.axis=1.7)
axis(1, at = BCADtoBP(c(-1200,-1000,-800, -600, -400, -200, 1, 200, 400, 600, 800, 1000, 1200)), labels=NA, tck=-0.01) #Minor tick marks
axis(3, at = seq(3200, 600, -50), labels=NA, tck=-0.01) #Minor tick marks
box()
#Legend
post.bar(c(3000,2900,2700,2600,2500,2300,2200), i=78.5, h=0.9, a=2350)
text(x=2850, y=79.5, "50% HPDI", cex=2)
text(x=2500, y=79.5,"95% HPDI", cex=2)
text(x=2600, y=77.7, "Median Posterior", cex=2)
text(x=2400, y=77, "Simulated value", cex=2)
rect(xleft=3050, xright=2150, ybottom=76, ytop=80.5, border="darkgrey", col=NA, lwd=2.5)
theme(legend.position = "none")
dev.off()


#===============================================================================
#####PAPER SECTION: The effect of the 800-400 BC calibration plateau
#===============================================================================
#SIMULATION 5: IN PLATEAU ICAR SIMULATION (with calibrated radiocarbon dates)

#Load data
load(here('data','tactical_sim_icar_withinplat.RData')) #simulated data
load(here('output', 'Womblemodel_tactical_withplat_errors.RData')) #inferred data

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
#Calculate accuracy and precision
ci_95 = credible_interval(out_womble_model, 0.95)
sim_a <- constants$true_a
sim5_model_accuracy = accuracy(sim_a, ci_95)
sim5_model_precision = precision(sim_a, ci_95)

#-------------------------------------------------------------------------------
##SUPPLEMENTARY Diagnostics: Traceplots and metric table 

#Create list to store recorded plots
plots <- vector("list", 81)

#Output
pdf(file = here('output','supplementary_figures',"traceplots_inplateau_cal.pdf"), width = 10, height = 15, onefile=TRUE)
par(mfrow = c(11,8), mar = c(2,0,2,0), oma = c(0, 0, 0, 0), mgp = c(1, 0.2, 0))
for (i in 1:81) {
  param_name <- paste0("a[", i, "]")
  plot.new()
  traceplot(out_womble_model[, param_name], main = TeX(paste0("$a[", i, "]$")), smooth = TRUE)
  plots[[i]] <- recordPlot()
}
dev.off()

#Diagnostics table
out.comb.tac_icar.model  <- do.call(rbind, out_womble_model)
post.a.model.i  <- out.comb.tac_icar.model[,paste0('a[',1:81,']')]  %>% round() 
med.model.i  <- apply(post.a.model.i, 2, median)

diagnostic_df <- data.frame(median_posterior = paste(med.model.i, "BP"),
                            HPDI95_low = paste(round(ci_95[1,]), "BP"),
                            HPDI95_high = paste(round(ci_95[2,]), "BP"),
                            rhat = round(rhat_womble_model$psrf[1:81,1],2),
                            ESS = round(ess_womble_model[1:81]))

write.csv(diagnostic_df,file=here('output','tables','diagnostics_inplateau_cal.csv'), row.names = TRUE)

#-------------------------------------------------------------------------------
##FIGURE 6C -- Map of sites, simulated arrival times and sampling window
#Plot
site_map <- ggplot(data = hex_area_win_proj) +
  geom_sf(data = sampling_win_proj, color = "grey50") +  # sampling window border
  geom_sf(aes(fill = constants$true_a)) +
  scale_fill_viridis_c(name = "Arrival time (in BP)    ", option="F", direction=-1,  limits = c(1750, 3400)) +
  guides(fill = guide_colorbar(direction = "horizontal", barwidth = 10)) + #horizontal legend
  geom_sf(data = sites_sf, size = 3, alpha = 0.5) +  # sites
  theme(
    panel.background = element_rect(fill = "lightblue", colour = "lightblue"),
    plot.margin = margin(t = 0.5 * 10, unit = "mm"), #extra space above panel for legend
    legend.position = c(0.999, 1.05),
    legend.justification = c(1, 1), 
    legend.title.position = "left",
    legend.text = element_text(size=14),
    legend.title = element_text(size=16, face="bold"),
    legend.margin = margin(t = 5, r = 30, b = 5, l = 20), 
    legend.background = element_rect(colour="grey70", size=0.5, linetype="solid"),
    axis.title = element_blank(),
    axis.text = element_text(size=15))

#Output
pdf(file=here('output','figures','sim5_sites2.pdf'), width=10, height=8.5)
grid.arrange(site_map, ncol=1, padding=0)
dev.off()

#------------
##FIGURE 6D -- Map of Inferred arrival times
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
  scale_fill_viridis_c(name = "Arrival time (in BP)    ", option="F", direction=-1,  limits = c(1750, 3400)) +
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
pdf(file=here('output','figures','sim5_arrivaltime.pdf'), width=10, height=8.5)
grid.arrange(modi, ncol=1, padding=0)
dev.off()

#-------------------------------------------------------------------------------
##Comparing Out Plateau simulation (sim 3) with In Plateau simulation (sim 5)

#Load Out Plateau data
load(here('data','tactical_sim_icar_withoutplat.RData')) #simulated data
load(here('output', 'Womblemodel_tactical_withoutplat_errors.RData')) #inferred data

#Extract values
sites_sf <- sim_dataset$sites_sf
siteInfo <- sim_dataset$siteInfo
sim_df <- sim_dataset$sim_df
constants <- sim_dataset$constants
sampling_win_proj <- sim_dataset$sampling_win_proj
hex_area_win_proj <-sim_dataset$hex_area_win_proj
tmp.a.sim3 = extract(out_womble_model)
out_womble_model3 <- out_womble_model
sim3_a <- constants$true_a
sites_sim3 <- sim_dataset$sites

#Load In Plateau data
load(here('data','tactical_sim_icar_withinplat.RData')) #simulated data
load(here('output', 'Womblemodel_tactical_withplat_errors.RData')) #inferred data

#Extract values
constants <- sim_dataset$constants
tmp.a.sim5 = extract(out_womble_model)
out_womble_model5 <- out_womble_model
sim5_a <- constants$true_a
sites_sim5 <- sim_dataset$sites

#-------------------------------------------------------------------------------
##FIGURE 7B -- Sample size per hexagon vs. precision for both Out and In Plateau simulations (sim 3 vs sim 5)

#Number of sites in area
sites5_in_areas_summarise <- sites_sim5 %>%
  group_by(area_id) %>%
  summarize(n_sites = n_distinct(site_id), .groups="drop") %>%
  complete(area_id = full_seq(min(area_id):max(area_id), 1), fill = list(n_sites = 0))

#Precision in each area
ci_95_sim3 = credible_interval(out_womble_model3, 0.95)
precision3_in_hex <- precision_in_each_area(out_womble_model3, ci_95_sim3)

ci_95_sim5= credible_interval(out_womble_model5, 0.95)
precision5_in_hex <- precision_in_each_area(out_womble_model5, ci_95_sim5)

df_precision <- data.frame(n_sites = sites5_in_areas_summarise$n_sites,
                           precision3 = precision3_in_hex,
                           precision5 = precision5_in_hex)
df_long <- df_precision %>% 
  rename("Out Plateau" = precision3, "In Plateau" = precision5) %>% 
  pivot_longer(cols = c("Out Plateau", "In Plateau"),
               names_to = "simulation",
               values_to = "precision") %>%
  mutate(simulation = factor(simulation))

#Plot
p0 <- ggplot(df_long, aes(x = n_sites, y = precision, color = simulation, shape = simulation)) +
  geom_segment(
    data = df_precision,
    aes(x = n_sites, xend = n_sites,
        y = precision5, yend = precision3),
    inherit.aes = FALSE,
    color = "grey50",
    lwd=1) +
  geom_point(size = 3) +
  labs(x = "Number of Sites",
    y = "Precision (in years)",
    color = "Model",
    shape = "Model") +
  scale_shape_manual(values = c("Out Plateau" = 16, "In Plateau" = 17)) +
  scale_color_manual(values = c("Out Plateau" = "darkblue", "In Plateau" = "darkorange")) +
  scale_x_continuous(breaks = round(seq(0, 30, by = 2),1)) +
  scale_y_continuous(breaks = round(seq(0, 1600, by = 200),1), limits =c(0,1600)) +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=2),
        axis.text = element_text(size=14),
        text = element_text(size=18),
        legend.position = c(0.96, 0.96),
        legend.justification = c(1, 1), 
        legend.title.position = "top",
        legend.text = element_text(size=14),
        legend.title = element_text(size=16, face="bold"),
        legend.margin = margin(t = 5, r = 30, b = 5, l = 20), 
        legend.background = element_rect(colour="grey70", size=0.5, linetype="solid"))

#Output
pdf(file=here('output','figures','precision_vs_nsites_sim35.pdf'), width=10, height=8)
grid.arrange(p0, ncol=1, padding=0)
dev.off()

#---------------
##SUPPLEMENTARY TABLE -- Average precision per hexagon for Out and In Plateau simulations

df_precision_diff <- df_precision %>% 
  mutate(precision_diff = precision5-precision3) %>% 
  group_by(n_sites) %>% 
  summarise(
    n_areas = n(),
    avg_precision_diff = paste(round(mean(precision_diff, na.rm = TRUE)), "years"),
    sd_precision_diff = paste(round(sd(precision_diff, na.rm = TRUE)),"years"))

write.csv(df_precision_diff,file=here('output','tables','in_out_plateau_hex_precis.csv'), row.names = FALSE)

#===============================================================================
#SIMULATION 6: IN PLATEAU ICAR SIMULATION (with uncalibrated/calendar dates)

#Load data
load(here('data','tactical_sim_icar_withinplat_noerror.RData')) #simulated data
load(here('output', 'Womblemodel_tactical_withplat_noerrors.RData')) #inferred data

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
#Calculate accuracy and precision
ci_95 = credible_interval(out_womble_model, 0.95)
sim_a <- constants$true_a
sim6_model_accuracy = accuracy(sim_a, ci_95)
sim6_model_precision = precision(sim_a, ci_95)

#-------------------------------------------------------------------------------
##SUPPLEMENTARY Diagnostics: Traceplots and metric table 

#Create list to store recorded plots
plots <- vector("list", 81)

#Output
pdf(file = here("output", "supplementary_figures", "traceplots_inplateau_noerror.pdf"), width = 10, height = 15, onefile=TRUE)
par(mfrow = c(11,8), mar = c(2,0,2,0), oma = c(0, 0, 0, 0), mgp = c(1, 0.2, 0))
for (i in 1:81) {
  param_name <- paste0("a[", i, "]")
  plot.new()
  traceplot(out_womble_model[, param_name], main = TeX(paste0("$a[", i, "]$")), smooth = TRUE)
  plots[[i]] <- recordPlot()
}
dev.off()

#Diagnostics table
out.comb.tac_icar.model  <- do.call(rbind, out_womble_model)
post.a.model.i  <- out.comb.tac_icar.model[,paste0('a[',1:81,']')]  %>% round() 
med.model.i  <- apply(post.a.model.i, 2, median)

diagnostic_df <- data.frame(median_posterior = paste(med.model.i, "BP"),
                            HPDI95_low = paste(round(ci_95[1,]), "BP"),
                            HPDI95_high = paste(round(ci_95[2,]), "BP"),
                            rhat = round(rhat_womble_model$psrf[1:81,1],2),
                            ESS = round(ess_womble_model[1:81]))

write.csv(diagnostic_df,file=here('output','tables','diagnostics_inplateau_noerror.csv'), row.names = TRUE)

#-------------------------------------------------------------------------------
##Comparing In Plateau simulation with radiocarbon dates (sim 5) with In Plateau simulation with calendar dates (sim 6)

#Load In Plateau with radiocabon dates
load(here('data','tactical_sim_icar_withinplat.RData')) #simulated data
load(here('output', 'Womblemodel_tactical_withplat_errors.RData')) #inferred data

#Extract values
sites_sf <- sim_dataset$sites_sf
siteInfo <- sim_dataset$siteInfo
sim_df <- sim_dataset$sim_df
constants <- sim_dataset$constants
sampling_win_proj <- sim_dataset$sampling_win_proj
hex_area_win_proj <-sim_dataset$hex_area_win_proj
tmp.a.sim5 = extract(out_womble_model)
out_womble_model5 <- out_womble_model
sim5_a <- constants$true_a
sites_sim5 <- sim_dataset$sites

#Load In Plateau with calendar dates
load(here('data','tactical_sim_icar_withinplat_noerror.RData')) #simulated data
load(here('output', 'Womblemodel_tactical_withplat_noerrors.RData')) #inferred data

#Extract values
constants <- sim_dataset$constants
tmp.a.sim6 = extract(out_womble_model)
out_womble_model6 <- out_womble_model
sim6_a <- constants$true_a
sites_sim6 <- sim_dataset$sites

#-------------------------------------------------------------------------------
##FIGURE 8 -- Plot and compare posteriors of arrival times for both models

pdf(file=here('output','figures','sim5_vs_sim6_posteriors_even.pdf'), width=10, height=15, pointsize=4)
par(mar = c(5, 5, 4, 2))   #pad space around plot
plot(NULL, xlim=c(3800, 1600), ylim=c(1.5, 40.5), xlab=paste('Arrival time (BP),', TeX('$a_k$')), ylab=paste('Area,', TeX('$k$')), cex.lab = 2, axes=F)
# Add Hallstatt Plateau smear
usr <- par("usr") #apply current plot limits
rect(
  xleft  = 2750, #~800BC
  xright = 2350, #~400BC
  ybottom = usr[3], ytop = usr[4], #lower and upper y limits
  col = rgb(1, 0, 0, 0.1),
  border = NA
)

iseq.a = seq(1,by=1,length.out=41)
abline(h=seq(1,by=1,length.out=41), col='lightgrey')

counter <- 1 #indexing counter
for (i in seq(1,81,2)) #all odd hex areas #for even hex areas: seq(2,80,2)
{
  #Plot bars from sim 5 and sim 6 in area i
  post.bar(tmp.a.sim5[,i], i=iseq.a[counter], h=0.5, a= sim5_a[[i]], barcolours=barcolours1)
  post.bar(tmp.a.sim6[,i], i=iseq.a[counter]+0.3, h=0.5, a= sim6_a[[i]], barcolours=barcolours2)
  counter <- counter + 1
}

axis(2, at=iseq.a, labels = paste0(seq(1,81,2)), las=2, cex.axis=1.7)
axis(1, at = BCADtoBP(c(-1900, -1700, -1500, -1300, -1100, -900, -700, -500, -300, -100, 100, 300)), labels=c('1900BC','1700BC','1500BC','1300BC', '1100BC', '900BC', '700BC', '500BC', '300BC', '100BC', '100AD', '300AD'), tck=-0.01, cex.axis=1.7)
axis(3, at = seq(3800, 1600, -200), labels=paste0(seq(3800, 1600, -200),'BP'), tck=-0.01, cex.axis=1.7)
axis(1, at = BCADtoBP(c(-1800, -1600, -1400, -1200, -1000, -800, -600, -400, -200, 1, 200)), labels=NA, tck=-0.01) #Minor tick marks
axis(3, at = seq(3800, 1600, -50), labels=NA, tck=-0.01) #Minor tick marks
box()
#Legend
post.bar(c(3900,3800,3600,3500,3400,3200,3100), i=34, h=0.9, a=3750)
text(x=3550, y=35, "50% HPDI", cex=2)
text(x=3300, y=35,"95% HPDI", cex=2)
text(x=3350, y=33, "Median Posterior", cex=2)
text(x=3700, y=33, "Simulated value", cex=2)
rect(xleft=3850, xright=3100, ybottom=31, ytop=35.5, border="darkgrey", col=NA, lwd=2.5)
text(x=3580, y=32, "Calendar dates", col="black", cex=2)
segments(x0=3800, x1=3750, y0=32, col="dodgerblue", lwd=4)
text(x=3450, y=31.5, "Calibrated radiocarbon dates", col="black", cex=2)
segments(x0=3800, x1=3750, y0=31.5, col="orchid", lwd=4)
theme(legend.position = "none")
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
##FIGURE 9A -- Accuracy vs.number of sites
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
##FIGURE 9B -- Precision vs. number of sites
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
##FIGURE 10 -- Examine clustering of sites
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
##FIGURE 11A -- Accuracy vs. sampling intensity
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
##FIGURE 11B -- Precision vs. sampling intensity
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

#===============================================================================
##SUPPLEMENTARY FIGURE: Duration parameter prior predictive check ---

nsim  <- 5000
gamma1  <- runif(nsim,1,20)
gamma2  <- rtruncnorm(nsim, mean=200, sd=100, 1, 500)
delta.mat = matrix(NA, ncol=1000, nrow=nsim)

for (i in 1:nsim) {
  delta.mat[i,] = dgamma(1:1000, gamma1[i], (gamma1[i]-1)/gamma2[i])
}

pdf(file=here('output','supplementary_figures','duration_prior_check.pdf'), height=6, width=6)
plot(NULL,xlab=TeX('$\\delta$'),ylab='Probability Density',xlim=c(1,1000),ylim=c(0,0.02))
polygon(x=c(1:1000, 1000:1), y=c(apply(delta.mat,2,quantile,prob=0.025), rev(apply(delta.mat,2,quantile,prob=0.975))), border=NA, col=rgb(0.67,0.84,0.9,0.5))
polygon(x=c(1:1000, 1000:1), y=c(apply(delta.mat,2,quantile,prob=0.25), rev(apply(delta.mat,2,quantile,prob=0.75))), border=NA, col=rgb(0.25,0.41,0.88,0.5))
legend('topright', legend=c('50% percentile range', '95% percentile range'), fill=c(rgb(0.67,0.84,0.9,0.5), rgb(0.25,0.41,0.88,0.5)))
dev.off()

# #Gamma Distribution of Simulated Duration Range
# x <- seq(0, 600, length.out = 1000)  # Range for duration
# y <- dgamma(x, shape = 5, rate = 0.04)
# plot(x, y, type = "l", lwd = 2, col = "blue", main = "Gamma Distribution of Simulated Duration", xlab = "Duration", ylab = "Density")


##SUPPLEMENTARY FIGURE: Tau prior predictive check ---
pdf(file=here('output','supplementary_figures','tau_prior_check.pdf'), height=6, width=6)
x <- seq(0, 100, length.out = 1000)
y <- dgamma(x, shape = 0.8, rate = 0.1)
plot(x, y, type = "l", lwd = 2, col = "blue", xlab = "Tau", ylab = "Probability Density")
dev.off()

