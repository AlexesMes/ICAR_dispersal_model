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
rect(xleft=6850, xright=6150, ybottom=75, ytop=79, border="darkgrey", col="white", lwd=2)
post.bar(c(6900,6800,6600,6500,6400,6200,6100), i=77, h=0.9, a=6750)
text(x=6550, y=78, "50% HPDI", cex=1.5)
text(x=6300, y=78,"95% HPDI", cex=1.5)
text(x=6450, y=76, "Median Posterior", cex=1.5)
text(x=6700, y=76, "Simulated value", cex=1.5)
theme(legend.position = "none")
dev.off()

#-------------------------------------------------------------------------------
##SUPPLEMENTARY FIGURE: Accuracy and precision in each hexagonal area vs. number of sites

#Number of sites in area
sites_in_areas_summarise <- sites %>% 
  group_by(area_id) %>% 
  summarize(n_sites = n_distinct(site_id), .groups="drop") %>%
  complete(area_id = full_seq(1:max(area_id), 1), fill = list(n_sites = 0))

#Accuracy and Precision in each hex
df_woa_precis_acc <- data.frame(n_sites = sites_in_areas_summarise$n_sites,
                                precision = sim1_precision_in_hex,
                                accuracy = sim1_accuracy_in_hex)

p1 <- ggplot(df_woa_precis_acc, aes(x = n_sites, y = precision, color = accuracy)) + 
  geom_point(size = 3.5) +
  labs(x = "Number of Sites",
       y = "Precision (in years)",
       color = "Accuracy (95% HPDI)") +
  scale_shape_manual(values = 16) +
  scale_color_manual(values = c("FALSE" = "firebrick2", "TRUE"= "blue")) +
  scale_x_continuous(breaks = round(seq(0, 35, by = 5),1)) +
  scale_y_continuous(breaks = round(seq(0, 1000, by = 200),1), limits =c(0,1000)) +
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

pdf(file=here('output','supplementary_figures','precision_vs_nsites_sim1_woa_acc.pdf'), width=10, height=8)
grid.arrange(p1, ncol=1, padding=0)
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
    panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
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
##SUPPLEMENTARY FIGURE -- Posterior distributions of arrival times

pdf(file=here('output','supplementary_figures','sim2_vs_sim7_posteriors_odd.pdf'), width=10, height=15, pointsize=4)
par(mar = c(5, 5, 4, 2))   #pad space around plot
plot(NULL, xlim=c(3250, 500), ylim=c(1.5, 40.5), xlab=paste('Arrival time (BP),', TeX('$a_k$')), ylab=paste('Area,', TeX('$k$')), cex.lab = 2, axes=F)

iseq.a = seq(1,by=1,length.out=41)
abline(h=seq(1,by=1,length.out=41), col='lightgrey')

counter <- 1 #indexing counter
for (i in seq(1,81,2)) #all even hex areas #for odd hex areas use: seq(1,81,2)
{
  #Plot bars from sim 2 and sim 7 in area i
  post.bar(tmp.a.sim2[,i], i=iseq.a[counter], h=0.5, a= sim2_a[[i]], barcolours=barcolours1)
  post.bar(tmp.a.sim7[,i], i=iseq.a[counter]+0.3, h=0.5, a= sim7_a[[i]], barcolours=barcolours2)
  counter <- counter + 1
}

axis(2, at=iseq.a, labels = paste0(seq(1,81,2)), las=2, cex.axis=1.7) #for odd hex areas use: seq(1,81,2)
axis(1, at = BCADtoBP(c(-1300, -1100, -900, -700, -500, -300, -100, 100, 300, 500, 700, 900, 1100, 1300)), labels=c('1300BC','1100BC','900BC','700BC','500BC','300BC', '100BC', '100AD', '300AD', '500AD', '700AD', '900AD', '1100AD','1300AD'), tck=-0.01, cex.axis=1.7)
axis(3, at = seq(3200, 500, -200), labels=paste0(seq(3200, 500, -200),'BP'), tck=-0.01, cex.axis=1.7)
axis(1, at = BCADtoBP(c(-1200,-1000,-800, -600, -400, -200, 1, 200, 400, 600, 800, 1000, 1200, 1400)), labels=NA, tck=-0.01) #Minor tick marks
axis(3, at = seq(3200, 500, -50), labels=NA, tck=-0.01) #Minor tick marks
box()
#Legend
rect(xleft=3150, xright=2300, ybottom=8, ytop=11.5, border="darkgrey", col="white", lwd=2.5)
post.bar(c(3100,3000,2800,2700,2600,2400,2300), i=10.5, h=0.9, a=2550)
text(x=2600, y=11, "50% HPDI", cex=2)
text(x=2950, y=11,"95% HPDI", cex=2)
text(x=2700, y=10, "Median Posterior", cex=2)
text(x=2550, y=9.5, "Simulated value", cex=2)
segments(y0=10.3, y1=9.6, x0=2550, col="darkgrey", lwd=2)
text(x=2750, y=9, "ICAR Model", col="black", cex=2)
segments(x0=3050, x1=3000, y0=9, col="dodgerblue", lwd=4)
text(x=2750, y=8.5, "Phasemodel", col="black", cex=2)
segments(x0=3050, x1=3000, y0=8.5, col="orchid", lwd=4)
theme(legend.position = "none")
dev.off()

#-------------------------------------------------------------------------------
##SUPPLEMENTARY FIGURE -- (A) ICAR Precision and accuracy vs. number of sites

#Number of sites in area
sites2_in_areas_summarise <- sites_sim2 %>% 
  group_by(area_id) %>% 
  summarize(n_sites = n_distinct(site_id), .groups="drop") %>%
  complete(area_id = full_seq(1:max(area_id), 1), fill = list(n_sites = 0))

#Accuracy (hit/miss) and precision in each area
ci_95_sim2 = credible_interval(out_womble_model2, 0.95)
accuracy2_in_hex <- accuracy_in_each_area(sim2_a, ci_95_sim2)
precision2_in_hex <- precision_in_each_area(out_womble_model2, ci_95_sim2)

df_ICAR_precis_acc <- data.frame(n_sites = sites2_in_areas_summarise$n_sites,
                                 precision = precision2_in_hex,
                                 accuracy = accuracy2_in_hex)

p1 <- ggplot(df_ICAR_precis_acc, aes(x = n_sites, y = precision, color = accuracy)) + 
  geom_point(size = 3.5) +
  labs(x = "Number of Sites",
       y = "Precision (in years)",
       color = "Accuracy (95% HPDI)") +
  scale_shape_manual(values = 16) +
  scale_color_manual(values = c("FALSE" = "firebrick2", "TRUE"= "blue")) +
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
        legend.margin = margin(t = 5, r = 30, b = 5, l = 20), 
        legend.background = element_rect(colour="grey70", size=0.5, linetype="solid"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=2))

pdf(file=here('output','supplementary_figures','precision_vs_nsites_sim27_icar_acc.pdf'), width=10, height=8)
grid.arrange(p1, ncol=1, padding=0)
dev.off()

#---------------
##SUPPLEMENTARY FIGURE -- (B) Phasemodel Precision and accuracy vs. number of sites

#Accuracy (hit/miss) and precision in each area
ci_95_sim7 = credible_interval(out_womble_model7, 0.95)
accuracy7_in_hex <- accuracy_in_each_area(sim7_a, ci_95_sim7)
precision7_in_hex <- precision_in_each_area(out_womble_model7, ci_95_sim7)

df_Phase_precis_acc <- data.frame(n_sites = sites2_in_areas_summarise$n_sites,
                                  precision = precision7_in_hex,
                                  accuracy = accuracy7_in_hex)

p2 <- ggplot(df_Phase_precis_acc, aes(x = n_sites, y = precision, color = accuracy)) + 
  geom_point(size = 3.5) +
  labs(x = "Number of Sites",
       y = "Precision (in years)",
       color = "Accuracy (95% HPDI)") +
  scale_shape_manual(values = 17) +
  scale_color_manual(values = c("FALSE" = "firebrick2", "TRUE"= "blue")) +
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
        legend.margin = margin(t = 5, r = 30, b = 5, l = 20),
        legend.background = element_rect(colour="grey70", size=0.5, linetype="solid"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=2))

pdf(file=here('output','supplementary_figures','precision_vs_nsites_sim27_phase_acc.pdf'), width=10, height=8)
grid.arrange(p2, ncol=1, padding=0)
dev.off()

#---------------
##SUPPLEMENTARY TABLE -- Average precision per hexagon for icar and phasemodel

#Precision
df_precision <- data.frame(n_sites = sites2_in_areas_summarise$n_sites,
                           precision2 = precision2_in_hex,
                           precision7 = precision7_in_hex)


df_precision_diff <- df_precision %>% 
  mutate(precision_diff = precision2-precision7) %>% 
  group_by(n_sites) %>% 
  summarise(
    n_areas = n(),
    avg_precision_diff = paste(round(mean(precision_diff, na.rm = TRUE)), "years"),
    sd_precision_diff = paste(round(sd(precision_diff, na.rm = TRUE)),"years"))

write.csv(df_precision_diff,file=here('output','tables','icar_phase_hex_precis.csv'), row.names = FALSE)

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