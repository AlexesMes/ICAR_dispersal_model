#Load Libraries ----
library(here)
library(coda)
library(nimbleCarbon)
library(rcarbon)
library(dplyr)
library(parallel)
library(spdep)
library(rlist)

rm(list = ls())
`%!in%` <- Negate(`%in%`)

set.seed(123)

##ICAR model tactical simulation with uncalibrated dates -- varying number of sites

#-------------------------------------------------------------------------------
## Data and Functions Setup ----
source(here('src', 'sim_data.R'))
load(here('data','trig.RData')) #nodes and edges between hex area centroids

#-------------------------------------------------------------------------------
output_varysitenumber_models <- list()

#Varying sample intensity (clustering in simulated data)
for (n in seq(200,2000,200)){
  #Generate simulated data set
  sim_dataset <- sim_data(with_calibration = FALSE, 
                          seed=123,
                          structure = 'CAR',
                          k = 0.5, 
                          n_sites = n, 
                          n_dates = 3*n,
                          x1_areas = c(7, 15, 23, 32, 39, 44, 6, 10, 18, 27, 35, 14, 31),
                          x2_areas = c(25, 39, 51, 44, 63, 57, 56, 52), 
                          beta1 = -400, 
                          beta2 = 300)
  
  sites <- sim_dataset$sites
  sites_sf <- sim_dataset$sites_sf
  siteInfo <- sim_dataset$siteInfo
  sim_df <- sim_dataset$sim_df
  constants <- sim_dataset$constants
  sampling_win_proj <- sim_dataset$sampling_win_proj
  hex_area_win_proj <-sim_dataset$hex_area_win_proj
  
  #Combine constants
  constants <- c(constants, constants_trig)
  #-------------------------------------------------------------------------------
  ##Simulated data summary
  #Number of dates in area
  dates_in_areas_summarise <- as.data.frame(table(sites$area_id))
  #Number of sites in area
  sites_in_areas_summarise <- sites %>% group_by(area_id) %>% summarize(n_sites =n_distinct(site_id))
  
  #-------------------------------------------------------------------------------
  ## Initialise Parameters ----
  
  # Initialise regional parameters ----
  buffer <- 100
  delta_init <- siteInfo$diff + buffer
  alpha_init <- siteInfo$earliest + buffer/2
  
  # Initialise hex areas which contain sites
  init_a  <- aggregate(earliest~area_id, FUN=max, data=siteInfo) #find earliest date in each region k
  init_b  <- aggregate(latest~area_id, FUN=min, data=siteInfo) #find latest date in each area k
  
  #Initialise hex areas which do not contain sites
  init_empty_area <- function(init_df) {
    for(i in 1:constants$n_area){
      
      area_ids <- init_df$area_id #List of hex areas ids with sites
      
      if (i %!in% area_ids){
        empty_hex_id <- i #Id of empty hex
        neighbour_hex <- which.min(abs(i - area_ids)) #Determine closest hex neighbor which has sites. If there are more than one neighbour hex with sites, it selects the first observation (i.e. hex with the smallest id, since ids are in ascending order). 
        neighbour_hex_id <- area_ids[neighbour_hex] #Determine area id of closest hex neighbor
        neighbour_date <- init_df[neighbour_hex , 2] #Select the date associated with the neighbour hex 
        init_df <- rbind(init_df, c(i, neighbour_date)) #Assign this date to the empty hex
      }
    }
    return(init_df)
  }
  init_a <- init_a %>% init_empty_area() %>%  arrange(area_id)
  init_b <- init_b %>% init_empty_area() %>%  arrange(area_id)
  
  ##Delta parameter initialisation
  init_nabla <- 0
  for (t in 1:constants$n_trans){
    m <- constants$edge_id1[t] #transition/edge t, select first area
    n <- constants$edge_id2[t] #transition/edge t, select second area
    
    init_nabla[t] <- (init_a$earliest[init_a$area_id==m] - init_a$earliest[init_a$area_id==n])
  }
  init_nabla_phi <- init_nabla*0.5
  
  #Add buffer
  init_a  <- init_a[ ,2] + buffer
  init_b  <- init_b[ ,2] - buffer
  
  # Initialise spatial residues
  init_phi <- init_a
  
  #-------------------------------------------------------------------------------
  #Spatial data ----
  
  library(spdep)
  nb_areas <- poly2nb(as(hex_area_win_proj, 'Spatial'), queen=FALSE, row.names = hex_area_win_proj$area_ID) #neighboring areas using sp library 
  #nb_areas <- st_intersects(hex_area_win, hex_area_win, remove_self = TRUE) #neighboring areas using sf library
  
  nbInfo <- nb2WB(nb_areas) #transform into iCAR inputs: adjacent matrix, weights, number of neighbors (for WinBUGS)
  
  #-------------------------------------------------------------------------------
  #Constants ----
  constants$adj <- nbInfo$adj
  constants$weights <- nbInfo$weights
  constants$num <- nbInfo$num
  constants$L <- length(nbInfo$adj)
  constants <- constants[names(constants) %!in% c("dist_mat", 
                                                  "dist_org", 
                                                  "center_coords",
                                                  "beta0",
                                                  "beta1",
                                                  "beta2",
                                                  "x1",
                                                  "x2",
                                                  "mu",
                                                  "mu2",
                                                  "tau.err",
                                                  "tau")] #remove constants which aren't used
  
  # #-------------------------------------------------------------------------------
  ## Model W: ICAR and wombling integrating sample interdependence, i.e. the addition of a hierarchical model ----
  
  modelW <- nimbleCode({
    #For Each Site
    for (j in 1:n_sites)
    {
      delta[j] ~ dgamma(gamma1, (gamma1-1)/gamma2)
      alpha[j] ~ dunif(max = a[id_areas[j]], min = b[id_areas[j]]);
    }
    
    #For Each Date
    for (i in 1:n_dates){
      theta[i] ~ dunif(min = (alpha[id_sites[i]] - (delta[id_sites[i]]+1)), max = alpha[id_sites[i]]);
    }
    
    #For Each Region
    for (k in 1:n_areas){
      a[k] <- phi[k];
      b[k] ~ dunif(1000, 6500);
      constraint_uniform[k] ~ dconstraint(b[k]<a[k]); #In each area, start date of occupation, a_k, must be greater than the end date of occupation, b_k (note: BP dates in the positive direction)
    }
    
    # ICAR Model prior to capture spatial random effects
    phi[1:n_areas] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n_areas], tau1, zero_mean=0)
    
    #For Each Boundary
    for (t in 1:n_trans){
      #nabla defines the difference in arrival time across a boundary
      nabla[t] <- abs(a[edge_id1[t]] - a[edge_id2[t]]) #edge t: select first area, m, and second area, n
      #nabla_phi defines the difference in spatial residues across a boundary
      nabla_phi[t] <- abs(phi[edge_id1[t]] - phi[edge_id2[t]])
    }
    
    #Priors
    tau1 ~ dgamma(0.8, 0.1);  #weak prior for ICAR model -- spatial autocorrelation precision parameter
    
    # Hyperprior for duration
    gamma1 ~ dunif(1,20); #Hyperprior for rate
    gamma2 ~ T(dnorm(mean=200, sd=100), 1, 500) #Hyperprior for mode
    
  })
  
  
  #Define initial values ----
  dW <- list(theta=sim_df$cra, 
             constraint_uniform = rep(1, constants$n_areas))
  
  
  initsW <- list(b=init_b,
                 alpha=alpha_init,
                 delta=delta_init,
                 phi=init_phi,
                 tau1=rgamma(1, shape = 0.8, rate = 0.1),
                 gamma1=10,
                 gamma2=200)
  
  
  #Run MCMC ----
  out_womble_model <- nimbleMCMC(code = modelW,
                                 constants = constants,
                                 data = dW,
                                 niter = 500000,
                                 nchains = 4,
                                 thin= 100,
                                 nburnin = 250000,
                                 monitors = c('a', 'b', 'theta', 'nabla', 'nabla_phi', 'delta', 'alpha', 'phi'), 
                                 inits = initsW,
                                 samplesAsCodaMCMC=TRUE)
  
  #Diagnostics ----
  rhat_womble_model  <- gelman.diag(out_womble_model, multivariate = FALSE)
  ess_womble_model  <- effectiveSize(out_womble_model)

  #-----------------------------------------------------------------------------
  #Output ----
  output <- list(out_womble_model,
                 rhat_womble_model,
                 ess_womble_model,
                 constants,
                 n) #number of sites 

  output_varysitenumber_models <- list.append(output_varysitenumber_models, output)
  
}

#-------------------------------------------------------------------------------
# Save output ----
save(output_varysitenumber_models, file=here('output','Womblemodel_tactical_icar_sitenumber.RData'))
