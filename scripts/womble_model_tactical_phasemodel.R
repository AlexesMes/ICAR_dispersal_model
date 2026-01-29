#Load Libraries ----
library(here)
library(coda)
library(nimbleCarbon)
library(rcarbon)
library(dplyr)
library(parallel)
library(spdep)

rm(list = ls())
`%!in%` <- Negate(`%in%`)

set.seed(123)

##Phase model tactical simulation with calibrated radiocarbon dates

#-------------------------------------------------------------------------------
## Data and Functions Setup ----
source(here('src', 'sim_data.R'))
load(here('data','trig.RData')) #nodes and edges between hex area centroids

# Generate simulated data set
sim_dataset <- sim_data(with_calibration = TRUE,
                        seed=10,
                        structure="CAR",
                        k = 0.3,
                        n_sites = 800,
                        n_dates = 3*800,
                        beta1=-500, 
                        beta2=0, 
                        x1_areas=c(1,2,6,10,11,15,19,20,24,25,29,38), 
                        x2_areas=0,
                        a_min=0,
                        a_max=3000,
                        mu1=1500)

# Save output
save(sim_dataset, file=here('data','tactical_sim_icar.RData'))

#-------------------------------------------------------------------------------
#Extract values
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

#Initialise regional parameters
buffer <- 100
delta_init <- siteInfo$diff + buffer
alpha_init <- siteInfo$earliest + buffer/2

#Initialise hex areas which contain sites
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

#Delta parameter initialisation
init_nabla <- 0
for (t in 1:constants$n_trans){
  m <- constants$edge_id1[t] #transition/edge t, select first area
  n <- constants$edge_id2[t] #transition/edge t, select second area
  
  init_nabla[t] <- (init_a$earliest[init_a$area_id==m] - init_a$earliest[init_a$area_id==n])
}

#Add buffer
init_a  <- init_a[ ,2] + buffer
init_b  <- init_b[ ,2] - buffer

# Initialise spatial residues
init_phi <- init_a

#-------------------------------------------------------------------------------
#Spatial data ----
library(spdep)
nb_areas <- poly2nb(as(hex_area_win_proj, 'Spatial'), queen=FALSE, row.names = hex_area_win_proj$area_ID) #neighboring areas using sp library 
nbInfo <- nb2WB(nb_areas) #transform into iCAR inputs: adjacent matrix, weights, number of neighbors (for WinBUGS)

#-------------------------------------------------------------------------------
##General Setup ----
#Data
dat <- list(cra=sim_df$cra,
            cra_error=sim_df$cra_error,
            constraint_uniform = rep(1, constants$n_areas),
            cra_constraint = rep(1, constants$n_dates)) # Set-up constraint for ignoring inference outside calibration range

theta_init <- medCal(calibrate(dat$cra, dat$cra_error, verbose=FALSE)) #initialize theta parameter

#Constants
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
#Calibration constants (with dummy extensions)
data("intcal20") 
constants$calBP <- c(1000000, intcal20$CalBP, -1000000)
constants$C14BP <- c(1000000, intcal20$C14Age, -1000000)
constants$C14err <- c(1000, intcal20$C14Age.sigma, 1000)

# #-------------------------------------------------------------------------------
## Model W: ICAR and wombling integrating sample interdependence, i.e. the addition of a hierarchical model ----

modelW <- function(seed, d, theta_init, alpha_init, delta_init, init_a, init_b, constants, nburnin, thin, niter)
{
  #Load Library
  library(nimbleCarbon)
  #Define Core Model
  model <- nimbleCode({
    #For Each Site
    for (j in 1:n_sites)
    {
      delta[j] ~ dgamma(gamma1, (gamma1-1)/gamma2)
      alpha[j] ~ dunif(max = a[id_areas[j]], min = b[id_areas[j]]);
    }
    
    #For Each Date
    for (i in 1:n_dates){
      theta[i] ~ dunif(min = (alpha[id_sites[i]] - (delta[id_sites[i]]+1)), max = alpha[id_sites[i]]);
      # Calibration
      mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]); #C14 age on the relevant calibration curve
      cra_constraint[i] ~ dconstraint(mu[i] < 50193 & mu[i] > 95) #C14 age must be within the calibration range
      sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]); #error on the calibration curve
      error[i] <- (cra_error[i]^2 + sigmaCurve[i]^2)^(1/2); #the samples' C14 error + error on calibration curve
      cra[i] ~ dnorm(mean=mu[i], sd=error[i]); #observed radiocarbon age of the sample
    }
    
    #For Each Region
    for (k in 1:n_areas){
      a[k] ~ dunif(50, 3000);
      b[k] ~ dunif(50, 3000);
      constraint_uniform[k] ~ dconstraint(a[k]>b[k]); #In each area, start date of occupation, a_k, must be greater than the end date of occupation, b_k (note: BP dates in the positive direction)
    }
    
    #For Each Boundary
    for (t in 1:n_trans){
      #nabla defines the difference in arrival time across a boundary
      nabla[t] <- abs(a[edge_id1[t]] - a[edge_id2[t]]) #edge t: select first area, m, and second area, n
    }
    
    # Hyperprior for duration
    gamma1 ~ dunif(1,20); #Hyperprior for rate
    gamma2 ~ T(dnorm(mean=200, sd=100), 1, 500) #Hyperprior for mode
    
  })
  
  #Define initial values ----
  inits <- list(a=init_a,
                b=init_b,
                alpha=alpha_init,
                delta=delta_init,
                theta=d$cra, #theta_init
                gamma1=10,
                gamma2=200)
  
  # Compile and Run model	----
  #Create model object
  model <- nimbleModel(model, constants=constants, data=d, inits=inits)
  cModel <- compileNimble(model)
  #Configure MCMC with conjugacy where possible
  conf <- configureMCMC(model, useConjugacy = TRUE, control = list(adaptInterval=5000, adaptFactorExponent=0.1))
  #Add monitors
  conf$addMonitors(c('a','b','theta','delta','alpha','nabla'))
  #Build, compile, and run MCMC
  MCMC <- buildMCMC(conf)
  cMCMC <- compileNimble(MCMC)
  results <- runMCMC(cMCMC, niter = niter, thin = thin, nburnin = nburnin, samplesAsCodaMCMC = T, setSeed = seed) 
}

#-------------------------------------------------------------------------------
# MCMC Setup
ncores  <-  4
cl <- makeCluster(ncores)
seeds <- c(12, 34, 56, 78)
niter  <- 2000000
nburnin  <- 1000000
thin  <-100

#Hierarchical Womble Model 
out_womble_model <-  parLapply(cl = cl, 
                               X = seeds, 
                               fun = modelW, 
                               d = dat, 
                               constants = constants, 
                               theta_init=dat$cra,
                               init_a = init_a, 
                               init_b = init_b, 
                               alpha_init = alpha_init, 
                               delta_init = delta_init,
                               niter = niter, 
                               nburnin = nburnin,
                               thin = thin)
out_womble_model <- mcmc.list(out_womble_model)

## Diagnostics ----
rhat_womble_model  <- gelman.diag(out_womble_model, multivariate = FALSE)
ess_womble_model  <- effectiveSize(out_womble_model)
agg_womble_model <- agreementIndex(dat$cra,
                                   dat$cra_error,
                                   calCurve = "intcal20",
                                   theta = out_womble_model[[1]][ , grep("theta", colnames(out_womble_model[[1]]))],
                                   verbose = F)

#-------------------------------------------------------------------------------
# Save output ----
save(out_womble_model, 
     rhat_womble_model, 
     ess_womble_model, 
     agg_womble_model, 
     file=here('output','Womblemodel_tactical_phasemodel.RData'))