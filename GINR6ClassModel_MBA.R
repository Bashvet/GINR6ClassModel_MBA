library(R6)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(plotly)
library(kableExtra)

#PARAMETERS CLASS
Parameters <- R6Class(
  "Parameters",
  public = list(
    # Parasite biology parameters
    ma = 0.0213,      # Mortality rate of adult worms
    ml3 = 0.099,      # Mortality rate of L3 larvae
    ml4 = 0.11,       # Mortality rate of L4 larvae
    jl3 = 2.5,        # Period from ingestion to L4 (days)
    jl4 = 12,         # Period from L4 to adult worm (days)
    el = 5,           # Period from egg to L3 larvae (days)
    S = 10000,        # Ewe egg deposition rate
    
    # Immune response parameters
    tau1 = 1.2,       # Decay rate of IgA
    z = 7,            # Immune response delay (days)
    rhoA = 0.971,     # IgA response factor
    lambda1 = 3.98,   # Coefficient for IgAm regulation
    lambda2 = 1.02,   # Coefficient for worm biomass regulation
    tau2 = 0.44,      # Decay rate of ECF
    zE = 5.5,         # ECF response delay (days)
    rhoE = 0.83,      # Establishment response factor
    Eearly = 0.4,     # Establishment rate for early larvae
    Elate = 0.6,      # Establishment rate for late larvae
    
    # Worm fecundity parameters
    alpha = 1.071,    # Parameter affecting worm fecundity
    beta = 0.65,      # Effect of IgA on fecundity
    gamma = 0.00052,  # Density-dependent effect on worm length
    epsilon = 1.7,    # Larval intake rate
    omega = 10.18,    # Exponent for fecundity formula
    
    # Host growth parameters
    theta = 0.000036, # Weight gain factor
    mu = 0.614,       # Growth rate
    kappa = 0.0471,   # Saturation factor for smoother growth
    phi = 10.18,      # Weight offset
    
    # Environmental parameters
    v = 0.109,        # Herbage intake proportionality constant
    D = 50,           # Stocking density (sheep/ha)
    H = 1000,         # Herbage mass (kg DM/ha)
    
    # Treatment parameters
    treatment_efficacy = 0.57,  # 95% reduction in worm burden
    treatment_coverage = 0.8,   # 80% of flock treated
    immune_suppression = 0.5,   # 50% reduction in immune markers post-treatment
    
    
    min_establishment = 0.05,    # Minimum establishment rate
    min_fecundity = 10,          # Minimum eggs per worm
    min_worm_length = 0.1,       # Minimum worm length
    env_contamination = 100,     # Baseline environmental contamination
    background_L3 = 20,          # Background L3 contamination
    lamb_introduction = 100       # Number of lambs to introduce annually
  )
)

#SHEEP POPULATION CLASS
SheepPopulation <- R6Class(
  "SheepPopulation",
  public = list(
    parameters = NULL,
    n_sheep = NULL,
    n_timesteps = NULL,
    environmental_reservoir = NULL,
    
    # State matrices(individual-based)
    age = NULL,
    sex = NULL,
    l = NULL,
    n = NULL,
    M = NULL,
    I = NULL,
    IgAp = NULL,
    IgAm = NULL,
    L = NULL,
    L4 = NULL,
    ECF = NULL,
    E = NULL,
    W = NULL,
    Q = NULL,
    c = NULL,
    alive = NULL,
    last_treatment = NULL,
    
    initialize = function(nSheep, nTimesteps) {
      self$parameters <- Parameters$new()
      self$n_sheep <- nSheep
      self$n_timesteps <- nTimesteps
      self$environmental_reservoir <- rep(1000, nTimesteps)  # Environmental reservoir over time
      
      # Initialize state matrices with biologically plausible values
      self$age <- matrix(sample(0:24, nSheep * nTimesteps, replace = TRUE), 
                         nrow = nTimesteps, ncol = nSheep)
      self$sex <- matrix(sample(c(0, 1), nSheep * nTimesteps, replace = TRUE), 
                         nrow = nTimesteps, ncol = nSheep)
      self$l <- matrix(pmax(0.5, rnorm(nSheep * nTimesteps, mean = 0.5, sd = 0.1)), 
                       nrow = nTimesteps, ncol = nSheep)
      self$n <- matrix(pmax(20, rnorm(nSheep * nTimesteps, mean = 20, sd = 5)), 
                       nrow = nTimesteps, ncol = nSheep)
      self$M <- matrix(pmax(0, rnorm(nSheep * nTimesteps, mean = 50, sd = 10)), 
                       nrow = nTimesteps, ncol = nSheep)
      self$I <- matrix(0, nrow = nTimesteps, ncol = nSheep)
      self$IgAp <- matrix(pmax(0.05, rnorm(nSheep * nTimesteps, mean = 0.05, sd = 0.01)), 
                          nrow = nTimesteps, ncol = nSheep)
      self$IgAm <- matrix(pmax(0.6, rnorm(nSheep * nTimesteps, mean = 0.6, sd = 0.1)), 
                          nrow = nTimesteps, ncol = nSheep)
      self$L <- matrix(pmax(100, rnorm(nSheep * nTimesteps, mean = 100, sd = 20)), 
                       nrow = nTimesteps, ncol = nSheep)
      self$L4 <- matrix(pmax(100, rnorm(nSheep * nTimesteps, mean = 100, sd = 20)), 
                        nrow = nTimesteps, ncol = nSheep)
      self$ECF <- matrix(pmax(0.05, rnorm(nSheep * nTimesteps, mean = 0.05, sd = 0.01)), 
                         nrow = nTimesteps, ncol = nSheep)
      self$E <- matrix(pmax(0.1, rnorm(nSheep * nTimesteps, mean = 1, sd = 0.2)), 
                       nrow = nTimesteps, ncol = nSheep)
      self$W <- matrix(pmax(20, rnorm(nSheep * nTimesteps, mean = 30, sd = 5)), 
                       nrow = nTimesteps, ncol = nSheep)
      self$Q <- matrix(pmax(100, rnorm(nSheep * nTimesteps, mean = 500, sd = 50)), 
                       nrow = nTimesteps, ncol = nSheep)
      self$c <- matrix(0, nrow = nTimesteps, ncol = nSheep)
      self$alive <- matrix(TRUE, nrow = nTimesteps, ncol = nSheep)
      self$last_treatment <- matrix(-365, nrow = nTimesteps, ncol = nSheep) # Days since last treatment
    },
    
    apply_treatment = function(t, sheep_indices) {
      # Only treat a proportion of the flock
      n_to_treat <- round(length(sheep_indices) * self$parameters$treatment_coverage)
      treated <- sample(sheep_indices, n_to_treat)
      
      # Apply treatment effects
      self$M[t, treated] <- pmax(0, self$M[t, treated] * (1 - self$parameters$treatment_efficacy))
      self$IgAm[t, treated] <- pmax(0.01, self$IgAm[t, treated] * self$parameters$immune_suppression)
      self$IgAp[t, treated] <- pmax(0.01, self$IgAp[t, treated] * self$parameters$immune_suppression)
      self$last_treatment[t, treated] <- 0  # Just treated
      
      message("Treatment applied to ", n_to_treat, " sheep at day ", t)
    },
    
    recruit_lambs = function(t, n_lambs) {
      if (t %% 365 == 0 && t < self$n_timesteps) {  # Annual lamb introduction
        new_indices <- sample(1:self$n_sheep, n_lambs)
        self$age[t, new_indices] <- 0
        self$M[t, new_indices] <- 0  # Naive to infection
        self$IgAm[t, new_indices] <- pmax(0.1, rnorm(n_lambs, mean = 0.1, sd = 0.02))
        self$IgAp[t, new_indices] <- pmax(0.1, rnorm(n_lambs, mean = 0.1, sd = 0.02))
        self$ECF[t, new_indices] <- pmax(0.01, rnorm(n_lambs, mean = 0.01, sd = 0.005))
        self$W[t, new_indices] <- pmax(15, rnorm(n_lambs, mean = 20, sd = 3))
        message(n_lambs, " new lambs added at day ", t)
      }
    },
    
    model_worm_burden = function(t) {
      if (t > 1) {
        alive_idx <- self$alive[t, ]
        self$M[t, alive_idx] <- pmax(0, self$M[t-1, alive_idx] * (1 - self$parameters$ma))
        
        # Update days since last treatment
        self$last_treatment[t, alive_idx] <- self$last_treatment[t-1, alive_idx] + 1
        
        delay_idx <- max(1, t - self$parameters$jl4)
        l4_to_adults <- self$L4[delay_idx, alive_idx] * (1 - self$parameters$ml4)^self$parameters$jl4
        self$M[t, alive_idx] <- self$M[t, alive_idx] + pmax(0, l4_to_adults)
      }
    },
    
    model_L = function(t) {
      if (t > 1) {
        alive_idx <- self$alive[t, ]
        
        # Enhanced seasonal effect
        season_factor <- 0.15 * (1 + sin(2 * pi * (t - 30) / 365))
        
        #Host immunity affects larval survival
        immunity_factor <- 1 / (1 + pmax(0, self$IgAm[t, alive_idx]) / 10)
        
        self$L[t, alive_idx] <- pmax(0, (self$L[t-1, alive_idx] - pmax(0, self$I[t-1, alive_idx])) * 
                                       (1 - self$parameters$ml3 * immunity_factor))
        
        #background L3 contamination
        self$L[t, alive_idx] <- self$L[t, alive_idx] + self$parameters$background_L3
        
        #Add transfer from environmental reservoir
        transfer_rate <- 0.01
        self$L[t, alive_idx] <- self$L[t, alive_idx] + pmax(0, self$environmental_reservoir[t-1] * transfer_rate)
        
        delay_idx <- max(1, t - self$parameters$el)
        if (t > delay_idx) {
          # Include environmental contamination from non-sheep sources
          new_larvae <- (self$parameters$S + pmax(0, self$M[delay_idx, alive_idx]) * 
                           pmax(0, self$n[delay_idx, alive_idx])) * 
            season_factor + self$parameters$env_contamination
          
          self$L[t, alive_idx] <- self$L[t, alive_idx] + pmax(0, new_larvae)
          
          #Update environmental reservoir
          self$environmental_reservoir[t] <- pmax(0, self$environmental_reservoir[t-1] * (1 - transfer_rate) + 
                                                    pmax(0, new_larvae) * 0.2)  # 20% of new larvae go to reservoir
        }
        
        #Ensure minimum pasture contamination
        self$L[t, alive_idx] <- pmax(10, self$L[t, alive_idx])
      }
    },
    
    model_L4 = function(t) {
      if (t > 1) {
        alive_idx <- self$alive[t, ]
        self$L4[t, alive_idx] <- pmax(0, self$L4[t-1, alive_idx] * (1 - self$parameters$ml4))
        delay_idx <- max(1, t - self$parameters$jl3)
        self$L4[t, alive_idx] <- self$L4[t, alive_idx] + 
          pmax(0, self$I[delay_idx, alive_idx] * pmax(0, self$E[delay_idx, alive_idx]))
      }
    },
    
    model_immune_responses = function(t) {
      if (t > 1) {
        alive_idx <- self$alive[t, ]
        
        #Slower IgA decay when antigen is present
        decay_rate <- ifelse(pmax(0, self$M[t-1, alive_idx]) > 0, 
                             0.5^(1/(self$parameters$tau1*2)), 
                             0.5^(1/self$parameters$tau1))
        
        #Recovery of immune response after treatment
        treatment_recovery <- ifelse(self$last_treatment[t, alive_idx] < 30, 
                                     0.5,  # Slower recovery immediately post-treatment
                                     1)    # Normal rate otherwise
        
        self$IgAm[t, alive_idx] <- pmax(0, decay_rate * self$IgAm[t-1, alive_idx] * treatment_recovery)
        
        delay_idx <- max(1, t - self$parameters$z)
        self$IgAm[t, alive_idx] <- self$IgAm[t, alive_idx] + 
          self$parameters$rhoA * pmax(0, self$L4[delay_idx, alive_idx])
        
        WM <- pmax(0, self$M[t, alive_idx]) * pmax(0, self$l[t, alive_idx])
        self$IgAp[t, alive_idx] <- pmax(0, self$parameters$lambda1 * self$IgAm[t, alive_idx] + 
                                          self$parameters$lambda2 * self$IgAm[t, alive_idx] * log10(WM + 1))
        
        #ECF response with minimum baseline
        self$ECF[t, alive_idx] <- pmax(0.1, 0.5^(1/self$parameters$tau2) * self$ECF[t-1, alive_idx])
        delay_idx <- max(1, t - self$parameters$zE)
        self$ECF[t, alive_idx] <- self$ECF[t, alive_idx] + 
          self$parameters$rhoE * pmax(0, self$I[delay_idx, alive_idx])
        
        #Ensure establishment rate keeps running
        self$E[t, alive_idx] <- pmax(self$parameters$min_establishment, 
                                     (self$parameters$Eearly - self$parameters$Elate) * 
                                       exp(-pmax(0, self$ECF[t, alive_idx])) + self$parameters$Elate)
      }
    },
    
    model_worm_fecundity = function(t) {
      if (t > 1) {
        alive_idx <- self$alive[t, ]
        # Ensure minimum worm length and fecundity
        self$l[t, alive_idx] <- pmax(self$parameters$min_worm_length, 
                                     self$parameters$alpha - 
                                       self$parameters$beta * log10(pmax(0.01, self$IgAm[t, alive_idx])) - 
                                       self$parameters$gamma * pmax(0, self$M[t, alive_idx]))
        
        #Density-dependent fecundity with minimum output
        self$n[t, alive_idx] <- pmax(self$parameters$min_fecundity, 
                                     (self$parameters$epsilon * 
                                        pmax(0, self$l[t, alive_idx])^self$parameters$omega - 1) * 500)
      }
    },
    
    model_host_growth = function(t) {
      if (t > 1) {
        alive_idx <- self$alive[t, ]
        #Growth affected by parasite burden
        parasite_effect <- 1 - (pmax(0, self$M[t, alive_idx]) / (100 + pmax(0, self$M[t, alive_idx])))
        
        self$W[t, alive_idx] <- (self$parameters$theta * 
                                   exp(self$parameters$mu * (1 - exp(-self$parameters$kappa * t))) / 
                                   self$parameters$kappa + self$parameters$phi) * parasite_effect
        
        self$Q[t, alive_idx] <- self$parameters$v * (self$W[t, alive_idx] - self$parameters$phi)
        
        #Age the sheep
        self$age[t, alive_idx] <- self$age[t-1, alive_idx] + 1
      }
    },
    
    model_infection = function(t) {
      alive_idx <- self$alive[t, ]
      self$I[t, alive_idx] <- pmax(0, self$L[t, alive_idx] * self$Q[t, alive_idx] * 
                                     self$parameters$D / max(1, self$parameters$H))
    },
    
    model_egg_count = function(t) {
      alive_idx <- self$alive[t, ]
      self$c[t, alive_idx] <- pmax(0, self$M[t, alive_idx] * self$n[t, alive_idx] / max(1, self$W[t, alive_idx] * 20))
    },
    
    run_timestep = function(t) {
      if (t > 1) {
        self$model_host_growth(t)
        self$model_infection(t)
        self$model_immune_responses(t)
        self$model_L4(t)
        self$model_worm_burden(t)
        self$model_worm_fecundity(t)
        self$model_L(t)
        self$model_egg_count(t)
        self$recruit_lambs(t, self$parameters$lamb_introduction)
      }
    },
    
    run_simulation = function(treatment_interval = NULL) {
      for (t in 2:self$n_timesteps) {
        self$run_timestep(t)
        
        if (!is.null(treatment_interval) && t %% treatment_interval == 0) {
          self$apply_treatment(t, 1:self$n_sheep)
        }
      }
    },
    
    #Helper function to get summary results
    get_results = function() {
      #Calculate means only for alive sheep
      list(
        WormBurden = sapply(1:self$n_timesteps, function(t) {
          if (any(self$alive[t, ])) {
            mean(self$M[t, self$alive[t, ]])
          } else {
            0
          }
        }),
        FEC = sapply(1:self$n_timesteps, function(t) {
          if (any(self$alive[t, ])) {
            mean(self$c[t, self$alive[t, ]])
          } else {
            0
          }
        }),
        PastureL3 = sapply(1:self$n_timesteps, function(t) {
          if (any(self$alive[t, ])) {
            mean(self$L[t, self$alive[t, ]])
          } else {
            0
          }
        }),
        Weight = sapply(1:self$n_timesteps, function(t) {
          if (any(self$alive[t, ])) {
            mean(self$W[t, self$alive[t, ]])
          } else {
            0
          }
        }),
        Immunity = sapply(1:self$n_timesteps, function(t) {
          if (any(self$alive[t, ])) {
            mean(self$IgAm[t, self$alive[t, ]])
          } else {
            0
          }
        }),
        Establishment = sapply(1:self$n_timesteps, function(t) {
          if (any(self$alive[t, ])) {
            mean(self$E[t, self$alive[t, ]])
          } else {
            0
          }
        }),
        Alive = sapply(1:self$n_timesteps, function(t) sum(self$alive[t, ]))
      )
    },
    
    #Individual sheep data for detailed analysis
    get_individual_data = function() {
      data.frame(
        SheepID = rep(1:self$n_sheep, each = self$n_timesteps),
        Day = rep(1:self$n_timesteps, self$n_sheep),
        Age = c(t(self$age)),
        WormBurden = c(t(self$M)),
        FEC = c(t(self$c)),
        Weight = c(t(self$W)),
        Immunity = c(t(self$IgAm)),
        LastTreatment = c(t(self$last_treatment)),
        Alive = c(t(self$alive))
      )
    }
  )
)

#SIMULATION SETUP
n_sheep <- 100
n_timesteps <- 730  #2 years of daily timesteps
treatment_intervals <- c(30, 60, 90)  #Treatment frequencies (30d, 60d, 90d)

#RUN SIMULATIONS
run_scenario <- function(interval) {
  pop <- SheepPopulation$new(nSheep = n_sheep, nTimesteps = n_timesteps)
  pop$run_simulation(treatment_interval = interval)
  return(list(
    summary = pop$get_results(),
    individual = pop$get_individual_data() %>% mutate(Group = paste0(interval, "d"))
  ))
}  

#Run control scenario (no treatment)
control_pop <- SheepPopulation$new(n_sheep, n_timesteps)
control_pop$run_simulation()
control_results <- list(
  summary = control_pop$get_results(),
  individual = control_pop$get_individual_data() %>% mutate(Group = "Control")
)

#Run treatment scenarios
treatment_results <- lapply(treatment_intervals, function(int) {
  run_scenario(int)
})
names(treatment_results) <- paste0(treatment_intervals, "d")

#Combine all results
all_summary_results <- c(list(Control = control_results$summary), 
                         lapply(treatment_results, function(x) x$summary))

all_individual_results <- bind_rows(
  control_results$individual,
  bind_rows(lapply(treatment_results, function(x) x$individual)))
