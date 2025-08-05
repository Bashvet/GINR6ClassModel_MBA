# Model calibration
library(R6)
library(dplyr)
library(tidyr)
library(ggplot2)
library(abc)
library(sensitivity)

target_stats <- list(
  mean_FEC = 5.77,
  mean_IgA = 36.31,
  mean_M   = 115.43
)

prior_ranges <- data.frame(
  rhoA = c(0.5, 1.5),
  rhoE = c(0.5, 1.5),
  beta = c(0.4, 1.0),
  gamma = c(0.0001, 0.001)
)

sample_priors <- function(n) {
  data.frame(
    rhoA = runif(n, prior_ranges$rhoA[1], prior_ranges$rhoA[2]),
    rhoE = runif(n, prior_ranges$rhoE[1], prior_ranges$rhoE[2]),
    beta = runif(n, prior_ranges$beta[1], prior_ranges$beta[2]),
    gamma = runif(n, prior_ranges$gamma[1], prior_ranges$gamma[2])
  )
}
# Parameters
set_particle_parameters <- function(particle, parameters) {
  parameters$rhoA <- particle$rhoA
  parameters$rhoE <- particle$rhoE
  parameters$beta <- particle$beta
  parameters$gamma <- particle$gamma
  return(parameters)
}

# Simulation
run_particle_simulation <- function(particle) {
  parameters <- Parameters$new()
  parameters <- set_particle_parameters(particle, parameters)
  
  pop <- SheepPopulation$new(nSheep = 100, nTimesteps = 365 * 20)
  pop$parameters <- parameters
  pop$run_simulation()
  
  results <- pop$get_results()
  last_year <- (length(results$FEC) - 365 + 1):length(results$FEC)
  
  stats <- list(
    mean_FEC = mean(results$FEC[last_year]),
    mean_IgA = mean(results$Immunity[last_year]),
    mean_M   = mean(results$WormBurden[last_year])
  )
  
  return(stats)
}

# Distance Function
compute_distance <- function(stats, target_stats) {
  sqrt(
    (stats$mean_FEC - target_stats$mean_FEC)^2 +
      (stats$mean_IgA - target_stats$mean_IgA)^2 +
      (stats$mean_M   - target_stats$mean_M)^2
  )
}

# ABC Calibration
abc_calibration <- function(n_particles = 100, tolerance = 0.01) {
  set.seed(123)
  priors <- sample_priors(n_particles)
  distances <- numeric(n_particles)
  summary_stats <- vector("list", n_particles)
  
  for (i in 1:n_particles) {
    stats <- tryCatch(run_particle_simulation(priors[i, ]), error = function(e) NULL)
    if (!is.null(stats)) {
      summary_stats[[i]] <- stats
      distances[i] <- compute_distance(stats, target_stats)
    } else {
      distances[i] <- Inf
    }
    cat("Completed particle", i, "\n")
  }
  
  accept_cutoff <- quantile(distances, tolerance)
  accepted_indices <- which(distances <= accept_cutoff)
  accepted_params <- priors[accepted_indices, ]
  accepted_stats <- do.call(rbind, summary_stats[accepted_indices])
  
  abc_fit <- abc(
    target = unlist(target_stats),
    param = accepted_params,
    sumstat = accepted_stats,
    method = "loclinear"
  )
  
  return(abc_fit$adj.values)
}

# Run abc
posterior_samples <- abc_calibration(n_particles = 100)  # Reduced for testing

# Visualize
posterior_df <- as.data.frame(posterior_samples)
posterior_long <- pivot_longer(posterior_df, everything(), names_to = "Parameter", values_to = "Value")

ggplot(posterior_long, aes(x = Value)) +
  geom_density(fill = "skyblue") +
  facet_wrap(~Parameter, scales = "free", ncol = 2) +
  theme_minimal() +
  labs(title = "Posterior Distributions from ABC Calibration")

# Sobol Sensitivity Analysis
sobol_input <- model.matrix(~ ., data = posterior_df)

sobol_result <- sobol2007(
  model = function(X) {
    apply(X, 1, function(x) {
      particle <- as.list(setNames(x, colnames(posterior_df)))
      stats <- tryCatch(run_particle_simulation(particle), error = function(e) NULL)
      if (is.null(stats)) return(NA) else return(stats$mean_FEC)
    })
  },
  X1 = sobol_input[1:50, ],
  X2 = sobol_input[51:100, ],
  nboot = 100
)

print(sobol_result)
plot(sobol_result)
