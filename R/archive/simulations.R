###############################################################################
# Notes: These simulations uses the trend model with weights for dispersion
# simulations. This will produce more variable simulated dispersion values.
#
# 01. Source data and functions
# 02. Simulations and models
#
#
###############################################################################


# 01. Source data and functions ###############################################
source("R/data-prep.R")
source("R/simulation-functions.R")


# 02. Simulations and models ##################################################

# Set the seed
set.seed(1)

# For troubleshooting data handling, dofit -> FALSE only filters data
# for model 1 and 2.
cores <- parallel::detectCores()

for(i in 1:10) {


  d <- simulate_datasets(nullgenes = 7500,
                          condB_true = 1250,
                          condB_time2_true = 1250,
                          dispersion_model = trend_model_observed,
                         dataset = i)

  if (!dir.exists("data_sim/simdata/raw/")) dir.create("data_sim/simdata/raw/", recursive = TRUE)
  if (!dir.exists("data_sim/simdata/clean/")) dir.create("data_sim/simdata/clean/", recursive = TRUE)
  if (!dir.exists("data_sim/simdata/popeffect/")) dir.create("data_sim/simdata/popeffect/", recursive = TRUE)

  if (!dir.exists("data_sim/estimates")) dir.create("data_sim/estimates", recursive = TRUE)
  if (!dir.exists("data_sim/evaluations")) dir.create("data_sim/evaluations", recursive = TRUE)

  # Save simulated data for later
  saveRDS(d$simdat, file = paste0("data_sim/simdata/raw/dataset_",i,".RDS"))
  saveRDS(d$combined_data, file = paste0("data_sim/simdata/clean/clean_dataset_",i,".RDS"))
  saveRDS(d$population_effects, file = paste0("data_sim/simdata/popeffect/population_effects_",i,".RDS"))


  # Model 1 and 2 ##########################
  # Fitting naive and informed Neg-Binom model. The informed model has
  # a wider prior for the mean-dispersion fit as we use weighted estimates
  # of the loess regression.
  m1_m2_results <- m1_m2_sim(d$combined_data,
                             dataset = i,
                             dofit = TRUE,
                             weighted_loess = TRUE,
                             CORES = cores)

  saveRDS(m1_m2_results$summaries_m1, file = paste0("data_sim/estimates/m1_estimates_", i, ".RDS"))
  saveRDS(m1_m2_results$summaries_m2, file = paste0("data_sim/estimates/m2_estimates_", i, ".RDS"))

  saveRDS(m1_m2_results$evaluations_m1, file = paste0("data_sim/evaluations/m1_evaluations_", i, ".RDS"))
  saveRDS(m1_m2_results$evaluations_m2, file = paste0("data_sim/evaluations/m2_evaluations_", i, ".RDS"))


  # Model 1b and 2b ##########################
  # Fitting naive and informed Neg-Binom model. The informed model has
  # a more narrow prior for the mean-dispersion fit as we use un-weighted estimates
  # of the loess regression.
 # m1_m2_results_b <- m1_m2_sim(d$combined_data,
 #                            dataset = i,
 #                            dofit = TRUE,
 #                            weighted_loess = FALSE,
 #                            CORES = cores)
#
 # saveRDS(m1_m2_results_b$summaries_m1, file = paste0("data_sim/estimates/m1b_estimates_", i, ".RDS"))
 # saveRDS(m1_m2_results_b$summaries_m2, file = paste0("data_sim/estimates/m2b_estimates_", i, ".RDS"))
#
 # saveRDS(m1_m2_results_b$evaluations_m1, file = paste0("data_sim/evaluations/m1b_evaluations_", i, ".RDS"))
 # saveRDS(m1_m2_results_b$evaluations_m2, file = paste0("data_sim/evaluations/m2b_evaluations_", i, ".RDS"))


  # Model 3 ################################################
  # A model for transformed counts.
  m3_results <- m3_sim(d$combined_data,
                        dataset = i,
                        dofit = TRUE,
                        CORES = cores)

  saveRDS(m3_results$summaries_m3, file = paste0("data_sim/estimates/m3_estimates_", i, ".RDS"))
  saveRDS(m3_results$evaluations_m3, file = paste0("data_sim/evaluations/m3_evaluations_", i, ".RDS"))



  # Model 4 and 5 ################################################
  # This is the Poisson model with observation-level random effects
  # Model 5 is the informed model (with priors).
  m4_m5_results <- m4_m5_sim(d$combined_data,
                             dataset = i,
                             dofit = TRUE,
                             CORES = cores)


  saveRDS(m4_m5_results$summaries_m4, file = paste0("data_sim/estimates/m4_estimates_", i, ".RDS"))
  saveRDS(m4_m5_results$summaries_m5, file = paste0("data_sim/estimates/m5_estimates_", i, ".RDS"))

  saveRDS(m4_m5_results$evaluations_m4, file = paste0("data_sim/evaluations/m4_evaluations_", i, ".RDS"))
  saveRDS(m4_m5_results$evaluations_m5, file = paste0("data_sim/evaluations/m5_evaluations_", i, ".RDS"))

 print(paste0("Simulation 1:", i, " is done."))
}

