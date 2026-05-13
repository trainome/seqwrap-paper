#' Run model 1: Uninformed NegBinom model on the Pillon data.
#'
#' This functions runs the model and store results in derived_data.
#'
#' @param CORES Integer. Number of cores to pass to `seqwrap`.
#' @param test Logical. For testing purposes, if TRUE the function will output
#' a list from `seqwrap_summarise(m1_results)` on a subset (1:24).
#' @export
run_model1 <- function(CORES, test = FALSE) {
  if (!dir.exists("analysis/data/derived_data/"))
    dir.create("analysis/data/derived_data/")

  dat <- seqwrappaper::pillon_counts
  all(dat$metadata$seq_sample_id == colnames(dat$countdata[, -1]))

  # Combine all gene counts after filtering
  keep <- edgeR::filterByExpr(
    dat$countdata[, -1],
    min.count = 10,
    min.total.count = 15,
    large.n = 10,
    min.prop = 0.7,
    group = paste(dat$metadata$group, dat$metadata$time)
  )

  countdat <- dat$countdata[keep, ]

  # Use EdgeR to calculate the TMM
  y <- edgeR::DGEList(countdat[, -1])
  y <- edgeR::calcNormFactors(y)

  # Store library sizes
  libsize <- y$samples |>
    rownames_to_column(var = "seq_sample_id") |>
    select(-group)

  # Combine all meta data
  metadat <- dat$metadata |>
    inner_join(libsize, by = "seq_sample_id") |>
    mutate(
      group = factor(group, levels = c("NGT", "T2D")),
      time = factor(time, levels = c("basal", "post", "rec")),
      efflibsize = (lib.size * norm.factors) /
        median(lib.size * norm.factors),
      ln_efflibsize = log(efflibsize)
    )

  # 03. A preliminary model #####################################################

  # A preliminary model is fitted using a conditional NB distribution.
  # The purpose of the conditional model is to estimate distributions
  # of parameter estimates. We will use these for simulations

  # A summary function to return the dispersion parameter with SE on the log scale
  # mean(predict(x, type = "link)) will give us the predicted log counts.
  # We will put this in the eval fun to also get estimates of the parameters in
  # the generic summary function.
  sigma_summary <- function(x) {
    out <- data.frame(
      dispersion = data.frame(summary(x$sdr))["betadisp", 1],
      dispersion.se = data.frame(summary(x$sdr))["betadisp", 2],
      log_mu = mean(predict(x, type = "link"))
    )
    return(out)
  }

  m1 <- seqwrap_compose(
    data = countdat,
    metadata = metadat,
    samplename = "seq_sample_id",
    modelfun = glmmTMB,
    eval_fun = sigma_summary,
    arguments = list(
      formula = y ~ time * group + offset(ln_efflibsize) + (1 | id),
      family = glmmTMB::nbinom2
    )
  )

  if (test) {
    m1_results <- seqwrap(
      m1,
      return_models = FALSE,
      subset = 1:24,
      cores = CORES
    )

    return(seqwrap_summarise(m1_results))
  }

  if (!file.exists("analysis/data/derived_data/m1_results.RDS")) {
    m1_results <- seqwrap(
      m1,
      return_models = FALSE,
      cores = CORES
    )

    saveRDS(m1_results, "analysis/data/derived_data/m1_results.RDS")
  } else {
    m1_results <- readRDS("analysis/data/derived_data/m1_results.RDS")
  }

  return(m1_results)
}


#' Run model 2: Informed NegBinom model on the Pillon data.
#'
#' @param CORES Integer. Number of cores to pass to `seqwrap`.
#' @param test Logical. For testing purposes, if TRUE the function will output
#' a list from `seqwrap_summarise(m2_results)` on a subset (1:24).
#' @export
run_model2 <- function(CORES, test = FALSE) {
  if (!dir.exists("analysis/data/derived_data/"))
    dir.create("analysis/data/derived_data/")

  # Run only if file does not exist
  if (!file.exists("analysis/data/derived_data/m1_results.RDS")) {
    run_model1(CORES = CORES)
  }

  dat <- seqwrappaper::pillon_counts
  all(dat$metadata$seq_sample_id == colnames(dat$countdata[, -1]))

  # Combine all gene counts after filtering
  keep <- edgeR::filterByExpr(
    dat$countdata[, -1],
    min.count = 10,
    min.total.count = 15,
    large.n = 10,
    min.prop = 0.7,
    group = paste(dat$metadata$group, dat$metadata$time)
  )

  countdat <- dat$countdata[keep, ]

  # Use EdgeR to calculate the TMM
  y <- edgeR::DGEList(countdat[, -1])
  y <- edgeR::calcNormFactors(y)

  # Store library sizes
  libsize <- y$samples |>
    rownames_to_column(var = "seq_sample_id") |>
    select(-group)

  # Combine all meta data
  metadat <- dat$metadata |>
    inner_join(libsize, by = "seq_sample_id") |>
    mutate(
      group = factor(group, levels = c("NGT", "T2D")),
      time = factor(time, levels = c("basal", "post", "rec")),
      efflibsize = (lib.size * norm.factors) /
        median(lib.size * norm.factors),
      ln_efflibsize = log(efflibsize)
    )

  # Load the available model
  m1_results <- readRDS("analysis/data/derived_data/m1_results.RDS")

  m1_sums <- seqwrap_summarise(m1_results)

  ## Model 2 ##
  # get successful targets
  targets <- m1_sums$evaluations |>
    filter(!is.na(dispersion.se)) |>
    distinct(target) |>
    pull(target)

  # 05. Prepare priors for model 2 #############################################

  ## Fitting a model for the mean-dispersion relationship

  # Fit a trend to the dispersion data from m1
  # save the data in a convenient format. Using the log_mu_raw (average
  # observed counts) allows for a prior on dispersion also for genes with
  # unsuccessful fits in model 1. First we gather all dispersion data.

  dispersion_dat <- m1_sums$evaluation |>
    filter(target %in% targets)

  # Calculate the raw observed counts from the data
  raw_log_counts <- data.frame(
    target = as.character(countdat[, 1]),
    log_mu_raw = log(rowMeans(countdat[, -1]))
  )

  # Adding log raw counts to the dispersion df for modeling.
  #
  dispersion_dat <- dispersion_dat |>
    inner_join(raw_log_counts)

  # Fit a loess model, using log_mu_raw as the predictor
  weighted_loess <- TRUE

  if (weighted_loess) {
    trend_model <- loess(
      dispersion ~ log_mu_raw,
      data = dispersion_dat,
      span = 0.7,
      weights = 1 / (dispersion.se^2)
    )
  } else {
    trend_model <- loess(
      dispersion ~ log_mu_raw,
      data = dispersion_dat,
      span = 0.7
    )
  }

  # Display the model estimate mean-dispersion parameter
  ## Preliminary plot ##
  #  dispersion_dat |>
  #    mutate(pred_dispersion = predict(trend_model)) |>
  #    ggplot(aes(log_mu_raw, dispersion)) +
  #
  #    geom_ribbon(aes(ymin =  pred_dispersion - trend_model$s,
  #                    ymax = pred_dispersion + trend_model$s),
  #                alpha = 0.2) +
  #
  #    geom_point() +
  #    geom_line(aes(log_mu_raw, pred_dispersion), color = "red",
  #              linewidth = 2)
  #

  # Predict dispersion for each gene based on log raw counts
  # and combine into a prior.

  dispersion_prior <- data.frame(
    gene = countdat[, 1],
    pred = round(
      predict(
        trend_model,
        newdata = data.frame(
          log_mu_raw = log(
            rowMeans(
              countdat[, -1]
            )
          )
        )
      ),
      3
    ),
    s = round(
      trend_model$s,
      3
    )
  ) |>
    mutate(prior = paste0("normal(", pred, ",", s, ")"))

  # Gather all distributions of estimates that can be used as prior information
  # in subsequent model. Fixed effects are centered on 0. We will use the SD
  # for priors

  ## Preliminary plot ##
  #  ms1_sum$summaries |>
  #  filter(target %in% targets) |>
  #  select(target, term, estimate) |>
  #  ggplot(aes(estimate)) + geom_density() +
  #  facet_wrap(~ term, scales = "free")

  estimate_distributions <- m1_sums$summaries |>
    filter(target %in% targets) |>
    select(target, term, estimate) |>

    summarise(.by = term, m = mean(estimate), s = sd(estimate)) |>
    filter(!term %in% c("(Intercept)", "sd__(Intercept)"))

  # Extract the random effects distribution to fit a gamma distribution
  random_sd_estimate <- m1_sums$summaries |>
    filter(target %in% targets, term == "sd__(Intercept)") |>
    select(target, term, estimate) |>
    pull(estimate)

  # The gamma distribution is parameterized using a shape and a rate
  # parameter. It looks like this prior will lead to a push towards 0,
  # consider adding a constant to push away from zero...
  # TODO this may needs testing.
  mean_sd <- mean(random_sd_estimate)
  var_sd <- var(random_sd_estimate)
  shape_param <- 2 # mean_sd^2 / var_sd

  # Here we prepare priors for the fixed effects, in this version all fixed
  # effects, except the intercept, will have regularizing priors corresponding to
  # the distributions of effects seen in the naive models.
  Priors_df <- bind_rows(
    data.frame(
      prior = paste0("normal(0,", round(estimate_distributions$s, 2), ")"),
      class = rep("fixef", 5),
      coef = estimate_distributions$term
    ),
    data.frame(
      prior = paste0(
        "gamma(",
        round(mean_sd, 2),
        ",",
        2,
        ")"
      ),
      class = "ranef",
      coef = "id"
    )
  )

  # We want to use the mean-dispersion relationship to add a prior for the
  # dispersion parameter. This means that we need a gene specific prior

  Priors_list <- list()
  for (j in 1:nrow(countdat)) {
    df <- bind_rows(
      Priors_df,
      data.frame(
        prior = dispersion_prior[dispersion_prior$gene == countdat[j, 1], 4],
        class = "fixef_disp",
        coef = "1"
      )
    )

    Priors_list[[j]] <- df
  }

  # 06. Fit model 2 ##############################################################

  m2 <- seqwrap_compose(
    data = countdat,
    metadata = metadat,
    samplename = "seq_sample_id",
    modelfun = glmmTMB::glmmTMB,
    eval_fun = sigma_summary2,
    targetdata = Priors_list,
    arguments = alist(
      formula = y ~ time * group + offset(ln_efflibsize) + (1 | id),
      family = glmmTMB::nbinom2,
      priors = data.frame(
        prior = prior,
        class = class,
        coef = coef
      )
    )
  )

  if (test) {
    m2_results <- seqwrap(
      m2,
      return_models = FALSE,
      subset = 1:24,
      cores = CORES
    )

    return(seqwrap_summarise(m2_results))
  }

  if (!file.exists("analysis/data/derived_data/m2_results.RDS")) {
    m2_results <- seqwrap(
      m2,
      return_models = FALSE,
      verbose = FALSE,
      cores = CORES
    )

    saveRDS(m2_results, "analysis/data/derived_data/m2_results.RDS")
  } else {
    m2_results <- readRDS("analysis/data/derived_data/m2_results.RDS")
  }

  return(m2_results)
}


#' Run model 3: VST-transformed linear mixed model on the Pillon data.
#'
#' @param CORES Integer. Number of cores to pass to `seqwrap`.
#' @param test Logical. For testing purposes, if TRUE the function will output
#' a list from `seqwrap_summarise(m3_results)` on a subset (1:24).
#' @export
run_model3 <- function(CORES, test = FALSE) {
  if (!dir.exists("analysis/data/derived_data/"))
    dir.create("analysis/data/derived_data/")

  dat <- seqwrappaper::pillon_counts
  all(dat$metadata$seq_sample_id == colnames(dat$countdata[, -1]))

  # Combine all gene counts after filtering
  keep <- edgeR::filterByExpr(
    dat$countdata[, -1],
    min.count = 10,
    min.total.count = 15,
    large.n = 10,
    min.prop = 0.7,
    group = paste(dat$metadata$group, dat$metadata$time)
  )

  countdat <- dat$countdata[keep, ]

  # Use EdgeR to calculate the TMM
  y <- edgeR::DGEList(countdat[, -1])
  y <- edgeR::calcNormFactors(y)

  # Store library sizes
  libsize <- y$samples |>
    rownames_to_column(var = "seq_sample_id") |>
    select(-group)

  # Combine all meta data
  metadat <- dat$metadata |>
    inner_join(libsize, by = "seq_sample_id") |>
    mutate(
      group = factor(group, levels = c("NGT", "T2D")),
      time = factor(time, levels = c("basal", "post", "rec")),
      efflibsize = (lib.size * norm.factors) /
        median(lib.size * norm.factors),
      ln_efflibsize = log(efflibsize)
    )

  dds <- DESeqDataSetFromMatrix(
    countData = countdat[, -1],
    colData = metadat,
    design = ~ time * group
  )

  dds <- DESeq(dds, quiet = TRUE)
  vst_mat <- assay(varianceStabilizingTransformation(
    dds,
    blind = FALSE,
    fitType = "parametric"
  ))
  vst_dat <- cbind(data.frame(gene = countdat[, 1], as.data.frame(vst_mat)))

  m3 <- seqwrap_compose(
    data = vst_dat,
    metadata = metadat,
    samplename = "seq_sample_id",
    modelfun = lmerTest::lmer,
    eval_fun = lmer_summary,
    arguments = list(
      formula = y ~ time * group + (1 | id)
    )
  )

  if (test) {
    m3_results <- seqwrap(
      m3,
      return_models = FALSE,
      subset = 1:24,
      verbose = FALSE,
      cores = CORES
    )

    return(seqwrap_summarise(m3_results))
  }

  if (!file.exists("analysis/data/derived_data/m3_results.RDS")) {
    m3_results <- seqwrap(
      m3,
      return_models = FALSE,
      verbose = FALSE,
      cores = CORES
    )

    saveRDS(m3_results, "analysis/data/derived_data/m3_results.RDS")
  } else {
    m3_results <- readRDS("analysis/data/derived_data/m3_results.RDS")
  }

  return(m3_results)
}


#' Run model 4: Poisson model with observation-level random effect on the Pillon data.
#'
#' @param CORES Integer. Number of cores to pass to `seqwrap`.
#' @param test Logical. For testing purposes, if TRUE the function will output
#' a list from `seqwrap_summarise(m4_results)` on a subset (1:24).
#' @export
run_model4 <- function(CORES, test = FALSE) {
  if (!dir.exists("analysis/data/derived_data/"))
    dir.create("analysis/data/derived_data/")

  dat <- seqwrappaper::pillon_counts
  all(dat$metadata$seq_sample_id == colnames(dat$countdata[, -1]))

  # Combine all gene counts after filtering
  keep <- edgeR::filterByExpr(
    dat$countdata[, -1],
    min.count = 10,
    min.total.count = 15,
    large.n = 10,
    min.prop = 0.7,
    group = paste(dat$metadata$group, dat$metadata$time)
  )

  countdat <- dat$countdata[keep, ]

  # Use EdgeR to calculate the TMM
  y <- edgeR::DGEList(countdat[, -1])
  y <- edgeR::calcNormFactors(y)

  # Store library sizes
  libsize <- y$samples |>
    rownames_to_column(var = "seq_sample_id") |>
    select(-group)

  # Combine all meta data
  metadat <- dat$metadata |>
    inner_join(libsize, by = "seq_sample_id") |>
    mutate(
      group = factor(group, levels = c("NGT", "T2D")),
      time = factor(time, levels = c("basal", "post", "rec")),
      efflibsize = (lib.size * norm.factors) /
        median(lib.size * norm.factors),
      ln_efflibsize = log(efflibsize)
    )

  m4 <- seqwrap_compose(
    data = countdat, # These are the filtered counts
    metadata = metadat,
    samplename = "seq_sample_id",
    modelfun = glmmTMB::glmmTMB,
    eval_fun = sigma_summary2,
    targetdata = NULL,
    arguments = list(
      formula = y ~
        time * group + offset(ln_efflibsize) + (1 | id) + (1 | seq_sample_id),
      family = stats::poisson
    )
  )

  if (test) {
    m4_results <- seqwrap(
      m4,
      return_models = FALSE,
      subset = 1:24,
      verbose = FALSE,
      cores = CORES
    )

    return(seqwrap_summarise(m4_results))
  }

  if (!file.exists("analysis/data/derived_data/m4_results.RDS")) {
    m4_results <- seqwrap(
      m4,
      return_models = FALSE,
      verbose = FALSE,
      cores = CORES
    )

    saveRDS(m4_results, "analysis/data/derived_data/m4_results.RDS")
  } else {
    m4_results <- readRDS("analysis/data/derived_data/m4_results.RDS")
  }

  return(m4_results)
}


#' Run model 5: Informed Poisson model with observation-level RE on the Pillon data.
#'
#' @param CORES Integer. Number of cores to pass to `seqwrap`.
#' @param test Logical. For testing purposes, if TRUE the function will output
#' a list from `seqwrap_summarise(m5_results)` on a subset (1:24).
#' @export
run_model5 <- function(CORES, test = FALSE) {
  if (!dir.exists("analysis/data/derived_data/"))
    dir.create("analysis/data/derived_data/")

  # Run only if file does not exist
  if (!file.exists("analysis/data/derived_data/m4_results.RDS")) {
    run_model4(CORES = CORES)
  }

  dat <- seqwrappaper::pillon_counts
  all(dat$metadata$seq_sample_id == colnames(dat$countdata[, -1]))

  # Combine all gene counts after filtering
  keep <- edgeR::filterByExpr(
    dat$countdata[, -1],
    min.count = 10,
    min.total.count = 15,
    large.n = 10,
    min.prop = 0.7,
    group = paste(dat$metadata$group, dat$metadata$time)
  )

  countdat <- dat$countdata[keep, ]

  # Use EdgeR to calculate the TMM
  y <- edgeR::DGEList(countdat[, -1])
  y <- edgeR::calcNormFactors(y)

  # Store library sizes
  libsize <- y$samples |>
    rownames_to_column(var = "seq_sample_id") |>
    select(-group)

  # Combine all meta data
  metadat <- dat$metadata |>
    inner_join(libsize, by = "seq_sample_id") |>
    mutate(
      group = factor(group, levels = c("NGT", "T2D")),
      time = factor(time, levels = c("basal", "post", "rec")),
      efflibsize = (lib.size * norm.factors) /
        median(lib.size * norm.factors),
      ln_efflibsize = log(efflibsize)
    )

  m4_results <- readRDS("analysis/data/derived_data/m4_results.RDS")
  m4_sum <- seqwrap_summarise(m4_results, verbose = FALSE)

  targets <- m4_sum$evaluations |>
    filter(convergence == 0) |>
    distinct(target) |>
    pull(target)

  ## Preliminary plot ##
  #  ms1_sum$summaries |>
  #  filter(target %in% targets) |>
  #  select(target, term, estimate) |>
  #  ggplot(aes(estimate)) + geom_density() +
  #  facet_wrap(~ term, scales = "free")

  estimate_distributions <- m4_sum$summaries |>
    filter(target %in% targets) |>
    select(target, term, group, estimate) |>

    summarise(.by = c(term, group), m = mean(estimate), s = sd(estimate)) |>
    filter(!term %in% c("(Intercept)", "sd__(Intercept)"))

  # Extract the random effects distribution to fit a gamma distribution
  # on the participant level intercepts
  random_sd_estimate <- m4_sum$summaries |>
    filter(target %in% targets, term == "sd__(Intercept)", group == "id") |>
    select(target, term, estimate) |>
    pull(estimate)

  random_sd_estimate_obs <- m4_sum$summaries |>
    filter(
      target %in% targets,
      term == "sd__(Intercept)",
      group == "seq_sample_id"
    ) |>
    select(target, term, estimate) |>
    pull(estimate)

  # The gamma distribution is parameterized using a shape and a rate
  # parameter. It looks like this prior will lead to a push towards 0,
  # consider adding a constant to push away from zero...
  # TODO this may needs testing.
  mean_sd <- mean(random_sd_estimate)
  var_sd <- var(random_sd_estimate)
  shape_param <- 2 # mean_sd^2 / var_sd
  # Observation level
  mean_sd_obs <- mean(random_sd_estimate_obs)
  var_sd_obs <- var(random_sd_estimate_obs)
  shape_param_obs <- 2 # mean_sd^2 / var_sd

  # Here we prepare priors for the fixed effects, in this version all fixed
  # effects, except the intercept, will have regularizing priors corresponding to
  # the distributions of effects seen in the naive models.

  # The priors for the random effects are added, glmmTMB accepts coefficient
  # index.
  Priors_df <- bind_rows(
    data.frame(
      prior = paste0("normal(0,", round(estimate_distributions$s, 2), ")"),
      class = rep("fixef", 5),
      coef = estimate_distributions$term
    ),
    data.frame(
      prior = paste0(
        "gamma(",
        c(round(mean_sd, 2), round(mean_sd_obs, 2)),
        ",",
        2,
        ")"
      ),
      class = "ranef",
      coef = c("1", "2")
    )
  )

  # We want to use the mean-dispersion relationship to add a prior for the
  # dispersion parameter. This means that we need a gene specific prior

  Priors_list <- list()
  for (j in 1:nrow(countdat)) {
    Priors_list[[j]] <- Priors_df
  }

  # seqwrap accepts target-wise data frames as a list,
  # this makes it easier to specify target-specific priors.

  # Here we specify priors based the results from ms1

  m5 <- seqwrap_compose(
    data = countdat,
    metadata = metadat,
    samplename = "seq_sample_id",
    modelfun = glmmTMB::glmmTMB,
    eval_fun = sigma_summary2,
    targetdata = Priors_list,
    arguments = alist(
      formula = y ~
        time * group + offset(ln_efflibsize) + (1 | id) + (1 | seq_sample_id),
      family = stats::poisson,
      priors = data.frame(
        prior = prior,
        class = class,
        coef = coef
      )
    )
  )

  if (test) {
    m5_results <- seqwrap(
      m5,
      return_models = FALSE,
      subset = 1:24,
      verbose = FALSE,
      cores = CORES
    )

    return(seqwrap_summarise(m5_results))
  }

  if (!file.exists("analysis/data/derived_data/m5_results.RDS")) {
    m5_results <- seqwrap(
      m5,
      return_models = FALSE,
      verbose = FALSE,
      cores = CORES
    )

    saveRDS(m5_results, "analysis/data/derived_data/m5_results.RDS")
  } else {
    m5_results <- readRDS("analysis/data/derived_data/m5_results.RDS")
  }

  return(m5_results)
}
