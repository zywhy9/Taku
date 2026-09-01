# Calculate the weekly standard deviation of the combined lake-stock proportion.
#
# The abundance target depends on stocks 2 and 3 through P_lake = P_2 + P_3.
# Consequently,
#
#   Var(P_lake) = Var(P_2) + Var(P_3) + 2 Cov(P_2, P_3).
#
# Under a Dirichlet distribution the covariance is negative. The distinction
# between sample-level and population-level uncertainty is essential here. Let
#
#   beta_hat   = 1 / (lambda + 1)
#   beta_tilde = 1 / n + (1 - 1 / n) * beta_hat
#
# denote the sample- and population-level Dirichlet variance coefficients. Their
# ratio is beta_tilde / beta_hat = 1 + lambda / n. The three methods use these
# coefficients as follows:
#
#   mom_naive: Follow the original Gazey calculation by treating the two lake
#              stocks as independent and summing their reported variances.
#   mom: Derive both marginal variances and covariance from the same
#        population-level coefficient beta_tilde.
#   mom_alt: Reconstruct the sample-level lake variance using beta_hat and
#            the reported marginals, then inflate that variance by
#            beta_tilde / beta_hat. Do not use beta_tilde inside the sample-
#            level reconstruction and inflate the covariance a second time.
#
# For mom_alt, the result is bounded below by the variance at correlation -1 so
# an approximate covariance cannot produce a negative total. No such bound is
# needed for mom_naive because it does not include a covariance term.
#
# Arguments:
#   data: One observed or simulated-data list containing mu, sd, n, and lam. Rows represent
#         weeks and columns represent stock categories.
#   method: One of "mom_naive", "mom", or "mom_alt".
#
# Returns:
#   A numeric vector with one standard deviation per week.
frequentist_weekly_sd <- function(data, method = c(
  "mom_naive",
  "mom",
  "mom_alt"
)) {
  method <- match.arg(method)
  mu <- as.matrix(data$mu)
  reported_sd <- as.matrix(data$sd)
  n <- as.numeric(data$n)
  lambda <- as.numeric(data$lam)

  if (nrow(mu) != length(n) || nrow(mu) != length(lambda) ||
      !identical(dim(mu), dim(reported_sd)) || ncol(mu) < 3L) {
    stop("data have incompatible dimensions")
  }

  # beta_hat describes the IA/sample-level variance. beta_tilde adds the
  # multinomial sampling contribution required at the population level.
  beta_hat <- 1 / (lambda + 1)
  beta_tilde <- 1 / n + (1 - 1 / n) * beta_hat
  inflation <- beta_tilde / beta_hat
  covariance_population <- -beta_tilde * mu[, 2] * mu[, 3]

  if (method == "mom") {
    # Here the marginal variances and covariance all come from one coherent
    # Dirichlet approximation rather than mixing reported and derived moments.
    marginal_variance <- (
      mu[, 2] * (1 - mu[, 2]) + mu[, 3] * (1 - mu[, 3])
    ) * beta_tilde
    return(sqrt(pmax(0, marginal_variance + 2 * covariance_population)))
  }

  sum_reported_variance <- reported_sd[, 2]^2 + reported_sd[, 3]^2

  if (method == "mom_naive") {
    # Gazey's original approximation treats stock-specific estimates as
    # independent. The reported stock variances are therefore added directly,
    # with no Dirichlet covariance correction and no correlation bound.
    return(sqrt(sum_reported_variance))
  }

  # The most negative admissible covariance is -sd_2 * sd_3. This supplies a
  # lower bound when the fitted Dirichlet covariance is more negative than the
  # reported marginal standard deviations permit.
  variance_lower_bound <- (reported_sd[, 2] - reported_sd[, 3])^2

  # First reconstruct s_lake^2 on the SAMPLE scale. The covariance here must
  # use beta_hat, not beta_tilde. The lower bound is also on the sample scale.
  covariance_sample <- -beta_hat * mu[, 2] * mu[, 3]
  sample_lake_variance <- pmax(
    sum_reported_variance + 2 * covariance_sample,
    variance_lower_bound
  )

  # Then inflate the VARIANCE to the population scale and take its square root.
  # Equivalently this is sqrt(max(f * sum_reported_variance +
  # 2 * covariance_population, f * variance_lower_bound)), because
  # f * beta_hat = beta_tilde. It is not f times the sample-level SD.
  sqrt(inflation * sample_lake_variance)
}

# Estimate total abundance and its frequentist standard error for one
# observed or simulated dataset.
#
# Let D be the weighted mean proportion belonging to stocks 2 and 3. The point
# estimate is N_hat = Nw / D. A first-order delta-method approximation propagates
# the independent weekly composition variances through this ratio.
#
# Arguments:
#   data: One observed or simulated-data list containing mu, sd, n, lam, w, and Nw.
#   method: Variance method passed to frequentist_weekly_sd().
#
# Returns:
#   A list with scalar elements estimate and standard_error.
estimate_frequentist_total <- function(data, method = c(
  "mom_naive",
  "mom",
  "mom_alt"
)) {
  method <- match.arg(method)
  lake_proportion <- rowSums(as.matrix(data$mu)[, 2:3, drop = FALSE])
  # Normalize here as well as during input preparation because externally
  # generated simulation files may store unnormalized weights.
  weights <- as.numeric(data$w)
  weights <- weights / sum(weights)
  denominator <- sum(lake_proportion * weights)

  if (!is.finite(denominator) || denominator <= 0) {
    stop("weighted lake-stock proportion must be positive")
  }

  estimate <- as.numeric(data$Nw) / denominator
  weekly_sd <- frequentist_weekly_sd(data, method = method)
  # The derivative of Nw / D with respect to each weekly lake proportion is
  # -Nw * weight / D^2 = -estimate * weight / D. Squaring and summing yields the
  # variance expression below.
  variance <- sum((weights * weekly_sd)^2) * (estimate / denominator)^2

  list(estimate = estimate, standard_error = sqrt(variance))
}

# Apply all three MoM methods to an empirical Taku input file.
#
# The input loader converts the
# reported percentages to proportions, normalizes the run weights, and fits
# the week-specific precisions using all four reporting units. The resulting
# data use exactly the same estimator functions as the simulation analysis.
# In particular, this does NOT generate synthetic data or use a simulated Nw.
#
# Arguments:
#   input_path: Private .RData file accepted by load_taku_inputs().
#   enumerated_total: Observed lake-type count (M in the paper). It is treated
#                     as fixed; uncertainty in weir counts is not modeled.
#   confidence_level: Two-sided normal confidence level, default 0.95.
#
# Returns: Three rows, with a common abundance estimate, method-specific
# standard errors, and normal-approximation confidence limits. Weeks are
# treated as independent in the delta-method variance. No MCMC or RNG is used.
estimate_taku_application <- function(
    input_path, enumerated_total, confidence_level = 0.95) {
  if (!is.numeric(enumerated_total) || length(enumerated_total) != 1L ||
      !is.finite(enumerated_total) || enumerated_total <= 0) {
    stop("enumerated_total must be one positive finite number")
  }
  if (!is.numeric(confidence_level) || length(confidence_level) != 1L ||
      !is.finite(confidence_level) || confidence_level <= 0 ||
      confidence_level >= 1) {
    stop("confidence_level must lie strictly between zero and one")
  }
  inputs <- load_taku_inputs(input_path)
  if (inputs$nstock != 4L) stop("Taku input must have four reporting units")
  data <- list(mu = inputs$p, sd = inputs$sd, n = inputs$n,
               lam = inputs$lambda, w = inputs$w, Nw = enumerated_total)
  methods <- c(mom = "MoM", mom_alt = "MoM(Alt)", mom_naive = "MoM(Naive)")
  critical <- stats::qnorm((1 + confidence_level) / 2)
  rows <- lapply(names(methods), function(method) {
    fit <- estimate_frequentist_total(data, method)
    data.frame(method = method, display_name = unname(methods[method]),
      estimate = fit$estimate, standard_error = fit$standard_error,
      lower = fit$estimate - critical * fit$standard_error,
      upper = fit$estimate + critical * fit$standard_error,
      confidence_level = confidence_level, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# Evaluate the frequentist estimators over a directory of simulation replicates.
#
# Arguments:
#   data_path: Directory containing data0000001.rds-style files written by the
#              simulation generator.
#   true_total: True abundance used to calculate bias, RMSE, and coverage.
#   methods: Variance methods to evaluate. Defaults to all three methods.
#
# Returns:
#   A data frame with one row per method and the following Monte Carlo summaries:
#     relative_bias: mean(estimate - truth) / truth;
#     relative_rmse: root mean squared error / truth;
#     coverage: proportion of normal 95% intervals containing the truth;
#     interval_length: mean width of those intervals;
#     simulations: number of successfully read replicate files.
summarize_frequentist_simulations <- function(
    data_path,
    true_total = 60000,
    methods = c("mom_naive", "mom", "mom_alt")) {
  if (!dir.exists(data_path)) {
    stop("simulation data directory does not exist: ", data_path)
  }

  data_files <- list.files(
    data_path,
    pattern = "^data[0-9]{7}[.]rds$",
    full.names = TRUE
  )
  if (!length(data_files)) {
    stop("no data0000001.rds-style files found in: ", data_path)
  }

  # Read each simulation object once. File I/O is common to every frequentist
  # estimator and can dominate these closed-form calculations, especially for
  # whichever method happens to run first on a cold filesystem cache. Excluding
  # shared loading time makes the per-method computation timings comparable.
  loading_started_at <- proc.time()[["elapsed"]]
  simulation_data <- lapply(data_files, readRDS)
  data_loading_seconds <- proc.time()[["elapsed"]] - loading_started_at

  method_results <- lapply(methods, function(method) {
    started_at <- proc.time()[["elapsed"]]
    estimates <- numeric(length(data_files))
    standard_errors <- numeric(length(data_files))

    for (index in seq_along(data_files)) {
      result <- estimate_frequentist_total(simulation_data[[index]], method)
      estimates[index] <- result$estimate
      standard_errors[index] <- result$standard_error
    }

    # The frequentist workflow reports symmetric normal-approximation
    # intervals. qnorm(0.975) is the two-sided 95% critical value.
    lower <- estimates - stats::qnorm(0.975) * standard_errors
    upper <- estimates + stats::qnorm(0.975) * standard_errors
    elapsed_seconds <- proc.time()[["elapsed"]] - started_at

    data.frame(
      scenario_id = paste0("frequentist_", method),
      display_name = switch(
        method,
        mom_naive = "MoM(Naive)",
        mom = "MoM",
        mom_alt = "MoM(Alt)",
        method
      ),
      relative_bias = mean(estimates - true_total) / true_total,
      relative_rmse = sqrt(mean((estimates - true_total)^2)) / true_total,
      coverage = mean(lower <= true_total & upper >= true_total),
      interval_length = mean(upper - lower),
      simulations = length(data_files),
      data_loading_seconds = data_loading_seconds,
      computation_seconds = elapsed_seconds,
      seconds_per_dataset = elapsed_seconds / length(data_files),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, method_results)
}
