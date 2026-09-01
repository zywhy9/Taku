# Fast regression checks. No empirical data and no MCMC are required.
source("R/taku.R")

stopifnot(
  identical(
    taku_scenarios()$scenario_id,
    c("rdm_ar", "rdm_dirichlet", "mmd_ar", "mmd_dirichlet")
  ),
  identical(taku_lake_categories(), c(2L, 3L))
)

paper <- taku_paper_settings()
stopifnot(
  paper$n_datasets == 1000L,
  paper$seed == 1234L,
  paper$n_chains == 3L,
  paper$n_adapt == 100000L,
  paper$n_save == 50000L,
  paper$rhat == 1.10,
  paper$ess == 100,
  paper$max_iter == 1000000000,
  paper$max_time_minutes == 1320,
  isTRUE(paper$save_p_means)
)

local({
  root <- tempfile("taku-core-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)

  n <- rep(100L, 12L)
  w <- seq_len(12L)
  p <- matrix(rep(c(.20, .30, .40, .10), 12L), nrow = 12L, byrow = TRUE)
  region1mean <- 100 * p
  region1sd <- matrix(2, 12L, 4L)
  input <- file.path(root, "input.RData")
  save(n, w, region1mean, region1sd, file = input)

  loaded <- load_taku_inputs(input)
  stopifnot(
    identical(dim(loaded$p), c(12L, 4L)),
    abs(sum(loaded$w) - 1) < 1e-12,
    all(loaded$lambda > 0)
  )

  generated <- generate_taku_datasets(
    input, file.path(root, "study"), nsims = 3L, seed = 1234L
  )
  checked <- validate_taku_simulation_run(
    loaded, file.path(root, "study"), true_total = 60000
  )
  stopifnot(
    generated$manifest$seed == 1234L,
    checked$manifest$nsims == 3L,
    length(list.files(
      generated$simulation_path, pattern = "^data[0-9]{7}[.]rds$"
    )) == 3L
  )

  data <- readRDS(file.path(generated$simulation_path, "data0000001.rds"))
  estimates <- lapply(c("mom", "mom_alt", "mom_naive"), function(method) {
    estimate_frequentist_total(data, method)
  })
  stopifnot(all(vapply(estimates, function(x) {
    is.finite(x$estimate) && is.finite(x$standard_error) &&
      x$standard_error >= 0
  }, logical(1))))

  # Exercise the automatic paper-table collector with compact fake MCMC rows.
  study <- file.path(root, "study")
  settings <- taku_paper_settings(n_datasets = 3L, n_save = 100L)
  taku_write_or_check_settings(study, settings)
  for (scenario in taku_scenarios()$scenario_id) {
    for (dataset_id in 1:3) {
      row <- data.frame(
        dataset_id = dataset_id,
        scenario_id = scenario,
        estimate = 60000 + dataset_id,
        lower = 56000,
        upper = 64000,
        converged = TRUE,
        sampling_seconds = dataset_id
      )
      taku_atomic_rds(
        list(row = row),
        taku_record_path(study, scenario, dataset_id)
      )
    }
  }
  summary <- taku_summarize_simulation(study)
  stopifnot(
    nrow(summary) == 7L,
    identical(
      names(summary),
      c(
        "scenario_id", "display_name", "method_family",
        "observation_model", "prior", "relative_bias", "relative_rmse",
        "coverage", "interval_length", "time_seconds_per_dataset",
        "simulations"
      )
    ),
    all(summary$simulations == 3L),
    file.exists(file.path(study, "results", "simulation_convergence.csv"))
  )
})

# Protect the corrected MoM(Alt) scale conversion and Gazey-style naive rule.
data <- list(
  mu = rbind(c(.2, .3, .4, .1), c(.1, .2, .5, .2)),
  sd = rbind(c(.03, .08, .09, .03), c(.03, .04, .05, .03)),
  n = c(50, 100), lam = c(20, 30)
)
beta_hat <- 1 / (data$lam + 1)
beta_tilde <- 1 / data$n + (1 - 1 / data$n) * beta_hat
inflation <- beta_tilde / beta_hat
reported <- data$sd[, 2]^2 + data$sd[, 3]^2
sample_lake <- pmax(
  reported - 2 * beta_hat * data$mu[, 2] * data$mu[, 3],
  (data$sd[, 2] - data$sd[, 3])^2
)
stopifnot(
  isTRUE(all.equal(frequentist_weekly_sd(data, "mom_naive")^2, reported)),
  isTRUE(all.equal(
    frequentist_weekly_sd(data, "mom_alt")^2,
    inflation * sample_lake
  ))
)

cat("Core checks passed.\n")
