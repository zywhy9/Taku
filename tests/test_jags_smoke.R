# Compile each model/prior combination on one synthetic dataset. This verifies
# software and JAGS wiring only; the short chains are not a convergence study.
source("R/taku.R")

required <- c("simanalyse", "mcmcr", "rjags")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("missing package(s): ", paste(missing, collapse = ", "))
if (as.character(rjags::jags.version()) != "4.3.2") {
  stop("JAGS 4.3.2 is required")
}

inputs <- list(
  n = rep(100L, 12L),
  w = rep(1 / 12, 12L),
  p = matrix(rep(c(.2, .3, .4, .1), 12L), 12L, 4L, byrow = TRUE),
  lambda = rep(25, 12L),
  nweek = 12L,
  nstock = 4L
)
set.seed(2468)
data <- simulate_taku_dataset(inputs, true_total = 60000)
jags_data <- data[c("mu", "lam", "n", "nweek", "nstock", "Nw", "w")]

for (index in seq_len(nrow(taku_scenarios()))) {
  scenario <- taku_scenarios()[index, ]
  code <- paste(
    "model {",
    analysis_model_code(scenario$observation_model),
    prior_model_code(scenario$prior),
    "}",
    sep = "\n"
  )
  connection <- textConnection(code)
  model <- rjags::jags.model(
    connection,
    data = jags_data,
    inits = scenario_initial_values(
      scenario$prior, data$nweek, data$nstock, n_chains = 2L
    ),
    n.chains = 2L,
    n.adapt = 1000L,
    quiet = TRUE
  )
  close(connection)
  samples <- rjags::coda.samples(
    model, variable.names = "N", n.iter = 100L, progress.bar = "none"
  )
  stopifnot(all(is.finite(unlist(samples))))
  cat("PASS: ", scenario$scenario_id, "\n", sep = "")
}

# Exercise the public simanalyse wrapper once, including the compact p means.
settings <- taku_paper_settings(
  n_datasets = 1L,
  n_chains = 2L,
  n_adapt = 1000L,
  n_save = 100L,
  rhat = 10,
  ess = 1,
  max_iter = 100L,
  max_time_minutes = 1,
  save_p_means = TRUE
)
task <- taku_baseline_tasks(1L, 1L, 1234L)
task <- task[task$scenario_id == "mmd_dirichlet", , drop = FALSE]
fit <- taku_fit_baseline(data, task, settings)
stopifnot(
  is.matrix(fit$p_mean),
  identical(dim(fit$p_mean), c(12L, 4L)),
  max(abs(rowSums(fit$p_mean) - 1)) < 1e-10,
  fit$row$converged,
  fit$row$batches >= 1L
)
cat("PASS: public simanalyse wrapper and compact p means\n")

cat("JAGS smoke checks passed.\n")
