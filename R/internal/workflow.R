# Reader-facing workflow helpers. The three scripts in the repository root call
# these functions; readers normally do not need to source internal files.

# Scientific controls used for the reported baseline study.  Arguments make a
# deliberately small smoke test possible without editing model code.
taku_paper_settings <- function(
    n_datasets = 1000L,
    truth = 60000,
    seed = 1234L,
    n_chains = 3L,
    n_adapt = 100000L,
    n_save = 50000L,
    rhat = 1.10,
    ess = 100,
    max_iter = 1000000000,
    max_time_minutes = 1320,
    save_p_means = TRUE,
    save_full_draws = FALSE) {
  settings <- list(
    n_datasets = as.integer(n_datasets),
    truth = as.numeric(truth),
    seed = as.integer(seed),
    n_chains = as.integer(n_chains),
    n_adapt = as.integer(n_adapt),
    n_save = as.integer(n_save),
    rhat = as.numeric(rhat),
    ess = as.numeric(ess),
    max_iter = as.numeric(max_iter),
    max_time_minutes = as.numeric(max_time_minutes),
    save_p_means = isTRUE(save_p_means),
    save_full_draws = isTRUE(save_full_draws)
  )
  positive <- c(
    "n_datasets", "truth", "n_chains", "n_adapt", "n_save", "rhat",
    "ess", "max_iter", "max_time_minutes"
  )
  if (any(!vapply(settings[positive], function(x) {
    length(x) == 1L && is.finite(x) && x > 0
  }, logical(1)))) {
    stop("all simulation and MCMC controls must be positive finite scalars")
  }
  if (is.na(settings$seed) || settings$rhat < 1) {
    stop("seed must be an integer and rhat must be at least one")
  }
  settings
}

taku_settings_equal <- function(x, y) {
  identical(names(x), names(y)) && isTRUE(all.equal(x, y, tolerance = 0))
}

taku_write_or_check_settings <- function(output_path, settings) {
  path <- file.path(output_path, "settings.rds")
  if (file.exists(path)) {
    previous <- readRDS(path)
    if (!taku_settings_equal(previous, settings)) {
      stop("output directory was created with different settings: ", path)
    }
  } else {
    dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
    taku_atomic_rds(settings, path)
  }
  invisible(path)
}

taku_scenario_ids <- function(scenario = "all") {
  available <- taku_scenarios()$scenario_id
  if (identical(scenario, "all")) return(available)
  scenario <- as.character(scenario)
  unknown <- setdiff(scenario, available)
  if (length(unknown)) {
    stop("unknown scenario: ", paste(unknown, collapse = ", "))
  }
  unique(scenario)
}

taku_record_path <- function(output_path, scenario_id, dataset_id) {
  file.path(
    output_path,
    "records",
    scenario_id,
    sprintf("data%07d.rds", as.integer(dataset_id))
  )
}

# Fit one scenario to one generated dataset and save only the quantities needed
# for the paper tables.  Full posterior draws are optional and off by default.
taku_fit_simulation_dataset <- function(
    data_file, task, output_path, settings) {
  destination <- taku_record_path(
    output_path, task$scenario_id, task$dataset_id
  )
  if (file.exists(destination)) return(readRDS(destination))

  data <- readRDS(data_file)
  result <- taku_fit_baseline(data, task, settings)
  record <- list(
    row = result$row,
    diagnostics = result$diagnostics,
    timing = result$timing,
    mode = result$mode,
    seed = task$seed,
    scenario_id = task$scenario_id,
    dataset_id = task$dataset_id
  )
  if (!is.null(result$p_mean)) record$p_mean <- result$p_mean
  if (isTRUE(settings$save_full_draws)) record$draws <- result$draws
  taku_atomic_rds(record, destination)
  record
}

taku_run_one_scenario <- function(
    scenario_id, simulation_path, output_path, workers, settings) {
  scenario <- get_taku_scenario(scenario_id)
  tasks <- taku_baseline_tasks(
    n = settings$n_datasets,
    shards = 1L,
    seed_base = settings$seed
  )
  tasks <- tasks[tasks$scenario_id == scenario_id, , drop = FALSE]
  data_files <- file.path(
    simulation_path,
    sprintf("data%07d.rds", tasks$dataset_id)
  )
  if (!all(file.exists(data_files))) {
    stop("simulation files are incomplete for ", scenario_id)
  }

  cat("\n", scenario$display_name, ": ", nrow(tasks), " datasets\n", sep = "")
  fit_index <- function(index) {
    cat(sprintf("[%s] %s dataset %d/%d\n", format(Sys.time()),
                scenario_id, index, nrow(tasks)))
    taku_fit_simulation_dataset(
      data_file = data_files[index],
      task = tasks[index, , drop = FALSE],
      output_path = output_path,
      settings = settings
    )
  }

  if (workers == 1L) {
    records <- lapply(seq_len(nrow(tasks)), fit_index)
  } else {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      stop("future.apply is required when workers > 1")
    }
    previous_plan <- future::plan()
    on.exit(future::plan(previous_plan), add = TRUE)
    future::plan(future::multisession, workers = workers)
    records <- future.apply::future_lapply(
      seq_len(nrow(tasks)), fit_index, future.seed = TRUE
    )
  }
  invisible(records)
}

# Generate and validate the shared baseline datasets without starting MCMC.
# Local and cluster workflows both call this function, so the data-generating
# process, seed, settings checks, and directory layout cannot drift apart.
taku_prepare_simulation <- function(
    input_path,
    output_path,
    settings = taku_paper_settings()) {
  output_path <- normalizePath(output_path, winslash = "/", mustWork = FALSE)
  taku_write_or_check_settings(output_path, settings)

  simulation_path <- file.path(output_path, "simdata", "baseline")
  if (!file.exists(file.path(simulation_path, "manifest.rds"))) {
    generate_taku_datasets(
      input_path = input_path,
      output_path = output_path,
      nsims = settings$n_datasets,
      seed = settings$seed,
      true_total = settings$truth
    )
  }
  validated <- validate_taku_simulation_run(
    inputs = load_taku_inputs(input_path),
    output_path = output_path,
    true_total = settings$truth
  )
  if (validated$manifest$nsims != settings$n_datasets ||
      validated$manifest$seed != settings$seed) {
    stop("generated data do not match n_datasets and seed in settings.rds")
  }

  invisible(validated)
}

# Run any or all four Bayesian scenarios.  Data are generated once and reused;
# completed per-dataset records make an interrupted run safely resumable.
taku_run_simulation <- function(
    input_path,
    output_path,
    scenario = "all",
    workers = 1L,
    settings = taku_paper_settings()) {
  workers <- as.integer(workers)
  if (length(workers) != 1L || is.na(workers) || workers < 1L) {
    stop("workers must be one positive integer")
  }
  scenario_ids <- taku_scenario_ids(scenario)
  output_path <- normalizePath(output_path, winslash = "/", mustWork = FALSE)
  prepared <- taku_prepare_simulation(input_path, output_path, settings)
  simulation_path <- prepared$simulation_path

  print(settings)
  for (scenario_id in scenario_ids) {
    taku_run_one_scenario(
      scenario_id = scenario_id,
      simulation_path = simulation_path,
      output_path = output_path,
      workers = workers,
      settings = settings
    )
  }
  taku_summarize_simulation(output_path)
}

taku_read_scenario_records <- function(output_path, scenario_id, expected) {
  paths <- taku_record_path(output_path, scenario_id, seq_len(expected))
  if (!all(file.exists(paths))) return(NULL)
  lapply(paths, readRDS)
}

# Create all paper-facing tables automatically after every complete run.
taku_summarize_simulation <- function(output_path) {
  settings <- readRDS(file.path(output_path, "settings.rds"))
  scenario_table <- taku_scenarios()
  bayesian_performance <- list()
  convergence <- list()
  dataset_rows <- list()
  p_rows <- list()

  for (index in seq_len(nrow(scenario_table))) {
    scenario <- scenario_table[index, , drop = FALSE]
    records <- taku_read_scenario_records(
      output_path, scenario$scenario_id, settings$n_datasets
    )
    if (is.null(records)) next
    rows <- do.call(rbind, lapply(records, `[[`, "row"))
    performance <- taku_baseline_performance(
      rows, truth = settings$truth, expected = settings$n_datasets
    )
    bayesian_performance[[length(bayesian_performance) + 1L]] <- data.frame(
      scenario_id = scenario$scenario_id,
      display_name = scenario$display_name,
      method_family = "Bayesian",
      observation_model = toupper(scenario$observation_model),
      prior = if (scenario$prior == "ar1") "AR(1)" else "Dirichlet",
      relative_bias = performance$relative_bias,
      relative_rmse = performance$relative_rmse,
      coverage = performance$coverage,
      interval_length = performance$interval_length,
      time_seconds_per_dataset = mean(rows$sampling_seconds),
      simulations = nrow(rows),
      stringsAsFactors = FALSE
    )
    convergence[[length(convergence) + 1L]] <- data.frame(
      scenario_id = scenario$scenario_id,
      datasets = nrow(rows),
      rhat_threshold = settings$rhat,
      ess_threshold = settings$ess,
      stopping_rule_passes = sum(rows$converged),
      stopping_rule_pass_rate = mean(rows$converged),
      stringsAsFactors = FALSE
    )
    dataset_rows[[length(dataset_rows) + 1L]] <- rows
    if (all(vapply(records, function(x) !is.null(x$p_mean), logical(1)))) {
      cells <- expand.grid(
        week = seq_len(nrow(records[[1L]]$p_mean)),
        stock = seq_len(ncol(records[[1L]]$p_mean))
      )
      for (record in records) {
        item <- cells
        item$scenario_id <- scenario$scenario_id
        item$dataset_id <- record$dataset_id
        item$posterior_mean <- as.vector(record$p_mean)
        p_rows[[length(p_rows) + 1L]] <- item
      }
    }
  }

  simulation_path <- file.path(output_path, "simdata", "baseline")
  frequentist <- summarize_frequentist_simulations(
    simulation_path, true_total = settings$truth
  )
  frequentist_performance <- data.frame(
    scenario_id = frequentist$scenario_id,
    display_name = frequentist$display_name,
    method_family = "Frequentist",
    observation_model = NA_character_,
    prior = NA_character_,
    relative_bias = frequentist$relative_bias,
    relative_rmse = frequentist$relative_rmse,
    coverage = frequentist$coverage,
    interval_length = frequentist$interval_length,
    time_seconds_per_dataset = frequentist$seconds_per_dataset,
    simulations = frequentist$simulations,
    stringsAsFactors = FALSE
  )

  performance <- frequentist_performance
  if (length(bayesian_performance)) {
    performance <- rbind(do.call(rbind, bayesian_performance), performance)
  }
  results_path <- file.path(output_path, "results")
  dir.create(results_path, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(
    performance,
    file.path(results_path, "simulation_performance.csv"),
    row.names = FALSE
  )
  if (length(convergence)) {
    utils::write.csv(
      do.call(rbind, convergence),
      file.path(results_path, "simulation_convergence.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      do.call(rbind, dataset_rows),
      file.path(results_path, "bayesian_datasets.csv"),
      row.names = FALSE
    )
  }
  if (length(p_rows)) {
    p_means <- do.call(rbind, p_rows)
    truth <- readRDS(file.path(simulation_path, "manifest.rds"))$p_true
    p_means$true_p <- truth[cbind(p_means$week, p_means$stock)]
    p_means <- p_means[, c(
      "scenario_id", "dataset_id", "week", "stock", "posterior_mean",
      "true_p"
    )]
    utils::write.csv(
      p_means,
      file.path(results_path, "posterior_p_means.csv"),
      row.names = FALSE
    )
  }
  cat("\nResults written to: ", results_path, "\n", sep = "")
  invisible(performance)
}

taku_application_data <- function(input_path, enumerated_total) {
  if (!is.numeric(enumerated_total) || length(enumerated_total) != 1L ||
      !is.finite(enumerated_total) || enumerated_total <= 0) {
    stop("enumerated_total must be one positive finite number")
  }
  inputs <- load_taku_inputs(input_path)
  structure(list(
    Nw = enumerated_total,
    lam = inputs$lambda,
    mu = inputs$p,
    sd = inputs$sd,
    n = inputs$n,
    w = inputs$w,
    nweek = inputs$nweek,
    nstock = inputs$nstock
  ), class = "nlist")
}

# Run the four Bayesian application fits and the three closed-form estimators.
taku_run_application <- function(
    input_path,
    enumerated_total,
    output_path,
    settings = taku_paper_settings(n_datasets = 1L)) {
  output_path <- normalizePath(output_path, winslash = "/", mustWork = FALSE)
  if (dir.exists(output_path) &&
      length(list.files(output_path, all.files = TRUE, no.. = TRUE))) {
    stop("refusing to overwrite a non-empty application directory: ", output_path)
  }
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  data <- taku_application_data(input_path, enumerated_total)
  tasks <- taku_baseline_tasks(1L, 1L, settings$seed)
  rows <- vector("list", nrow(tasks))

  for (index in seq_len(nrow(tasks))) {
    task <- tasks[index, , drop = FALSE]
    cat("Fitting ", task$scenario_id, "...\n", sep = "")
    result <- taku_fit_baseline(data, task, settings)
    taku_atomic_rds(
      result,
      file.path(output_path, paste0(task$scenario_id, ".rds"))
    )
    rows[[index]] <- data.frame(
      scenario_id = task$scenario_id,
      estimate = result$row$estimate,
      standard_error = result$row$posterior_sd,
      lower = result$row$lower,
      upper = result$row$upper,
      time_seconds = result$row$sampling_seconds,
      stringsAsFactors = FALSE
    )
  }

  methods <- c(mom = "MoM", mom_alt = "MoM(Alt)", mom_naive = "MoM(Naive)")
  frequentist <- do.call(rbind, lapply(names(methods), function(method) {
    started <- proc.time()[["elapsed"]]
    fit <- estimate_frequentist_total(data, method)
    elapsed <- proc.time()[["elapsed"]] - started
    critical <- stats::qnorm(.975)
    data.frame(
      scenario_id = paste0("frequentist_", method),
      estimate = fit$estimate,
      standard_error = fit$standard_error,
      lower = fit$estimate - critical * fit$standard_error,
      upper = fit$estimate + critical * fit$standard_error,
      time_seconds = elapsed,
      stringsAsFactors = FALSE
    )
  }))
  combined <- rbind(
    do.call(rbind, rows),
    frequentist
  )
  utils::write.csv(
    combined,
    file.path(output_path, "application_results.csv"),
    row.names = FALSE
  )
  invisible(combined)
}
