# Shared R entry point for the generic Slurm wrappers.
#
# Readers normally use cluster/submit.sh. This file keeps preparation, one-fit
# execution, and final summarization in one place while delegating all
# scientific calculations to the same R/taku.R functions used locally.
arguments <- commandArgs(trailingOnly = TRUE)
if (!length(arguments)) stop("missing cluster command")

command <- arguments[1L]

load_repository <- function(path) {
  repository <- normalizePath(path, winslash = "/", mustWork = TRUE)
  setwd(repository)
  source(file.path("R", "taku.R"))
  repository
}

check_packages <- function() {
  required <- c("coda", "mcmcr", "nlist", "rjags", "simanalyse")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("missing R packages: ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

if (identical(command, "check")) {
  if (length(arguments) != 4L) {
    stop("check requires repository, input file, and output directory")
  }
  load_repository(arguments[2L])
  check_packages()
  load_taku_inputs(arguments[3L])
  output_parent <- dirname(normalizePath(
    arguments[4L], winslash = "/", mustWork = FALSE
  ))
  if (!dir.exists(output_parent)) {
    stop("output parent does not exist: ", output_parent)
  }
  cat("PASS: packages and empirical input verified; no MCMC run.\n")
} else if (identical(command, "prepare")) {
  if (length(arguments) != 4L) {
    stop("prepare requires repository, input file, and output directory")
  }
  load_repository(arguments[2L])
  check_packages()
  prepared <- taku_prepare_simulation(
    input_path = arguments[3L],
    output_path = arguments[4L],
    settings = taku_paper_settings()
  )
  cat(
    "Prepared ", prepared$manifest$nsims,
    " baseline datasets in ", prepared$simulation_path, "\n",
    sep = ""
  )
} else if (identical(command, "fit")) {
  if (length(arguments) != 5L) {
    stop("fit requires repository, output directory, scenario, and dataset ID")
  }
  load_repository(arguments[2L])
  check_packages()
  output_path <- normalizePath(arguments[3L], winslash = "/", mustWork = TRUE)
  scenario_id <- taku_scenario_ids(arguments[4L])
  if (length(scenario_id) != 1L) stop("fit requires one scenario")
  dataset_id <- as.integer(arguments[5L])
  settings <- readRDS(file.path(output_path, "settings.rds"))
  if (!taku_settings_equal(settings, taku_paper_settings())) {
    stop("cluster output does not use the paper settings")
  }
  if (is.na(dataset_id) || dataset_id < 1L ||
      dataset_id > settings$n_datasets) {
    stop("dataset ID is outside the configured range")
  }

  tasks <- taku_baseline_tasks(
    n = settings$n_datasets,
    shards = 1L,
    seed_base = settings$seed
  )
  task <- tasks[
    tasks$scenario_id == scenario_id & tasks$dataset_id == dataset_id,
    , drop = FALSE
  ]
  stopifnot(nrow(task) == 1L)
  data_file <- file.path(
    output_path, "simdata", "baseline", sprintf("data%07d.rds", dataset_id)
  )
  if (!file.exists(data_file)) stop("missing simulated dataset: ", data_file)

  record <- taku_fit_simulation_dataset(
    data_file = data_file,
    task = task,
    output_path = output_path,
    settings = settings
  )
  cat(
    "Completed ", scenario_id, " dataset ", dataset_id,
    "; iterations per chain = ", record$row$raw_iterations_per_chain,
    "\n",
    sep = ""
  )
} else if (identical(command, "summarize")) {
  if (length(arguments) != 3L) {
    stop("summarize requires repository and output directory")
  }
  load_repository(arguments[2L])
  output_path <- normalizePath(arguments[3L], winslash = "/", mustWork = TRUE)
  taku_summarize_simulation(output_path)
  cat("Simulation summaries completed in ", output_path, "\n", sep = "")
} else {
  stop("unknown cluster command: ", command)
}
