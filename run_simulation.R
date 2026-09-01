# Reproduce the baseline simulation study.
#
# Usage from the repository root:
#   Rscript run_simulation.R path/to/region1.RData path/to/output [scenario] [workers]
#
# `scenario` may be `all` (default), rdm_ar, rdm_dirichlet, mmd_ar, or
# mmd_dirichlet.  `workers` defaults to 1; increasing it fits independent
# datasets in parallel on the current computer.
arguments <- commandArgs(trailingOnly = TRUE)
if (!length(arguments) %in% 2:4) {
  stop(paste(
    "usage: Rscript run_simulation.R",
    "<input.RData> <output-directory> [scenario=all] [workers=1]"
  ))
}

source("R/taku.R")

input_path <- arguments[1L]
output_path <- arguments[2L]
scenario <- if (length(arguments) >= 3L) arguments[3L] else "all"
workers <- if (length(arguments) >= 4L) as.integer(arguments[4L]) else 1L

# These are the settings used for the paper results.  To conduct a small wiring
# check, pass a modified object explicitly rather than editing internal code,
# for example: taku_paper_settings(n_datasets = 2L, n_save = 100L).
settings <- taku_paper_settings()

taku_run_simulation(
  input_path = input_path,
  output_path = output_path,
  scenario = scenario,
  workers = workers,
  settings = settings
)
