# Reproduce the Taku data application.
#
# Usage from the repository root:
#   Rscript run_application.R path/to/region1.RData <lake-count> path/to/output
#
# The private input and generated output remain outside the repository.  The
# script runs all four Bayesian scenarios and the three frequentist estimators.
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 3L) {
  stop(paste(
    "usage: Rscript run_application.R",
    "<input.RData> <enumerated-lake-total> <output-directory>"
  ))
}

source("R/taku.R")

taku_run_application(
  input_path = arguments[1L],
  enumerated_total = as.numeric(arguments[2L]),
  output_path = arguments[3L],
  settings = taku_paper_settings(n_datasets = 1L, save_p_means = FALSE)
)
