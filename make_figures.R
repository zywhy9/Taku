# Recreate Figures 1 and 6--8 from posterior means of the weekly
# stock proportions.
#
# This is the compact-data equivalent of the original plotting script in
# `simanalyse/pplot.R`. The original script first averaged
# `p[chain, draw, week, stock]` within each simulated dataset. The refactored
# simulation saves exactly those posterior means, so plotting does not require
# retaining or reloading the full MCMC draws.
#
# Usage from the repository root:
#   Rscript make_figures.R <simulation-or-review-directory> <figure-directory>

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) {
  stop(paste(
    "usage: Rscript make_figures.R",
    "<simulation-or-review-directory> <figure-directory>"
  ))
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("package ggplot2 is required; run renv::restore()")
}

input_directory <- normalizePath(
  arguments[1L], winslash = "/", mustWork = TRUE
)
figure_directory <- normalizePath(
  arguments[2L], winslash = "/", mustWork = FALSE
)

# A standard repository run writes the CSV below `results/`. The compact
# review archive produced on the cluster places the same CSV at its root.
input_candidates <- file.path(
  input_directory,
  c("results/posterior_p_means.csv", "posterior_p_means.csv")
)
input_file <- input_candidates[file.exists(input_candidates)]
if (length(input_file) != 1L) {
  stop(
    "expected exactly one posterior_p_means.csv in the input directory"
  )
}

posterior_means <- utils::read.csv(
  input_file, stringsAsFactors = FALSE
)
required_columns <- c(
  "scenario_id", "dataset_id", "week", "stock",
  "posterior_mean", "true_p"
)
if (!identical(names(posterior_means), required_columns)) {
  stop("posterior_p_means.csv has unexpected columns")
}

scenarios <- c(
  "rdm_ar", "rdm_dirichlet", "mmd_ar", "mmd_dirichlet"
)
expected_keys <- expand.grid(
  scenario_id = scenarios,
  dataset_id = seq_len(1000L),
  week = seq_len(12L),
  stock = seq_len(4L),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
key_columns <- c("scenario_id", "dataset_id", "week", "stock")

# These checks prevent incomplete or duplicated simulation results from being
# turned into manuscript figures.
observed_keys <- posterior_means[key_columns]
key_text <- function(data) {
  do.call(paste, c(data, sep = ":"))
}
stopifnot(
  nrow(posterior_means) == nrow(expected_keys),
  !anyDuplicated(observed_keys),
  identical(sort(key_text(observed_keys)), sort(key_text(expected_keys))),
  identical(sort(unique(posterior_means$scenario_id)), sort(scenarios)),
  all(is.finite(posterior_means$posterior_mean)),
  all(posterior_means$posterior_mean >= 0),
  all(posterior_means$posterior_mean <= 1)
)

posterior_means$week <- factor(
  posterior_means$week, levels = seq_len(12L)
)
posterior_means$stock <- factor(
  posterior_means$stock, levels = seq_len(4L)
)

# Each panel has one true proportion shared by all 1,000 datasets. This table
# supplies the black dashed horizontal reference lines described in the
# manuscript caption. The line does not rely on colour and remains clear in
# grayscale printouts.
truth <- unique(
  posterior_means[c("scenario_id", "week", "stock", "true_p")]
)
stopifnot(nrow(truth) == length(scenarios) * 12L * 4L)

# Figure numbering follows the manuscript: MMD--AR(1) is Figure 1 and the
# remaining specifications appear as Figures 6--8 in the supplementary file.
figure_names <- c(
  mmd_ar = "Figure_1",
  rdm_dirichlet = "Figure_6",
  rdm_ar = "Figure_7",
  mmd_dirichlet = "Figure_8"
)

dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)
destinations <- file.path(
  figure_directory, paste0(unname(figure_names), ".pdf")
)
if (any(file.exists(destinations))) {
  stop("refusing to overwrite existing PDF files in: ", figure_directory)
}

for (scenario in names(figure_names)) {
  values <- posterior_means[
    posterior_means$scenario_id == scenario, , drop = FALSE
  ]
  reference <- truth[
    truth$scenario_id == scenario, , drop = FALSE
  ]

  # Preserve the original plot: reporting units in rows, weeks in columns,
  # ggplot2's default boxplot theme, and a black dashed line for the true value.
  figure <- ggplot2::ggplot(
    values, ggplot2::aes(y = posterior_mean)
  ) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_hline(
      data = reference,
      ggplot2::aes(yintercept = true_p),
      colour = "black",
      linetype = "dashed"
    ) +
    ggplot2::facet_grid(stock ~ week) +
    # Plotmath treats pi[k, t] as array-style indexing and silently displays
    # only k. Explicit concatenation is required for the intended pi_{k,t}.
    ggplot2::labs(y = expression(pi[k * "," * t])) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )

  ggplot2::ggsave(
    file.path(
      figure_directory, paste0(figure_names[[scenario]], ".pdf")
    ),
    figure,
    width = 8,
    height = 8,
    units = "in",
    device = grDevices::cairo_pdf
  )
}

cat(
  "Created Figures 1 and 6--8 in: ",
  figure_directory,
  "\n",
  sep = ""
)
