# Generate the baseline simulation datasets used in the paper.
#
# The data-generating hierarchy is
#
#   p_true -> X -> mu_observed,
#
# where X contains the latent multinomial sample counts and mu_observed is the
# reported IA composition generated from the RDM working model.

# Draw one Dirichlet vector without requiring an MCMC package.
draw_taku_dirichlet <- function(alpha) {
  if (!is.numeric(alpha) || anyNA(alpha) || any(!is.finite(alpha)) ||
      any(alpha <= 0)) {
    stop("Dirichlet concentrations must be positive and finite")
  }

  draws <- stats::rgamma(length(alpha), shape = alpha, rate = 1)
  total <- sum(draws)
  if (!is.finite(total) || total <= 0) {
    stop("Dirichlet draw underflowed for all categories")
  }

  draws / total
}

# Convert one week-by-category count matrix into reported IA summaries.
simulate_taku_ia_summary <- function(
    observed_counts,
    sample_sizes,
    input_precision) {
  observed_counts <- as.matrix(observed_counts)
  sample_sizes <- as.integer(sample_sizes)
  input_precision <- as.numeric(input_precision)

  if (anyNA(observed_counts) || any(!is.finite(observed_counts)) ||
      any(observed_counts < 0) || any(observed_counts != floor(observed_counts))) {
    stop("observed_counts must contain non-negative integers")
  }
  if (length(sample_sizes) != nrow(observed_counts) ||
      length(input_precision) != nrow(observed_counts) ||
      anyNA(sample_sizes) || any(sample_sizes <= 0L) ||
      anyNA(input_precision) || any(!is.finite(input_precision)) ||
      any(input_precision <= 0)) {
    stop("sample sizes and input precisions must match the count rows")
  }
  if (any(rowSums(observed_counts) != sample_sizes)) {
    stop("each observed-count row must sum to its sample size")
  }

  nweek <- nrow(observed_counts)
  nstock <- ncol(observed_counts)
  mu <- matrix(NA_real_, nrow = nweek, ncol = nstock)
  standard_deviation <- matrix(NA_real_, nrow = nweek, ncol = nstock)
  fitted_precision <- numeric(nweek)

  for (week in seq_len(nweek)) {
    center <- pmax(1e-10, observed_counts[week, ] / sample_sizes[week])
    variance <- center * (1 - center) / (input_precision[week] + 1)
    standard_deviation[week, ] <- sqrt(variance)

    raw_mu <- draw_taku_dirichlet(input_precision[week] * center)
    mu[week, ] <- pmin(pmax(1e-10, raw_mu), 1 - 1e-7)

    a <- mu[week, ] * (1 - mu[week, ])
    denominator <- sum(a^2)
    numerator <- sum(a * variance)
    fitted_precision[week] <- max(1e-10, denominator / numerator - 1)
  }

  if (any(!is.finite(fitted_precision)) || any(fitted_precision <= 0)) {
    stop("simulated moment-matched precisions must be positive and finite")
  }

  list(mu = mu, sd = standard_deviation, lam = fitted_precision)
}

# Simulate one baseline replicate.
simulate_taku_dataset <- function(inputs, true_total = 60000) {
  required_inputs <- c("n", "w", "p", "lambda", "nweek", "nstock")
  if (!is.list(inputs) || !all(required_inputs %in% names(inputs))) {
    stop("inputs must be returned by load_taku_inputs()")
  }
  if (length(true_total) != 1L || !is.finite(true_total) || true_total <= 0) {
    stop("true_total must be a positive finite value")
  }
  if (inputs$nstock != nrow(taku_reporting_units())) {
    stop("the Taku DGP requires exactly four reporting categories")
  }

  sample_counts <- t(vapply(
    seq_len(inputs$nweek),
    function(week) {
      as.integer(stats::rmultinom(
        1L,
        size = inputs$n[week],
        prob = inputs$p[week, ]
      )[, 1L])
    },
    integer(inputs$nstock)
  ))
  summary <- simulate_taku_ia_summary(
    sample_counts,
    inputs$n,
    inputs$lambda
  )

  lake_proportion <- rowSums(
    inputs$p[, taku_lake_categories(), drop = FALSE]
  )
  lake_abundance <- true_total * sum(inputs$w * lake_proportion)

  structure(
    list(
      Nw = lake_abundance,
      lam = summary$lam,
      mu = summary$mu,
      sd = summary$sd,
      n = inputs$n,
      lambda = inputs$lambda,
      w = inputs$w,
      nweek = inputs$nweek,
      nstock = inputs$nstock
    ),
    class = "nlist"
  )
}

# Write the metadata marker required by file-based simanalyse tools.
write_taku_sims_metadata <- function(path, nsims, seed_state, true_total) {
  metadata_path <- file.path(path, ".sims.rds")
  if (file.exists(metadata_path)) {
    stop("refusing to overwrite simanalyse metadata: ", metadata_path)
  }

  metadata <- list(
    code = "# Custom baseline Taku generator",
    constants = structure(list(), class = "nlist"),
    parameters = structure(list(N = true_total), class = "nlist"),
    monitor = ".*",
    nsims = as.integer(nsims),
    seed = as.integer(seed_state)
  )
  saveRDS(metadata, metadata_path)
  invisible(metadata_path)
}

# Validate the subset of sims metadata used by the file-based analysis API.
validate_taku_sims_metadata <- function(path, nsims, true_total) {
  metadata_path <- file.path(path, ".sims.rds")
  if (!file.exists(metadata_path)) {
    stop("simulation directory is missing .sims.rds: ", path)
  }
  metadata <- readRDS(metadata_path)
  required <- c("code", "constants", "parameters", "monitor", "nsims", "seed")
  valid <- is.list(metadata) && all(required %in% names(metadata)) &&
    inherits(metadata$constants, "nlist") &&
    inherits(metadata$parameters, "nlist") &&
    identical(metadata$nsims, as.integer(nsims)) &&
    isTRUE(all.equal(metadata$parameters$N, true_total, tolerance = 0)) &&
    is.integer(metadata$seed) && length(metadata$seed) > 0L
  if (!valid) {
    stop("invalid simanalyse metadata: ", metadata_path)
  }
  invisible(metadata)
}

# Validate an existing baseline run before fitting or reusing it.
validate_taku_simulation_run <- function(
    inputs,
    output_path,
    true_total = 60000) {
  output_path <- normalizePath(output_path, winslash = "/", mustWork = FALSE)
  simulation_path <- file.path(output_path, "simdata", "baseline")
  manifest_path <- file.path(simulation_path, "manifest.rds")

  if (!file.exists(manifest_path)) {
    stop("simulation manifest does not exist: ", manifest_path)
  }
  manifest <- readRDS(manifest_path)
  required_manifest <- c(
    "generator", "nsims", "seed", "true_total", "p_true", "variant_id"
  )
  if (!is.list(manifest) || !all(required_manifest %in% names(manifest))) {
    stop("simulation manifest is incomplete: ", manifest_path)
  }

  matches <- c(
    identical(manifest$generator, "baseline_rdm"),
    isTRUE(all.equal(manifest$true_total, true_total, tolerance = 0)),
    isTRUE(all.equal(manifest$p_true, inputs$p, tolerance = 1e-15)),
    identical(manifest$variant_id, "baseline")
  )
  if (!all(matches)) {
    stop("simulation manifest does not match the requested inputs: ", manifest_path)
  }

  expected_names <- sprintf("data%07d.rds", seq_len(manifest$nsims))
  actual_names <- list.files(
    simulation_path,
    pattern = "^data[0-9]{7}[.]rds$"
  )
  if (!identical(sort(actual_names), expected_names)) {
    stop("simulation files are incomplete or do not match the manifest")
  }
  validate_taku_sims_metadata(
    simulation_path,
    nsims = manifest$nsims,
    true_total = manifest$true_total
  )

  invisible(list(
    output_path = output_path,
    simulation_path = simulation_path,
    variant_id = "baseline",
    manifest = manifest
  ))
}

# Create baseline simulation files without fitting a Bayesian model.
generate_taku_datasets <- function(
    input_path,
    output_path,
    nsims = 1000L,
    seed = 1234L,
    true_total = 60000) {
  inputs <- load_taku_inputs(input_path)
  nsims <- as.integer(nsims)
  seed <- as.integer(seed)

  if (length(nsims) != 1L || is.na(nsims) || nsims <= 0L) {
    stop("nsims must be a positive integer")
  }
  if (length(seed) != 1L || is.na(seed)) {
    stop("seed must be a single integer")
  }

  output_path <- normalizePath(output_path, winslash = "/", mustWork = FALSE)
  simulation_path <- file.path(output_path, "simdata", "baseline")
  if (dir.exists(simulation_path) &&
      length(list.files(simulation_path, all.files = TRUE, no.. = TRUE))) {
    stop("refusing to overwrite or mix simulation files in: ", simulation_path)
  }
  dir.create(simulation_path, recursive = TRUE, showWarnings = FALSE)

  set.seed(seed)
  seed_state <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  for (simulation in seq_len(nsims)) {
    dataset <- simulate_taku_dataset(inputs, true_total = true_total)
    suffix <- sprintf("%07d.rds", simulation)
    saveRDS(dataset, file.path(simulation_path, paste0("data", suffix)))
  }

  write_taku_sims_metadata(
    path = simulation_path,
    nsims = nsims,
    seed_state = seed_state,
    true_total = true_total
  )
  manifest <- list(
    generator = "baseline_rdm",
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3%z"),
    nsims = nsims,
    seed = seed,
    true_total = true_total,
    p_true = inputs$p,
    lake_categories = taku_lake_categories(),
    variant_id = "baseline"
  )
  saveRDS(manifest, file.path(simulation_path, "manifest.rds"))

  invisible(list(
    output_path = output_path,
    simulation_path = simulation_path,
    variant_id = "baseline",
    manifest = manifest
  ))
}
