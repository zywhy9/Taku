# Estimate a separate Dirichlet precision parameter for every week.
#
# For a Dirichlet composition with mean p and total concentration lambda,
# Var(P_j) = p_j(1 - p_j) / (lambda + 1). For each row, this function regresses
# the reported marginal variances on p_j(1 - p_j), with no intercept, and
# converts the fitted slope to lambda = 1 / slope - 1.
estimate_dirichlet_precision <- function(proportions, standard_deviations) {
  proportions <- as.matrix(proportions)
  standard_deviations <- as.matrix(standard_deviations)

  if (!identical(dim(proportions), dim(standard_deviations))) {
    stop("proportions and standard_deviations must have identical dimensions")
  }
  if (anyNA(proportions) || anyNA(standard_deviations) ||
      any(!is.finite(proportions)) || any(!is.finite(standard_deviations))) {
    stop("proportions and standard_deviations must be finite and non-missing")
  }
  if (any(proportions < 0 | proportions > 1) ||
      any(standard_deviations < 0)) {
    stop("invalid proportion or standard-deviation value")
  }

  x <- proportions * (1 - proportions)
  y <- standard_deviations^2
  denominator <- rowSums(x^2)

  if (any(denominator <= 0)) {
    stop("Dirichlet precision is undefined for a degenerate row")
  }

  slope <- rowSums(x * y) / denominator
  precision <- 1 / slope - 1

  if (any(!is.finite(precision)) || any(precision <= 0)) {
    stop("estimated Dirichlet precision must be positive and finite")
  }

  unname(precision)
}

# Load and validate the empirical inputs without attaching them to the global
# R environment.
load_taku_inputs <- function(path, percentages = TRUE) {
  if (length(path) != 1L || !file.exists(path)) {
    stop("input file does not exist: ", path)
  }

  input_environment <- new.env(parent = emptyenv())
  loaded_names <- load(path, envir = input_environment)
  required_names <- c("n", "w", "region1mean", "region1sd")
  missing_names <- setdiff(required_names, loaded_names)

  if (length(missing_names)) {
    stop("input file is missing: ", paste(missing_names, collapse = ", "))
  }

  n <- get("n", envir = input_environment, inherits = FALSE)
  w <- get("w", envir = input_environment, inherits = FALSE)
  proportions <- as.matrix(
    get("region1mean", envir = input_environment, inherits = FALSE)
  )
  standard_deviations <- as.matrix(
    get("region1sd", envir = input_environment, inherits = FALSE)
  )

  if (percentages) {
    proportions <- proportions / 100
    standard_deviations <- standard_deviations / 100
  }

  if (length(n) != nrow(proportions) || length(w) != nrow(proportions)) {
    stop("n and w must have one value per week")
  }
  if (!identical(dim(proportions), dim(standard_deviations))) {
    stop("region1mean and region1sd must have identical dimensions")
  }
  if (anyNA(n) || anyNA(w) || any(!is.finite(n)) || any(!is.finite(w)) ||
      any(n <= 0) || any(n != floor(n)) || any(w <= 0)) {
    stop("n must contain positive integers and w must contain positive values")
  }
  if (any(abs(rowSums(proportions) - 1) > 1e-6)) {
    stop("each row of region1mean must sum to one after scaling")
  }

  weights <- as.numeric(w) / sum(w)

  list(
    n = as.integer(n),
    w = weights,
    p = unname(proportions),
    sd = unname(standard_deviations),
    lambda = estimate_dirichlet_precision(proportions, standard_deviations),
    nweek = nrow(proportions),
    nstock = ncol(proportions)
  )
}
