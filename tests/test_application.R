# Application checks use synthetic summaries and do not run MCMC.
source("R/taku.R")

local({
  path <- tempfile("taku-application-", fileext = ".RData")
  on.exit(unlink(path), add = TRUE)
  p <- matrix(rep(c(.2, .3, .4, .1), 12L), nrow = 12L, byrow = TRUE)
  n <- rep(100L, 12L)
  w <- seq_len(12L)
  region1mean <- 100 * p
  region1sd <- matrix(2, 12L, 4L)
  save(n, w, region1mean, region1sd, file = path)

  set.seed(1234)
  before <- .Random.seed
  result <- estimate_taku_application(path, 35000)
  stopifnot(
    identical(before, .Random.seed),
    nrow(result) == 3L,
    length(unique(result$estimate)) == 1L,
    all(result$standard_error > 0),
    all(result$lower < result$estimate),
    all(result$upper > result$estimate)
  )

  data <- taku_application_data(path, 35000)
  stopifnot(
    inherits(data, "nlist"), data$nweek == 12L, data$nstock == 4L,
    data$Nw == 35000, isTRUE(taku_validate_baseline_data(data))
  )
})

cat("Application checks passed.\n")
