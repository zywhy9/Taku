# One-dataset baseline fitting helpers. The locked simanalyse sampler owns the
# adaptation/sampling/stopping loop; these helpers add validation, timing and
# reporting only. Source model_code.R and scenarios.R before using this file.

# Atomic completion: readers see either no result or a fully written RDS.
# Existing results are never silently overwritten (including on Windows).
taku_atomic_rds <- function(object, path) {
  if (file.exists(path)) stop("Refusing to overwrite: ", path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(object, temporary, compress = FALSE)
  if (!file.rename(temporary, path)) stop("Cannot finalize: ", path)
  invisible(path)
}

# These are the ORIGINAL archived datasets, not regenerated baseline inputs.
# In particular, lam (fitted precision) is used unchanged, not lambda (truth).
taku_validate_baseline_data <- function(d) {
  stopifnot(inherits(d, "nlist"), d$nweek == 12, d$nstock == 4,
    identical(dim(d$mu), c(12L, 4L)), identical(dim(d$sd), c(12L, 4L)),
    length(d$n) == 12L, length(d$w) == 12L, length(d$lam) == 12L,
    all(is.finite(unlist(d))), all(d$n > 0), all(d$n == round(d$n)),
    all(d$lam > 0), all(d$sd >= 0), all(d$mu > 0), all(d$mu < 1),
    all(d$w >= 0), abs(sum(d$w) - 1) < 1e-10, d$Nw > 0)
  # Archived generation clips mu at the boundary. Do not renormalize it here:
  # that would silently change the data used by the original study.
  fitted <- pmax(1e-10, rowSums((d$mu * (1 - d$mu))^2) /
    rowSums(d$mu * (1 - d$mu) * d$sd^2) - 1)
  stopifnot(max(abs(fitted - d$lam)) < 1e-8)
  invisible(TRUE)
}

# Six balanced shards per scenario. Dataset identity and seeds never depend
# on node number, job ID, scheduling order or resubmission count.
# 1234 is the base seed, not a seed reused for every fit. The scenario's fixed
# row in taku_scenarios() supplies a 10,000-wide offset, then dataset_id is added.
# Record the resulting seed with each fit so that an isolated retry is identical.
taku_baseline_tasks <- function(n = 1000L, shards = 6L, seed_base = 1234L) {
  stopifnot(n >= shards, shards >= 1L)
  scenarios <- taku_scenarios()
  do.call(rbind, lapply(seq_len(nrow(scenarios)), function(i) {
    data.frame(scenario_id = scenarios$scenario_id[i],
      observation_model = scenarios$observation_model[i], prior = scenarios$prior[i],
      dataset_id = seq_len(n), shard = (seq_len(n) - 1L) %% shards + 1L,
      seed = as.integer(seed_base + i * 10000L + seq_len(n)),
      key = sprintf("%s_%07d", scenarios$scenario_id[i], seq_len(n)),
      stringsAsFactors = FALSE)
  }))
}

# Exactly the native simanalyse controls. The paper-facing stopping targets are
# N, plus phi for AR(1). ESS=100 is a modest Monte Carlo precision floor.
# Weekly p is output-only and deviance is not monitored; neither controls
# stopping. Only compact p means are retained. No transformed diagnostics are
# used.
taku_baseline_mode <- function(settings, prior) {
  target <- if (prior == "ar1") c("N", "phi") else "N"
  simanalyse::sma_set_mode("paper", n.chains = settings$n_chains,
    n.adapt = settings$n_adapt, n.save = settings$n_save,
    max.iter = settings$max_iter, max.time = settings$max_time_minutes,
    units = "mins", ess = settings$ess, r.hat = settings$rhat,
    ess.nodes = target, r.hat.nodes = target,
    normalize = FALSE)
}

# Observe rjags entry/exit without replacing its sampler or modifying the
# simanalyse loop. Each worker is an independent R process running ONE fit.
# Instrumentation is removed even if sampling throws an error. Compilation +
# adaptation, sampling calls, and complete fit elapsed time are separate.
taku_jags_clock <- function() {
  clock <- new.env(parent = emptyenv())
  clock$events <- list()
  clock$starts <- list()
  clock$start <- function(stage, iterations = NA_real_) {
    clock$starts[[stage]] <- list(time = proc.time(), iterations = iterations)
    cat(sprintf("[%s] %s started; iterations=%s\n", format(Sys.time()),
      stage, format(iterations, scientific = FALSE)))
    flush.console()
  }
  clock$end <- function(stage) {
    start <- clock$starts[[stage]]
    delta <- proc.time() - start$time
    clock$events[[length(clock$events) + 1L]] <- data.frame(stage = stage,
      iterations = as.numeric(start$iterations), elapsed_seconds = delta[["elapsed"]],
      cpu_seconds = unname(delta[["user.self"]] + delta[["sys.self"]]))
    cat(sprintf("[%s] %s finished; elapsed=%.2fs\n", format(Sys.time()),
      stage, delta[["elapsed"]]))
    flush.console()
  }
  clock
}

taku_call_simanalyse <- function(arguments, analyser = NULL, namespace = NULL) {
  # Optional dependency injection is for tests; production uses both defaults.
  if (is.null(analyser)) analyser <- getFromNamespace("analyse_dataset_bayesian", "simanalyse")
  if (!setequal(names(arguments), names(formals(analyser))))
    stop("The locked simanalyse interface does not match this wrapper")
  clock <- taku_jags_clock()
  if (is.null(namespace)) namespace <- asNamespace("rjags")
  installed <- character()
  on.exit(for (name in installed) untrace(name, where = namespace), add = TRUE)
  for (name in c("jags.model", "jags.samples")) {
    if (inherits(get(name, envir = namespace), "functionWithTrace"))
      stop("rjags already has instrumentation: ", name)
    stage <- if (name == "jags.model") "compile_adapt" else "sampling"
    entry <- if (name == "jags.model")
      substitute(CLOCK$start(STAGE, n.adapt), list(CLOCK = clock, STAGE = stage)) else
      substitute(CLOCK$start(STAGE, n.iter), list(CLOCK = clock, STAGE = stage))
    leave <- substitute(CLOCK$end(STAGE), list(CLOCK = clock, STAGE = stage))
    trace(name, where = namespace, tracer = entry, exit = leave, print = FALSE)
    installed <- c(installed, name)
  }
  start <- proc.time()
  result <- do.call(analyser, arguments)
  list(draws = result, events = do.call(rbind, clock$events),
    fit_elapsed_seconds = unname((proc.time() - start)[["elapsed"]]))
}

# Finite diagnostics for every stopping target are required. simanalyse uses
# na.rm=TRUE internally; this audit must not silently accept a missing or
# undefined target R-hat. Such fits are retained but not labelled converged.
taku_baseline_diagnostics <- function(result, prior, settings) {
  target <- if (prior == "ar1") c("N", "phi") else "N"
  stopifnot(all(target %in% names(result)))
  selected <- subset(result, pars = target)
  rh <- mcmcr::rhat(selected, as_df = TRUE, by = "term")
  es <- mcmcr::ess(selected, as_df = TRUE, by = "term")
  stopifnot(identical(as.character(rh$term), as.character(es$term)))
  terms <- data.frame(term = as.character(rh$term), rhat = rh$rhat, ess = es$ess)
  expected <- sum(vapply(selected, function(x) prod(dim(x)[-(1:2)]), 0.0))
  stopifnot(nrow(terms) == expected, !anyDuplicated(terms$term),
    all(target %in% terms$term))
  terms$ess_required <- TRUE
  terms$rhat_pass <- is.finite(terms$rhat) & terms$rhat <= settings$rhat
  terms$ess_pass <- is.finite(terms$ess) & terms$ess >= settings$ess
  terms$rhat_excess <- ifelse(is.finite(terms$rhat), pmax(0, terms$rhat - settings$rhat), NA_real_)
  terms$ess_shortfall <- ifelse(is.finite(terms$ess), pmax(0, settings$ess - terms$ess), NA_real_)
  list(terms = terms, converged = all(terms$rhat_pass) &&
    all(terms$ess_pass[terms$ess_required]))
}

# One independent fit; returns the last retained simanalyse batch. Earlier
# batches are not concatenated, and n.save is not the total raw iteration count.
# 'sampler' is injectable solely for no-MCMC unit tests.
taku_fit_baseline <- function(data, task, settings, sampler = taku_call_simanalyse) {
  taku_validate_baseline_data(data)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(task$seed)
  mode <- taku_baseline_mode(settings, task$prior)
  inits <- scenario_initial_values(task$prior, data$nweek, data$nstock, mode$n.chains)
  code <- getFromNamespace("prepare_code", "simanalyse")(
    analysis_model_code(task$observation_model), prior_model_code(task$prior), NULL)
  # sd and lambda are retained in each generated dataset for the frequentist
  # estimators, but they are not JAGS inputs. Omitting them avoids distracting
  # unused-variable warnings without changing the fitted model.
  jags_data <- data[c("mu", "lam", "n", "nweek", "nstock", "Nw", "w")]
  output_nodes <- if (isTRUE(settings$save_p_means)) "p" else character()
  arguments <- c(list(nlistdata = jags_data, code = code,
    monitor = unique(c(mode$r.hat.nodes, mode$ess.nodes, output_nodes)), deviance = FALSE,
    inits = inits, quiet = FALSE), mode)
  sampled <- sampler(arguments)
  result <- sampled$draws
  stopifnot(inherits(result, "mcmcr"), all(is.finite(unlist(result))))
  expected_dim <- c(as.integer(mode$n.chains), as.integer(mode$n.save))
  stopifnot(identical(as.integer(dim(result$N)[1:2]), expected_dim))
  N <- as.numeric(result$N)
  p_mean <- NULL
  if (isTRUE(settings$save_p_means)) {
    stopifnot(
      "p" %in% names(result),
      identical(as.integer(dim(result$p)[1:2]), expected_dim),
      identical(
        as.integer(dim(result$p)[3:4]),
        as.integer(c(data$nweek, data$nstock))
      )
    )
    # Monitoring p changes neither the sampler nor the stopping targets. This
    # identity confirms that its compact posterior mean belongs to the same
    # draws as N.
    denominator <- numeric(length(N))
    for (week in seq_len(data$nweek)) {
      denominator <- denominator + data$w[week] *
        as.numeric(result$p[, , week, 2L] + result$p[, , week, 3L])
    }
    stopifnot(max(abs(N - data$Nw / denominator)) < 1e-5)
    p_mean <- apply(result$p, c(3L, 4L), mean)
    stopifnot(
      identical(
        as.integer(dim(p_mean)),
        as.integer(c(data$nweek, data$nstock))
      ),
      all(is.finite(p_mean)),
      max(abs(rowSums(p_mean) - 1)) < 1e-10
    )
  }
  diagnostics <- taku_baseline_diagnostics(result, task$prior, settings)
  events <- sampled$events
  sampling <- events[events$stage == "sampling", , drop = FALSE]
  batches <- nrow(sampling)
  stopifnot(batches >= 1L,
    identical(as.numeric(sampling$iterations), as.numeric(mode$n.save * seq_len(batches))))
  raw <- sum(sampling$iterations)
  iteration_limit <- raw + mode$n.save * (batches + 1) > mode$max.iter
  # Time is predictive inside simanalyse, not a hard measured timeout. Avoid
  # asserting an exact cause when the package exposes no structured reason.
  reason <- if (diagnostics$converged) "converged" else if
    (iteration_limit) "iteration_budget_rule_unmet" else "time_budget_or_rule_unmet"
  q <- unname(quantile(N, c(.025, .975)))
  t <- diagnostics$terms
  row <- data.frame(dataset_id = task$dataset_id, scenario_id = task$scenario_id,
    shard = task$shard, seed = task$seed, estimate = mean(N), median = median(N),
    posterior_sd = sd(N), lower = q[1], upper = q[2],
    converged = diagnostics$converged, stop_reason = reason,
    max_rhat = if (all(is.finite(t$rhat))) max(t$rhat) else NA_real_,
    min_target_ess = if (all(is.finite(t$ess[t$ess_required]))) min(t$ess[t$ess_required]) else NA_real_,
    N_rhat = t$rhat[match("N", t$term)], N_ess = t$ess[match("N", t$term)],
    rhat_failed_terms = sum(!t$rhat_pass), ess_failed_targets = sum(!t$ess_pass[t$ess_required]),
    batches = batches, raw_iterations_per_chain = raw, final_thin = batches,
    sampling_seconds = sum(sampling$elapsed_seconds), sampling_cpu_seconds = sum(sampling$cpu_seconds),
    compile_adapt_seconds = sum(events$elapsed_seconds[events$stage == "compile_adapt"]),
    fit_elapsed_seconds = sampled$fit_elapsed_seconds)
  # Scalar target draws suffice to reconstruct all abundance performance
  # measures and keep the production archive compact.
  kept <- if (isTRUE(settings$save_full_draws)) result else subset(result, pars = mode$ess.nodes)
  list(row = row, diagnostics = t, p_mean = p_mean, draws = kept,
    timing = events, mode = mode,
    initial_values = inits, precision_used = data$lam, rng_kind = RNGkind())
}

# Summary denominators always expose the number expected and available.
# Budget-limited fits remain in performance. Missing fits prevent finalization.
taku_baseline_performance <- function(rows, truth = 60000, expected = 1000L) {
  stopifnot(nrow(rows) > 0L, all(is.finite(rows$estimate)),
    all(is.finite(rows$lower)), all(is.finite(rows$upper)))
  covered <- rows$lower <= truth & rows$upper >= truth
  data.frame(n_expected = expected, n_available = nrow(rows),
    relative_bias = mean(rows$estimate / truth - 1),
    relative_rmse = sqrt(mean((rows$estimate / truth - 1)^2)),
    coverage = mean(covered), coverage_mcse = sqrt(mean(covered) * (1 - mean(covered)) / nrow(rows)),
    interval_length = mean(rows$upper - rows$lower))
}
