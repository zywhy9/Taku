# Build the JAGS likelihood for the observed weekly compositions.
#
# Arguments:
#   observation_model:
#     "rdm" uses a random-Dirichlet-multinomial representation. Sequential
#     binomial draws reconstruct multinomial counts X, after which mu is modeled
#     around X / n using weekly precision lam.
#     "mmd" uses a moment-matched Dirichlet approximation centered directly at
#     the latent composition p. Its effective precision lambdat incorporates
#     both multinomial sampling variation and the reported Dirichlet variation.
#
# Returns:
#   A character scalar containing the selected JAGS likelihood block. It must be
#   combined with one block from prior_model_code().
analysis_model_code <- function(observation_model = c("rdm", "mmd")) {
  observation_model <- match.arg(observation_model)

  if (observation_model == "rdm") {
    return(paste(
      "for(week in 1:nweek){",
      "  X[week,1] ~ dbin(p[week,1],n[week])",
      "  for(s in 2:(nstock-1)){",
      "    X[week,s] ~ dbin(p[week,s]/(1-sum(p[week,1:(s-1)])),n[week]-sum(X[week,1:(s-1)]))",
      "  }",
      "  X[week,nstock] <- n[week] - sum(X[week,1:(nstock-1)])",
      "  for(stock in 1:nstock){",
      "    mn[week,stock] <- max(1e-10,X[week,stock]/n[week])",
      "  }",
      "  mu[week,1:nstock] ~ ddirch(mn[week,1:nstock]*lam[week])",
      "}",
      sep = "\n"
    ))
  }

  paste(
    "for(week in 1:nweek){",
    "  lambdat[week] <- 1/(1/n[week]+(1-1/n[week])/(lam[week]+1))-1",
    "  mu[week,1:nstock] ~ ddirch(p[week,1:nstock]*lambdat[week])",
    "}",
    sep = "\n"
  )
}

# Build the JAGS prior and abundance-transformation block.
#
# Arguments:
#   prior:
#     "ar1" models an unconstrained latent process Z separately for each stock,
#     links adjacent weeks through the shared autocorrelation phi, and converts
#     Z to a composition using a softmax transformation.
#     "dirichlet" assigns each week's composition an independent uniform
#     Dirichlet(1, ..., 1) prior.
#
# In both cases, total abundance N is derived by dividing the observed lake-stock
# abundance Nw by the weighted latent proportion belonging to stocks 2 and 3.
#
# Returns:
#   A character scalar containing the selected JAGS prior block plus the formula
#   for N.
prior_model_code <- function(prior = c("ar1", "dirichlet")) {
  prior <- match.arg(prior)

  if (prior == "ar1") {
    return(paste(
      "for(stock in 1:nstock){",
      "  Z[1,stock] ~ dnorm(0,1/2^2)",
      "  expZ[1,stock] <- exp(Z[1,stock])",
      "  for(week in 2:nweek){",
      "    Z[week,stock] <- phi*Z[(week-1),stock]+epsilon[week,stock]",
      "    expZ[week,stock] <- exp(Z[week,stock])",
      "    epsilon[week,stock] ~ dnorm(0,1/(2^2*(1-phi^2)))",
      "  }",
      "  for(week in 1:nweek){",
      "    p[week,stock] <- expZ[week,stock]/sum(expZ[week,])",
      "  }",
      "}",
      "phi ~ dunif(-1,1)",
      "N <- Nw/sum(w*p[,2],w*p[,3])",
      sep = "\n"
    ))
  }

  paste(
    "for(week in 1:nweek){",
    "  p[week,1:nstock] ~ ddirch(rep(1,nstock))",
    "}",
    "N <- Nw/sum(w*p[,2],w*p[,3])",
    sep = "\n"
  )
}

# Construct dispersed initial-value lists for the selected prior.
#
# Arguments:
#   prior: "ar1" or "dirichlet".
#   nweek: Number of weekly compositions in the model.
#   nstock: Number of stock categories per week.
#   n_chains: Number of chain-specific initial-value lists to return.
#
# Details:
#   Four AR(1) chains start phi at 0.3, 0.5, 0.7, and 0.9. Three-chain runs
#   retain the historical 0.5, 0.7, and 0.9 values; other multi-chain counts use
#   evenly spaced values from 0.3 to 0.9. Dirichlet chains start from independent
#   draws on the composition simplex. Other latent quantities are left to JAGS.
#
# Returns:
#   A list containing one chain-specific initial-value list per requested chain.
scenario_initial_values <- function(prior, nweek, nstock, n_chains = 3L) {
  prior <- match.arg(prior, c("ar1", "dirichlet"))
  n_chains <- as.integer(n_chains)
  if (length(n_chains) != 1L || is.na(n_chains) || n_chains <= 0L) {
    stop("n_chains must be a positive integer")
  }

  if (prior == "ar1") {
    phi_values <- if (n_chains == 3L) {
      c(0.5, 0.7, 0.9)
    } else if (n_chains == 1L) {
      0.7
    } else {
      seq(0.3, 0.9, length.out = n_chains)
    }
    return(lapply(phi_values, function(phi) list(phi = phi)))
  }

  if (!requireNamespace("MCMCpack", quietly = TRUE)) {
    stop("MCMCpack is required to generate Dirichlet initial values")
  }

  replicate(
    n_chains,
    list(p = MCMCpack::rdirichlet(nweek, rep(1, nstock))),
    simplify = FALSE
  )
}
