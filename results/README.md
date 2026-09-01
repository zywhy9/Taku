# Reference results

The two tracked tables report the 1,000-replicate baseline simulation. Individual
simulation datasets, fitted objects, and posterior draws are not included.

- `simulation_performance.csv` reports relative bias, relative RMSE, empirical
  95% interval coverage, mean interval length, and mean computation time per
  dataset for the four Bayesian scenarios and three frequentist MoM
  implementations.
- `simulation_convergence.csv` reports the number and proportion of Bayesian
  fits satisfying `R-hat <= 1.10` and ESS at least 100 for the stopping
  targets (`N`, plus `phi` for AR(1) models).

Running `run_simulation.R` writes tables with the same schema under the selected
output directory, so the tracked reference tables can be compared directly.
