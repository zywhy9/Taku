# Taku compositional-uncertainty analysis

This repository contains the code for the baseline simulation study and the
Taku salmon data application. It compares four Bayesian model combinations
with three frequentist method-of-moments (MoM) estimators.

## Start here

There are only three analysis entry points:

| File | Purpose |
| --- | --- |
| `run_simulation.R` | Generate the 1,000 baseline datasets, fit the four Bayesian scenarios, and summarize all seven methods. |
| `run_application.R` | Fit all seven methods to the Taku data. |
| `make_figures.R` | Recreate Figures 1 and 6--8 from a completed simulation run. |

Both scripts load `R/taku.R`. The files under `R/internal/` contain the tested
implementation; they do not need to be run individually.

## Software

The analysis was run with R 4.6.1, JAGS 4.3.2, and the package versions in
`renv.lock`.

```sh
Rscript -e "install.packages('renv')"
Rscript -e "renv::restore()"
Rscript tests/test_core.R
Rscript tests/test_application.R
Rscript tests/test_jags_smoke.R
```

The first two tests use synthetic data and do not run MCMC. The JAGS smoke test
uses deliberately short chains to check model compilation only.

## Private input

The empirical data are not distributed. To reproduce the analysis, supply one
`.RData` file containing:

- `n`: weekly sample sizes;
- `w`: weekly fish-wheel weights;
- `region1mean`: week-by-reporting-unit estimated percentages;
- `region1sd`: corresponding standard deviations in percentage points.

The four columns must be ordered as Taku River, Other Lakes, Tatsamenie Lake,
and U.S. River. Other Lakes and Tatsamenie Lake form the lake-stock proportion
in the abundance denominator. Input data and generated analysis directories are
excluded by `.gitignore`.

## Simulation study

From the repository root, run:

```sh
Rscript run_simulation.R path/to/region1.RData path/to/output
```

This command generates the baseline datasets once, fits all four Bayesian
scenarios, evaluates the three frequentist estimators, and writes the summary
tables under `path/to/output/results/`. Existing per-dataset records are reused,
so an interrupted run can be resumed with the same command.

Independent datasets can be fitted in parallel on one computer by supplying a
scenario and worker count:

```sh
Rscript run_simulation.R path/to/region1.RData path/to/output rdm_ar 4
```

Valid scenario IDs are `rdm_ar`, `rdm_dirichlet`, `mmd_ar`, and
`mmd_dirichlet`; use `all` to run all four. The default is one worker because a
safe worker count depends on the reader's available memory and CPU allocation.

The paper settings are defined once in `taku_paper_settings()`:

- 1,000 datasets and seed 1234;
- three chains;
- 100,000 adaptation iterations;
- 50,000 retained draws per diagnostic batch and chain;
- stopping targets `N`, plus `phi` for AR(1) models;
- `R-hat <= 1.10` and ESS at least 100;
- maximum 1,000,000,000 iterations or 1,320 minutes per dataset.

Weekly composition posterior means are retained as a compact 12-by-4 summary
for each fit. They are output only: they do not enter the stopping rule, and the
full posterior draws are not saved.

The generated data use the RDM data-generating process. Weekly multinomial
sample counts are converted to the IA summaries used by every fitted method.
All scenarios reuse exactly the same generated datasets.

After all four scenarios finish, create the four 12-by-4 boxplot figures:

```sh
Rscript make_figures.R path/to/output path/to/figure-output
```

### Optional cluster execution

The commands above are the default workflow and can run sequentially or with
multiple local workers. For larger runs on a Slurm system, the optional
[`cluster/`](cluster/) wrapper submits one independent dataset per array task
while calling the same R functions and paper settings. It is scheduler-generic:
account, partition, modules, resources, walltime, and maximum concurrency are
defined in an untracked local config file rather than hard-coded for one named
cluster.

## Taku data application

Supply the private input, the observed lake-type count, and a new output
directory:

```sh
Rscript run_application.R path/to/region1.RData 34351 path/to/application-output
```

Replace `34351` if reproducing an application with a different observed count.
The script writes one compact CSV plus one RDS fit for each Bayesian scenario.
It refuses to overwrite a non-empty output directory.

## Models and estimators

The Bayesian scenarios combine either the reverse Dirichlet-multinomial model
(RDM) or the moment-matching Dirichlet approximation (MMD) with either an AR(1)
logistic-normal prior or an independent Dirichlet prior.

The three frequentist methods share the same MoM abundance estimate but use
different variance calculations:

- `MoM` derives marginal variances and covariance from the population-level
  coefficient;
- `MoM(Alt)` reconstructs the combined lake variance at the sample level and
  then inflates that variance to the population level;
- `MoM(Naive)` follows the original Gazey calculation by summing the two
  reported marginal variances and treating the lake stocks as independent.

## Repository layout

```text
run_simulation.R       complete simulation workflow
run_application.R      complete data-application workflow
make_figures.R          Figures 1 and 6--8
R/taku.R               single function loader
R/internal/             model and estimator implementation
cluster/                optional generic Slurm wrapper
tests/                  synthetic regression and JAGS smoke tests
results/                aggregate reference results
renv.lock               exact R package versions
```

The tracked tables in `results/` contain aggregate results only. Private data,
generated datasets, individual fit records, and posterior draws are not
included.
