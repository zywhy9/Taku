# Running the simulation on a Slurm cluster

The repository is designed to run locally by default. This directory contains
an optional, scheduler-only layer for clusters using Slurm. It does not contain
a second copy of the models, estimators, simulation generator, or scientific
settings.

The cluster workflow uses four Slurm job arrays, one per Bayesian scenario.
Each array task fits one independent dataset with one R process. The array
concurrency limit prevents all 4,000 fits from starting simultaneously and can
be adapted to the user's allocation and local fair-use policy.

## 1. Configure the environment

Clone the repository on a filesystem visible to compute nodes, restore the R
environment, and create the local scheduler configuration:

```sh
Rscript -e "install.packages('renv')"
Rscript -e "renv::restore()"
cp cluster/config.example.sh cluster/config.sh
```

Edit only `cluster/config.sh`. At minimum, review the account, partition,
module names, memory, walltime, and `TAKU_MAX_CONCURRENT`. The local config is
ignored by Git.

Leave `TAKU_R_LIBRARY` empty after `renv::restore()`. If packages are instead
installed in a cluster-specific user library, set its absolute path there; the
wrapper will disable renv autoload for those jobs and use that library.

The defaults request one CPU and 4 GB per dataset for up to 24 hours. Those are
portable starting values, not universal recommendations. Some clusters allocate
whole nodes, restrict array size, impose shorter walltimes, or require a
different memory syntax. Follow the local scheduler documentation and reduce
the concurrency limit when the allocation cannot support the requested number
of simultaneous tasks.

## 2. Check and prepare the shared datasets

The check is lightweight and does not run MCMC:

```sh
bash cluster/submit.sh check path/to/region1.RData /shared/path/taku-run
```

Submit one small job to generate and validate the 1,000 shared datasets:

```sh
bash cluster/submit.sh prepare path/to/region1.RData /shared/path/taku-run
```

Wait for that job to finish successfully before fitting. Its ID is stored in
`/shared/path/taku-run/jobs/prepare.jobid`.

## 3. Submit the Bayesian fits

```sh
bash cluster/submit.sh run /shared/path/taku-run
```

This submits four arrays and one dependent summary job. Dataset identity and
MCMC seeds depend only on the fixed paper seed, scenario, and dataset ID; they
do not depend on node, array job ID, or scheduling order. Existing completed
records are reused, so the command can be run again after failed tasks have left
the queue. The submission script refuses to create a duplicate batch while any
job recorded for the same output directory is still active.

Monitor without creating additional jobs:

```sh
bash cluster/submit.sh status /shared/path/taku-run
```

If the dependent summary did not run because one or more array tasks failed,
resubmit the fits after diagnosing the scheduler logs. Once all records exist,
the summary alone can be submitted with:

```sh
bash cluster/submit.sh summarize /shared/path/taku-run
```

Final tables are written below `/shared/path/taku-run/results/`. The same
`make_figures.R` used locally recreates Figures 1 and 6--8.
