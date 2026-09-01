#!/usr/bin/env bash

# Copy this file to cluster/config.sh and adapt it to the local Slurm system.
# Scientific settings are not defined here; all jobs use taku_paper_settings().

# Leave account or partition empty when the scheduler supplies a default.
TAKU_ACCOUNT=""
TAKU_PARTITION=""

# Optional environment modules. Examples differ across clusters, so no module
# name or version is assumed by the repository.
TAKU_MODULES=()

# Rscript used after the modules above are loaded. Leave TAKU_R_LIBRARY empty
# to activate the project renv library. Setting an external library disables
# renv autoload for the cluster commands so that the requested library is used.
TAKU_RSCRIPT="Rscript"
TAKU_R_LIBRARY=""

# Each array task fits one independent dataset with one R process. Increase
# CPUs per task only if the local JAGS/R installation is deliberately threaded.
TAKU_CPUS_PER_TASK=1
TAKU_MEMORY_PER_TASK="4G"
TAKU_WALLTIME_PER_DATASET="24:00:00"

# Slurm creates 1,000 tasks for each of the four scenarios but runs no more
# than this many at once. Choose the value with the cluster allocation and fair
# use policy in mind.
TAKU_MAX_CONCURRENT=50

# Small preparation and summary jobs do not run MCMC.
TAKU_PREPARE_MEMORY="2G"
TAKU_PREPARE_TIME="00:30:00"
TAKU_SUMMARY_MEMORY="4G"
TAKU_SUMMARY_TIME="00:30:00"
