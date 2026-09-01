#!/usr/bin/env bash
set -euo pipefail

REPOSITORY="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG="$REPOSITORY/cluster/config.sh"

die() {
  echo "ERROR: $*" >&2
  exit 1
}

[[ -f "$CONFIG" ]] || die \
  "copy cluster/config.example.sh to cluster/config.sh and edit it first"
source "$CONFIG"

for variable in TAKU_RSCRIPT TAKU_CPUS_PER_TASK TAKU_MEMORY_PER_TASK \
  TAKU_WALLTIME_PER_DATASET TAKU_MAX_CONCURRENT TAKU_PREPARE_MEMORY \
  TAKU_PREPARE_TIME TAKU_SUMMARY_MEMORY TAKU_SUMMARY_TIME; do
  [[ -n "${!variable:-}" ]] || die "missing setting in config.sh: $variable"
done
[[ "$TAKU_CPUS_PER_TASK" =~ ^[1-9][0-9]*$ ]] || die \
  "TAKU_CPUS_PER_TASK must be a positive integer"
[[ "$TAKU_MAX_CONCURRENT" =~ ^[1-9][0-9]*$ ]] || die \
  "TAKU_MAX_CONCURRENT must be a positive integer"

scheduler_options=()
[[ -z "${TAKU_ACCOUNT:-}" ]] || scheduler_options+=(--account="$TAKU_ACCOUNT")
[[ -z "${TAKU_PARTITION:-}" ]] || scheduler_options+=(--partition="$TAKU_PARTITION")

run_r() {
  if (( ${#TAKU_MODULES[@]} > 0 )); then
    command -v module >/dev/null 2>&1 || die \
      "TAKU_MODULES is set but the module command is unavailable"
    module load "${TAKU_MODULES[@]}"
  fi
  if [[ -n "${TAKU_R_LIBRARY:-}" ]]; then
    export R_LIBS_USER="$TAKU_R_LIBRARY"
    export RENV_CONFIG_AUTOLOADER_ENABLED=false
  fi
  (cd "$REPOSITORY" && "$TAKU_RSCRIPT" "$@")
}

active_job_ids() {
  local jobs_directory="$1"
  local record job_id
  [[ -d "$jobs_directory" ]] || return 0
  for record in "$jobs_directory"/*.jobid; do
    [[ -f "$record" ]] || continue
    job_id="$(<"$record")"
    [[ "$job_id" =~ ^[0-9]+$ ]] || die "invalid job ID in $record"
    if squeue --noheader --jobs="$job_id" 2>/dev/null | grep -q .; then
      printf '%s\n' "$job_id"
    fi
  done
}

command="${1:-}"
case "$command" in
  check)
    [[ $# == 3 ]] || die \
      "usage: bash cluster/submit.sh check INPUT.RData OUTPUT_DIRECTORY"
    run_r cluster/run.R check "$REPOSITORY" "$2" "$3"
    ;;
  prepare)
    [[ $# == 3 ]] || die \
      "usage: bash cluster/submit.sh prepare INPUT.RData OUTPUT_DIRECTORY"
    input="$(realpath -- "$2")"
    output="$(realpath -m -- "$3")"
    mkdir -p "$output/logs" "$output/jobs"
    job_id="$(sbatch --parsable "${scheduler_options[@]}" \
      --job-name=taku-prepare --cpus-per-task=1 \
      --mem="$TAKU_PREPARE_MEMORY" --time="$TAKU_PREPARE_TIME" \
      --output="$output/logs/prepare-%j.out" \
      "$REPOSITORY/cluster/job.sbatch" prepare "$REPOSITORY" "$CONFIG" \
      "$input" "$output" "")"
    job_id="${job_id%%;*}"
    printf '%s\n' "$job_id" > "$output/jobs/prepare.jobid"
    echo "Submitted preparation job: $job_id"
    ;;
  run)
    [[ $# == 2 ]] || die \
      "usage: bash cluster/submit.sh run OUTPUT_DIRECTORY"
    output="$(realpath -- "$2")"
    [[ -f "$output/settings.rds" && \
       -f "$output/simdata/baseline/manifest.rds" ]] || die \
      "prepared data are missing; wait for the prepare job"
    mapfile -t active_jobs < <(active_job_ids "$output/jobs")
    (( ${#active_jobs[@]} == 0 )) || die \
      "jobs are already active for this output: ${active_jobs[*]}"
    n_datasets="$(run_r -e \
      'cat(readRDS(commandArgs(TRUE)[1])$n_datasets)' \
      "$output/settings.rds")"
    [[ "$n_datasets" =~ ^[1-9][0-9]*$ ]] || die \
      "could not read n_datasets from settings.rds"

    scenarios=(rdm_ar rdm_dirichlet mmd_ar mmd_dirichlet)
    fit_jobs=()
    for scenario in "${scenarios[@]}"; do
      job_id="$(sbatch --parsable "${scheduler_options[@]}" \
        --job-name="taku-$scenario" \
        --array="1-${n_datasets}%${TAKU_MAX_CONCURRENT}" \
        --cpus-per-task="$TAKU_CPUS_PER_TASK" \
        --mem="$TAKU_MEMORY_PER_TASK" \
        --time="$TAKU_WALLTIME_PER_DATASET" \
        --output="$output/logs/${scenario}-%A_%a.out" \
        "$REPOSITORY/cluster/job.sbatch" fit "$REPOSITORY" "$CONFIG" \
        "" "$output" "$scenario")"
      job_id="${job_id%%;*}"
      fit_jobs+=("$job_id")
      printf '%s\n' "$job_id" > "$output/jobs/${scenario}.jobid"
      echo "Submitted $scenario array: $job_id"
    done

    dependency="$(IFS=:; echo "${fit_jobs[*]}")"
    summary_job="$(sbatch --parsable "${scheduler_options[@]}" \
      --dependency="afterok:$dependency" --job-name=taku-summary \
      --cpus-per-task=1 --mem="$TAKU_SUMMARY_MEMORY" \
      --time="$TAKU_SUMMARY_TIME" \
      --output="$output/logs/summary-%j.out" \
      "$REPOSITORY/cluster/job.sbatch" summarize "$REPOSITORY" "$CONFIG" \
      "" "$output" "")"
    summary_job="${summary_job%%;*}"
    printf '%s\n' "$summary_job" > "$output/jobs/summary.jobid"
    echo "Submitted dependent summary job: $summary_job"
    ;;
  summarize)
    [[ $# == 2 ]] || die \
      "usage: bash cluster/submit.sh summarize OUTPUT_DIRECTORY"
    output="$(realpath -- "$2")"
    mkdir -p "$output/logs" "$output/jobs"
    mapfile -t active_jobs < <(active_job_ids "$output/jobs")
    (( ${#active_jobs[@]} == 0 )) || die \
      "jobs are already active for this output: ${active_jobs[*]}"
    job_id="$(sbatch --parsable "${scheduler_options[@]}" \
      --job-name=taku-summary --cpus-per-task=1 \
      --mem="$TAKU_SUMMARY_MEMORY" --time="$TAKU_SUMMARY_TIME" \
      --output="$output/logs/summary-%j.out" \
      "$REPOSITORY/cluster/job.sbatch" summarize "$REPOSITORY" "$CONFIG" \
      "" "$output" "")"
    job_id="${job_id%%;*}"
    printf '%s\n' "$job_id" > "$output/jobs/summary.jobid"
    echo "Submitted summary job: $job_id"
    ;;
  status)
    [[ $# == 2 ]] || die \
      "usage: bash cluster/submit.sh status OUTPUT_DIRECTORY"
    output="$(realpath -- "$2")"
    squeue --user="$(id -un)" \
      --format='%.18i %.28j %.10T %.12M %.12l %R'
    echo "Completed records by scenario:"
    for scenario in rdm_ar rdm_dirichlet mmd_ar mmd_dirichlet; do
      count=0
      [[ ! -d "$output/records/$scenario" ]] || \
        count="$(find "$output/records/$scenario" -maxdepth 1 \
          -type f -name 'data*.rds' | wc -l)"
      printf '  %-16s %s\n' "$scenario" "$count"
    done
    ;;
  *)
    die "usage: bash cluster/submit.sh check|prepare|run|status|summarize ..."
    ;;
esac
