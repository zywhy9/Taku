#!/usr/bin/env bash
#SBATCH --no-requeue
set -euo pipefail

if [[ $# != 6 ]]; then
  echo "Use cluster/submit.sh rather than submitting this file directly." >&2
  exit 2
fi

MODE="$1"
REPOSITORY="$(realpath -- "$2")"
CONFIG="$(realpath -- "$3")"
INPUT_PATH="$4"
OUTPUT_PATH="$5"
SCENARIO="$6"

source "$CONFIG"

if (( ${#TAKU_MODULES[@]} > 0 )); then
  command -v module >/dev/null 2>&1 || {
    echo "TAKU_MODULES is set but the module command is unavailable." >&2
    exit 1
  }
  module load "${TAKU_MODULES[@]}"
fi

if [[ -n "${TAKU_R_LIBRARY:-}" ]]; then
  export R_LIBS_USER="$TAKU_R_LIBRARY"
  export RENV_CONFIG_AUTOLOADER_ENABLED=false
fi

cd "$REPOSITORY"

case "$MODE" in
  prepare)
    exec "$TAKU_RSCRIPT" cluster/run.R prepare \
      "$REPOSITORY" "$INPUT_PATH" "$OUTPUT_PATH"
    ;;
  fit)
    [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || {
      echo "A fit job must be submitted as a Slurm array." >&2
      exit 1
    }
    exec "$TAKU_RSCRIPT" cluster/run.R fit \
      "$REPOSITORY" "$OUTPUT_PATH" "$SCENARIO" \
      "$SLURM_ARRAY_TASK_ID"
    ;;
  summarize)
    exec "$TAKU_RSCRIPT" cluster/run.R summarize \
      "$REPOSITORY" "$OUTPUT_PATH"
    ;;
  *)
    echo "Unknown cluster job mode: $MODE" >&2
    exit 2
    ;;
esac
