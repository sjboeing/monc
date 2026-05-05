#!/bin/bash
set -Eeuo pipefail

# Usage:
#   submit_next.sh <job_script> <checkpoint_file> <chain_index>

# The content of this file has been produced with the assistance of
# University of Leeds, m365 co-pilot
#
JOB_SCRIPT="$1"
CHECKPOINT_FILE="$2"
CHAIN_INDEX="$3"

if [[ -z "${JOB_SCRIPT:-}" || -z "${CHECKPOINT_FILE:-}" || -z "${CHAIN_INDEX:-}" ]]; then
  echo "[FATAL] Usage: submit_next.sh <job_script> <checkpoint_file> <chain_index>"
  exit 1
fi

if [[ ! -f "$JOB_SCRIPT" ]]; then
  echo "[FATAL] Job script not found: $JOB_SCRIPT"
  exit 1
fi

if [[ ! -f "$CHECKPOINT_FILE" ]]; then
  echo "[FATAL] submit_next.sh: checkpoint missing: $CHECKPOINT_FILE"
  exit 1
fi

echo "[INFO] Submitting continuation job"
echo " Chain index : $CHAIN_INDEX"
echo " Checkpoint  : $CHECKPOINT_FILE"
echo " Job script  : $JOB_SCRIPT"

if [[ ! -x "$JOB_SCRIPT" ]]; then
  echo "[FATAL] Job script is not executable: $JOB_SCRIPT"
  exit 1
fi

qsub -v \
  INPUT_CP_FILE="$CHECKPOINT_FILE",CHAIN_INDEX="$CHAIN_INDEX" \
  "$JOB_SCRIPT"
