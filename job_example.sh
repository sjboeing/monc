#!/bin/bash --login
#PBS -N monc_craymods
#PBS -l select=8:coretype=milan:ncpus=128:mpiprocs=128
#PBS -l walltime=03:00:00
#PBS -q collab
#PBS -j oe

########################################################################
## Strict mode
########################################################################

set -Eeuo pipefail
cd "$PBS_O_WORKDIR"

########################################################################
## Job configuration (EDIT HERE WHEN COPYING SCRIPT)
########################################################################

PROJECT_ROOT="$PBS_O_WORKDIR"

# MONC / scientific configuration
MONC_EXEC="$PROJECT_ROOT/build/bin/monc_driver.exe"
TESTCASE="$PROJECT_ROOT/testcases/arm_diurnal_cycle/arm_sb.mcf"

# Continuation control
MAX_CHAIN_LENGTH=99

# Case‑specific run directory root
RUNS_BASE="$PROJECT_ROOT/runs"

########################################################################
## Environment
########################################################################

#source /projects/parachute/steven.boeing.ext/cce15-stack/env_monc_craymods.sh
module load  cray-hdf5-parallel
module load cray-netcdf-hdf5parallel
ulimit -c unlimited
ulimit -s unlimited

# These should help to stabilise restarts
export MPICH_USE_SHM=0
export MPICH_SMP_SINGLE_COPY_MODE=NONE
export MPICH_RANK_REORDER_METHOD=0
export MPICH_CPU_BINDING=1

########################################################################
## Chain state (PBS-safe)
########################################################################

CHAIN_INDEX="${CHAIN_INDEX:-0}"
INPUT_CP_FILE="${INPUT_CP_FILE:-}"

echo "[INFO] Chain index        : $CHAIN_INDEX"
echo "[INFO] Input checkpoint  : ${INPUT_CP_FILE:-<none>}"
echo "[INFO] Max chain length  : $MAX_CHAIN_LENGTH"

if [[ "$CHAIN_INDEX" -ge "$MAX_CHAIN_LENGTH" ]]; then
  echo "[INFO] Max chain length reached — stopping"
  exit 0
fi

########################################################################
## Job-scoped run directory (case‑scoped namespace)
########################################################################

JOB_ID="$PBS_JOBID"
RUN_DIR="$RUNS_BASE/$JOB_ID"
mkdir -p "$RUN_DIR"

OUTFILE="$RUN_DIR/output.log"
RUN_CHECKPOINT="$RUN_DIR/checkpoint.nc"

echo "[INFO] Job ID        : $JOB_ID"
echo "[INFO] Run directory : $RUN_DIR"
echo "[INFO] Output log    : $OUTFILE"
echo "[INFO] Run checkpoint: $RUN_CHECKPOINT"

########################################################################
## Extract live checkpoint path from MCF (authoritative during run)
########################################################################

CHECKPOINT_REL=$(python3 "$PROJECT_ROOT/get_checkpoint_from_mcf.py" "$TESTCASE")
if [[ -z "$CHECKPOINT_REL" ]]; then
  echo "[FATAL] Failed to extract checkpoint_file from MCF"
  exit 1
fi

LIVE_CHECKPOINT="$PROJECT_ROOT/$CHECKPOINT_REL"
LIVE_CKPT_DIR="$(dirname "$LIVE_CHECKPOINT")"
mkdir -p "$LIVE_CKPT_DIR"

echo "[INFO] Live checkpoint path (MCF): $LIVE_CHECKPOINT"

########################################################################
## PBS diagnostics and config sanity
########################################################################

NP=$(wc -l < "$PBS_NODEFILE")
NODES=$(sort "$PBS_NODEFILE" | uniq | wc -l)

echo "[INFO] MPI ranks (PBS) : $NP"
echo "[INFO] Nodes           : $NODES"

if [[ "$NP" -le 1 ]]; then
  echo "[FATAL] PBS_NODEFILE reports NP <= 1"
  exit 1
fi

if [[ ! -f "$TESTCASE" ]]; then
  echo "[FATAL] MONC config not found: $TESTCASE"
  exit 1
fi

########################################################################
## Restart input sanity
########################################################################

if [[ -n "$INPUT_CP_FILE" ]]; then
  if [[ ! -s "$INPUT_CP_FILE" ]]; then
    echo "[FATAL] Restart requested but input checkpoint is missing or empty:"
    echo "  $INPUT_CP_FILE"
    exit 1
  fi
fi

########################################################################
## Run MONC
########################################################################

if [[ -n "$INPUT_CP_FILE" ]]; then
  echo "[INFO] Restarting from checkpoint: $INPUT_CP_FILE"
  mpirun "$MONC_EXEC" \
    --checkpoint="$INPUT_CP_FILE" \
    >"$OUTFILE" 2>&1
else
  echo "[INFO] Cold start from config"
  mpirun "$MONC_EXEC" \
    --config="$TESTCASE" \
    >"$OUTFILE" 2>&1
fi

MONC_RC=$?

if [[ "$MONC_RC" -ne 0 ]]; then
  echo "[FATAL] MONC exited with code $MONC_RC"
  exit 1
fi

########################################################################
## Verify output log
########################################################################

if [[ ! -s "$OUTFILE" ]]; then
  echo "[FATAL] Output log was not created or is empty:"
  echo "  $OUTFILE"
  exit 1
fi

########################################################################
## Verify live checkpoint was produced, then MOVE it into run dir
##
## PROGRESS INVARIANT:
##   A run made progress iff it produced a checkpoint file
########################################################################

if [[ ! -f "$LIVE_CHECKPOINT" ]]; then
  echo "[INFO] No checkpoint produced by this run"
  echo "[INFO] Chain will not continue"
  exit 0
fi

# Atomically move checkpoint into run directory (ownership transfer)
mv "$LIVE_CHECKPOINT" "$RUN_CHECKPOINT"

echo "[INFO] Checkpoint captured:"
echo "       $RUN_CHECKPOINT"

########################################################################
## Scientific termination decision
########################################################################

RESULT=$(python3 "$PROJECT_ROOT/check_continuation.py" "$OUTFILE")
echo "[INFO] Continuation decision: $RESULT"

if [[ "$RESULT" = "continue" ]]; then
  NEXT_CHAIN_INDEX=$((CHAIN_INDEX + 1))
  echo "[INFO] Continuation required (next index = $NEXT_CHAIN_INDEX)"

  "$PROJECT_ROOT/submit_next.sh" \
    "$0" \
    "$RUN_CHECKPOINT" \
    "$NEXT_CHAIN_INDEX"
else
  echo "[INFO] Scientific termination reached — chain ends"
fi
