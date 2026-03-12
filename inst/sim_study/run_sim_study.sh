#!/bin/bash
#SBATCH --job-name=PPDisentangle_Sim
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=128G

set -e

# ---- Resolve paths ----
if [ -z "$PP_SCRIPT_DIR" ]; then
    PP_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
if [ -z "$PP_PKG_ROOT" ]; then
    PP_PKG_ROOT="$(cd "$PP_SCRIPT_DIR/../.." && pwd)"
fi

PP_TEST="${PP_TEST:-}"
PP_SIMS="${PP_SIMS:-}"
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --sims) PP_SIMS="$2"; shift 2 ;;
        --test) PP_TEST="--test"; shift ;;
        *) shift ;;
    esac
done
if [ -n "$PP_TEST" ]; then
    PP_SIMS="${PP_SIMS:-2}"
else
    PP_SIMS="${PP_SIMS:-100}"
fi
PP_CPUS="$PP_SIMS"

# ---- If on login node, submit via sbatch and exit ----
if [ -z "$SLURM_JOB_ID" ]; then
    cd "$PP_PKG_ROOT"

    echo "Pulling latest code..."
    git pull origin main

    export PP_SCRIPT_DIR PP_PKG_ROOT PP_SIMS="$PP_SIMS" PP_CPUS="$PP_CPUS" PP_TEST="$PP_TEST"
    [ -n "$ATE_SEQUENTIAL" ] && export ATE_SEQUENTIAL
    [ -n "$PP_LOG_MEMORY" ] && export PP_LOG_MEMORY
    [ -n "$PP_SKIP_CRAZY_PARAMS" ] && export PP_SKIP_CRAZY_PARAMS

    SBATCH_EXPORT="PP_SCRIPT_DIR,PP_PKG_ROOT,PP_SIMS,PP_CPUS,PP_TEST"
    [ -n "$ATE_SEQUENTIAL" ] && SBATCH_EXPORT="${SBATCH_EXPORT},ATE_SEQUENTIAL"
    [ -n "$PP_LOG_MEMORY" ] && SBATCH_EXPORT="${SBATCH_EXPORT},PP_LOG_MEMORY"
    [ -n "$PP_SKIP_CRAZY_PARAMS" ] && SBATCH_EXPORT="${SBATCH_EXPORT},PP_SKIP_CRAZY_PARAMS"

    EXTRA_SBATCH=""
    if [ "$PP_CPUS" -gt 72 ]; then
        if [ "$PP_CPUS" -gt 256 ]; then
            echo "ERROR: max 256 CPUs per node (milan). Reduce --sims or split jobs."
            exit 1
        fi
        EXTRA_SBATCH="--partition=milan"
        echo "Note: >72 CPUs requested, using Milan partition"
    fi

    mkdir -p "$PP_PKG_ROOT/cluster_output"

    JOB_ID=$(sbatch --parsable \
        --cpus-per-task="$PP_CPUS" \
        $EXTRA_SBATCH \
        --output="$PP_PKG_ROOT/cluster_output/%j_slurm.out" \
        --error="$PP_PKG_ROOT/cluster_output/%j_slurm.err" \
        --export="$SBATCH_EXPORT" \
        "$PP_SCRIPT_DIR/run_sim_study.sh")

    echo "Job $JOB_ID submitted ($PP_SIMS sims, $PP_CPUS cores${PP_TEST:+, test mode})"
    echo "  Results: cluster_output/$JOB_ID.rds"
    echo "  Log:     cluster_output/$JOB_ID.log"
    echo "  SLURM:   cluster_output/${JOB_ID}_slurm.out"
    exit 0
fi

# ---- Inside SLURM job ----
cd "$PP_PKG_ROOT"
mkdir -p "$PP_PKG_ROOT/cluster_output"

echo "=== PPDisentangle Sim Study ==="
echo "Job ${SLURM_JOB_ID:-local} | $(date)"
echo "Sims: $PP_SIMS | Cores: $PP_CPUS"
echo ""

# R-Geo bundles R + GDAL + GEOS + PROJ + UDUNITS
module purge 2>/dev/null || true
module load R-Geo/4.3.2-foss-2023a 2>/dev/null || module load R 2>/dev/null || true

echo "Rscript: $(which Rscript)"
echo ""

Rscript "$PP_SCRIPT_DIR/sim_study.R" --cluster --sims "$PP_SIMS" $PP_TEST

echo "=== Done $(date) ==="
