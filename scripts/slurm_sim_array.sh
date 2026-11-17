#!/bin/bash
#SBATCH --job-name=bw_sim_bench
#SBATCH --partition=all
#SBATCH --array=0-179
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --output=logs/sim_%A_%a.out
#SBATCH --error=logs/sim_%A_%a.err

set -euo pipefail

# Repository checkout to run the worker from. Defaults to the repo containing
# this script; override at submit time with, e.g.:
#   sbatch --export=ALL,BUSHWHACKER_REPO=/path/to/checkout scripts/slurm_sim_array.sh
# The %A/%a log paths above are relative to the directory the job runs in.
REPO_DIR="${BUSHWHACKER_REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$REPO_DIR"
mkdir -p logs

export RAYON_NUM_THREADS=4
python3 scripts/benchmarks/run_silverback_sim_worker.py "${SLURM_ARRAY_TASK_ID}"
