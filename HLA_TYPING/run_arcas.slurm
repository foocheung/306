#!/bin/bash
#SBATCH --job-name=arcasHLA
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=72:00:00
#SBATCH --array=0-6
#SBATCH --output=logs/arcasHLA_%A_%a.out
#SBATCH --error=logs/arcasHLA_%A_%a.err

set -euo pipefail

module load arcashla

mkdir -p logs arcasHLA_results

SAMPLES=(
  Manthiram_T10_S39
  Manthiram_T22_S44
  Manthiram_T70_S45
  Manthiram_T101_S40
  Manthiram_T103_S41
  Manthiram_T126_S42
  Manthiram_T133_S43
)

SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

R1=${SAMPLE}_R1_001.fastq.gz
R2=${SAMPLE}_R2_001.fastq.gz

echo "Running arcashla on ${SAMPLE}"

arcasHLA genotype \
"$R1" "$R2" \
-g A,B,C,DQA1,DQB1,DRB1 \
-t "${SLURM_CPUS_PER_TASK}" \
-v \
-o "arcasHLA_results/${SAMPLE}"

