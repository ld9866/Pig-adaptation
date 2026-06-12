#!/bin/bash
# =============================================================================
# 02_featurecounts.sh
#
# Step 2: per-sample, gene-level read quantification with featureCounts.
#
# This reproduces the original command:
#   featureCounts -p -T 50 -t exon -g gene_id \
#       -a Sus_scrofa.Sscrofa11.1.111.chr.gtf \
#       -o ./05.result/<i>.counts.txt  ./04.bam/<i>.group.sorted.dedup.bam
# run once per sample (the originals were launched concurrently with '&').
#
#   -p           : paired-end input. NOTE: in subread >= 2.0.2, -p only
#                  declares paired-end input; counting fragments requires
#                  --countReadPairs. The original run used -p alone, so keep it
#                  identical and record the subread version in the README.
#   -t exon      : count over 'exon' features.
#   -g gene_id   : aggregate to gene level -> row IDs are Ensembl gene IDs
#                  (ENSSSCG...). No --extraAttributes were used, so each output
#                  is the standard 7-column table (Geneid..Length, then count).
#   (no -Q)      : no mapping-quality filter (matches the .summary, where
#                  Unassigned_MappingQuality = 0 while multimappers remain).
#
# INPUTS : ${BAM_DIR}/<i>.group.sorted.dedup.bam   for each <i> in $LIST_FILE
#          Ensembl GTF given by $GTF_FILE
# OUTPUT : ${OUTPUT_DIR}/<i>.counts.txt (+ .summary) per sample
#
# USAGE  : bash 02_featurecounts.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GTF_FILE="/home/lidong/Reference/02.Ensemble.Duroc/Sus_scrofa.Sscrofa11.1.111.chr.gtf"
BAM_DIR="${SCRIPT_DIR}/04.bam"
OUTPUT_DIR="${SCRIPT_DIR}/05.result"
LOG_DIR="${SCRIPT_DIR}/logs"
LIST_FILE="${SCRIPT_DIR}/samplelist"
THREADS=16   # the original run used -T 50

mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"

while IFS= read -r i; do
    [[ -z "${i}" ]] && continue
    bam="${BAM_DIR}/${i}.group.sorted.dedup.bam"
    if [[ ! -f "${bam}" ]]; then
        echo "WARNING: missing BAM ${bam}, skipping." >&2
        continue
    fi
    echo "[$(date '+%H:%M:%S')] featureCounts: ${i}"
    featureCounts -p -T "${THREADS}" \
        -t exon \
        -g gene_id \
        -a "${GTF_FILE}" \
        -o "${OUTPUT_DIR}/${i}.counts.txt" \
        "${bam}" \
        1> "${LOG_DIR}/${i}.featureCounts.log" 2>&1
done < "${LIST_FILE}"

# To reproduce the original parallel launch instead, replace the loop body's
# featureCounts call with a backgrounded one ('... &') and add 'wait' after the
# loop. Mind total load: -T 50 x 4 samples = 200 threads.

echo "Done. Per-sample count tables in ${OUTPUT_DIR}/"
