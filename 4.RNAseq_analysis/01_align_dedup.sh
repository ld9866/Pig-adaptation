#!/bin/bash
# =============================================================================
# 01_align_dedup.sh
#
# Step 1 (TEMPLATE): align cleaned paired-end RNA-seq reads to the Sus scrofa
# reference, add read groups, sort, index, and mark duplicates, producing
#   04.bam/<i>.group.sorted.dedup.bam   (the input expected by step 2).
#
# NOTE: this is a representative HISAT2-based template. Replace the aligner /
# index block with the exact command you used and confirm it reproduces your
# BAMs before publishing.
#
#   cleaned FASTQ  --HISAT2-->  SAM  --samtools sort-->  sorted BAM
#                  --samtools index-->  .bai
#                  --GATK MarkDuplicates-->  deduplicated BAM (+ metrics)
#
# INPUTS  (per sample <i>, listed in $LIST_FILE):
#   ${INPUT_DIR}/<i>.clean_1.fq.gz
#   ${INPUT_DIR}/<i>.clean_2.fq.gz
#   HISAT2 index given by $INDEX
#
# OUTPUTS:
#   ${OUTPUT_DIR}/02.sorted_bam/<i>.sorted.bam(.bai)
#   ${OUTPUT_DIR}/04.bam/<i>.group.sorted.dedup.bam(.bai) + metrics
#   ${OUTPUT_DIR}/logs/<i>.hisat2.log, <i>.gatk.log
#
# USAGE:
#   bash 01_align_dedup.sh
# =============================================================================

set -euo pipefail

# ----------------------------- configuration ---------------------------------
# Edit these paths to match your environment.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"  # repo root
INPUT_DIR="${SCRIPT_DIR}/01.cleanfq"          # cleaned FASTQ files
OUTPUT_DIR="${SCRIPT_DIR}"                     # where 02/03/logs are written
INDEX="${SCRIPT_DIR}/00.index/Sus_scrofa.Sscrofa11.1"  # HISAT2 index prefix
LIST_FILE="${SCRIPT_DIR}/samplelist"           # one sample ID per line
THREADS=16

# ----------------------------- preparation -----------------------------------
mkdir -p "${OUTPUT_DIR}/02.sorted_bam" \
         "${OUTPUT_DIR}/04.bam" \
         "${OUTPUT_DIR}/logs"

if [[ ! -f "${LIST_FILE}" ]]; then
    echo "ERROR: sample list not found: ${LIST_FILE}" >&2
    exit 1
fi

echo "=========================================="
echo " Input  : ${INPUT_DIR}"
echo " Output : ${OUTPUT_DIR}"
echo " Index  : ${INDEX}"
echo "=========================================="

# ----------------------------- main loop -------------------------------------
while IFS= read -r i; do
    [[ -z "${i}" ]] && continue
    echo "[$(date '+%H:%M:%S')] Processing sample: ${i}"

    FQ1="${INPUT_DIR}/${i}.clean_1.fq.gz"
    FQ2="${INPUT_DIR}/${i}.clean_2.fq.gz"

    if [[ ! -f "${FQ1}" || ! -f "${FQ2}" ]]; then
        echo "WARNING: missing FASTQ for ${i}, skipping." >&2
        continue
    fi

    # --- HISAT2 alignment piped directly into coordinate sorting -------------
    # --dta            : report alignments tailored for transcript assemblers
    # --rg-id / --rg   : read-group tags (required by GATK downstream)
    hisat2 -p "${THREADS}" \
           -x "${INDEX}" \
           -1 "${FQ1}" \
           -2 "${FQ2}" \
           --dta \
           --rg-id "${i}" \
           --rg "SM:${i}" \
           --rg "LB:lib1" \
           --rg "PL:illumina" \
           2> "${OUTPUT_DIR}/logs/${i}.hisat2.log" \
    | samtools sort -@ "${THREADS}" \
           -o "${OUTPUT_DIR}/02.sorted_bam/${i}.group.sorted.bam" -

    samtools index "${OUTPUT_DIR}/02.sorted_bam/${i}.group.sorted.bam"

    # --- mark duplicates -----------------------------------------------------
    gatk MarkDuplicates \
        -I "${OUTPUT_DIR}/02.sorted_bam/${i}.group.sorted.bam" \
        -O "${OUTPUT_DIR}/04.bam/${i}.group.sorted.dedup.bam" \
        -M "${OUTPUT_DIR}/04.bam/${i}.dedup.metrics.txt" \
        --CREATE_INDEX true \
        --VALIDATION_STRINGENCY SILENT \
        1> "${OUTPUT_DIR}/logs/${i}.gatk.log" 2>&1

    echo "Sample ${i} done."
done < "${LIST_FILE}"

echo "All samples processed."
