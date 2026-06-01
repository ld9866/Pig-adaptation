#!/usr/bin/env bash
# ============================================================
# BayPass GEA Pipeline — Full Workflow
# Author : Lidong
# Version: 1.0
#
# Description:
#   End-to-end pipeline for BayPass genome-environment
#   association analysis. Covers:
#     Step 1  — Create directory structure
#     Step 2  — Split VCF by chromosome
#     Step 3  — Generate per-population sample lists
#     Step 4  — Subset VCF by population
#     Step 5  — Compute allele frequencies (AF/AC/AN)
#     Step 6  — Convert AF VCF to plain-text genotype files
#     Step 7  — Merge per-population txt into per-chromosome geno files
#     Step 8  — (Optional) Count columns in a geno file
#     Step 9  — Rename files containing underscores in population names
#     Step 10 — Extract unique population list from file names
#     Step 11 — Build eco covariate files and transpose
#     Step 12 — Install / test BayPass binary
#     Step 13 — Run BayPass covariate analysis (AnaCovis), then final run
#     Step 14 — Initialise per-eco merged result files (chr1 header)
#     Step 15 — Append chr2–18 results (header stripped) to merged files
#     Step 16 — Annotate merged results with SNP positions
#
# Usage:
#   bash baypass_pipeline.sh >> baypass_pipeline.log 2>&1
#
# Requirements:
#   - bcftools >= 1.13 on $PATH
#   - python + vcf2geno.py in working directory
#   - BayPass binary compiled at BAYPASS path below
#   - Input: 919qcldbeagle.vcf.gz + .tbi in working directory
#            919.sample.group.list  (col1=sampleID, col2=population)
#            temp_pop.txt           (one population name per line)
#            env.txt                (environmental variable table)
#            chr919qcldbeagle.bim   (PLINK BIM for SNP positions)
# ============================================================

set -euo pipefail

# ---- Configuration ---------------------------------------------------
BAYPASS="/home/Lidong/software/baypass_public-master/sources/g_baypass"
N_CHR=18
NTHREADS_COVIS=64
NTHREADS_FINAL=128


# ============================================================
# Step 1: Create directory structure
# ============================================================

echo "[Step 1] Creating directories..."

mkdir -p \
    00.original.data \
    01.group.list    \
    02.group.vcf     \
    03.group.af.vcf  \
    04.group.af.txt  \
    05.geno          \
    06.eco           \
    07.covresult     \
    08.omega         \
    09.final.result

echo "[Step 1] Done."


# ============================================================
# Step 2: Split merged VCF into per-chromosome files
# ============================================================
# bcftools view -r N restricts output to chromosome N.
# Output is bgzipped (.vcf.gz) and then indexed.
#
# NOTE: chromosome identifiers in the VCF must match the loop
# variable exactly (e.g. "1" not "chr1"). Check with:
#   bcftools view -h 919qcldbeagle.vcf.gz | grep "^##contig"

echo "[Step 2] Splitting VCF by chromosome..."

for i in $(seq 1 "$N_CHR"); do
    bcftools view \
        -r "$i" \
        919qcldbeagle.vcf.gz \
        -Oz -o "00.original.data/beagle.chr${i}.vcf.gz"
done

for i in $(seq 1 "$N_CHR"); do
    bcftools index "00.original.data/beagle.chr${i}.vcf.gz"
done

echo "[Step 2] Done."


# ============================================================
# Step 3: Generate per-population sample lists
# ============================================================
# Reads 919.sample.group.list (col1=sampleID, col2=population).
# For each population in temp_pop.txt, writes a file containing
# the sample IDs belonging to that population.

echo "[Step 3] Generating per-population sample lists..."

while IFS= read -r pop; do
    awk -v pop="$pop" '$2 == pop {print $1}' \
        919.sample.group.list \
        >> "01.group.list/group.${pop}.txt"
done < temp_pop.txt

echo "[Step 3] Done."


# ============================================================
# Step 4: Subset VCF by population
# ============================================================
# For each population × chromosome combination, extract only
# the samples listed in the corresponding group file.
# Jobs are submitted to the background; all finish before Step 5.
#
# NOTE: bcftools -S expects one sample ID per line with no
# extra whitespace. Verify your group list files accordingly.

echo "[Step 4] Subsetting VCF by population..."

while IFS= read -r pop; do
    for j in $(seq 1 "$N_CHR"); do
        nohup bcftools view \
            -S "./01.group.list/group.${pop}.txt" \
            "./00.original.data/beagle.chr${j}.vcf.gz" \
            -Oz -o "./02.group.vcf/group.${pop}.chr${j}.vcf.gz" \
            > "nohup_view_${pop}_chr${j}.out" 2>&1 &
    done
done < temp_pop.txt

wait
echo "[Step 4] Done."


# ============================================================
# Step 5: Compute allele frequency tags (AF, AC, AN)
# ============================================================
# bcftools +fill-tags recalculates AF/AC/AN from genotype data.
# This is required because subsetting removes samples and the
# original tags no longer reflect the subpopulation allele counts.

echo "[Step 5] Computing allele frequencies..."

while IFS= read -r pop; do
    for j in $(seq 1 "$N_CHR"); do
        nohup bcftools +fill-tags \
            "./02.group.vcf/group.${pop}.chr${j}.vcf.gz" \
            -Oz -o "./03.group.af.vcf/group.af.${pop}.chr${j}.vcf.gz" \
            -- -t AF,AC,AN \
            > "nohup_fill_${pop}_chr${j}.out" 2>&1 &
    done
done < temp_pop.txt

wait
echo "[Step 5] Done."


# ============================================================
# Step 6: Convert AF VCF to plain-text genotype files
# ============================================================
# vcf2geno.py extracts the AF field from each variant and
# writes allele counts per SNP as a space-separated text file.
# Ensure vcf2geno.py is in the working directory.

echo "[Step 6] Converting VCF to geno text files..."

while IFS= read -r pop; do
    for j in $(seq 1 "$N_CHR"); do
        nohup python vcf2geno.py \
            "./03.group.af.vcf/group.af.${pop}.chr${j}.vcf.gz" \
            "./04.group.af.txt/group.af.${pop}.chr${j}.txt" \
            > "nohup_vcf2geno_${pop}_chr${j}.out" 2>&1 &
    done
done < temp_pop.txt

wait
echo "[Step 6] Done."


# ============================================================
# Step 7: Merge per-population txt into per-chromosome geno files
# ============================================================
# paste joins files column-wise (space-delimited).
# The glob group.af.*.chr{N}.txt must expand in a consistent
# population order across all chromosomes — verify with ls
# before running if you are unsure.
#
# NOTE: >> appends; if chr{N}.geno.txt already exists from a
# previous run this will duplicate data. Delete old files first.

echo "[Step 7] Merging per-population txt files by chromosome..."

cd 04.group.af.txt || exit 1

for i in $(seq 1 "$N_CHR"); do
    paste -d ' ' group.af.*.chr${i}.txt \
        >> "../05.geno/chr${i}.geno.txt"
done

cd ..
echo "[Step 7] Done."


# ============================================================
# Step 8: Count columns in a geno file (diagnostic check)
# ============================================================
# Reports the number of columns (= number of populations × 2)
# in a representative geno file. Adjust the filename as needed.

echo "[Step 8] Column count in chr6.geno.txt:"
awk '{print NF; exit}' 05.geno/chr6.geno.txt


# ============================================================
# Step 9: Rename files — strip underscores from population names
# ============================================================
# If a population name contains underscores (e.g. Nera_Siciliana),
# bcftools glob patterns may not match reliably. This renames
# affected files to remove the underscore.
#
# Extend this block for any other population names with underscores.

echo "[Step 9] Renaming files with underscores in population names..."

for file in 04.group.af.txt/group.af.Nera_Siciliana.chr*; do
    [[ -e "$file" ]] || continue
    newfile="${file/Nera_Siciliana/NeraSiciliana}"
    mv -v "$file" "$newfile"
done

echo "[Step 9] Done."


# ============================================================
# Step 10: Extract unique population list from file names
# ============================================================
# Derives the canonical population order directly from the
# file names produced in Step 6, and saves it to population_list.txt.
# This file is used as the key in Step 11.

echo "[Step 10] Extracting population list from file names..."

ls 04.group.af.txt/group.af.*.chr*.txt \
    | sed -E 's|^04\.group\.af\.txt/group\.af\.(.*)\.chr[0-9]+\.txt$|\1|' \
    | sort -u \
    > population_list.txt

echo "[Step 10] Done. Populations:"
cat population_list.txt


# ============================================================
# Step 11: Build eco covariate files and split by variable
# ============================================================
# Stage A: join population_list.txt with env.txt to produce eco.txt,
#   ordered to match the population order in the geno files.
#   eco.txt has populations as rows and eco1–eco18 as columns.
#
# Stage B: strip the header row and ecotype column, transpose so
#   each row is one eco variable across all populations, then
#   split into 18 separate files: eco1.txt … eco18.txt.
#
# NOTE: awk inline comments (#) inside the awk program are not
#   standard POSIX awk — they work in gawk but may fail in mawk.
#   Run with gawk if your system default is mawk.

echo "[Step 11] Building eco covariate files..."

# Stage A: build ordered eco.txt
awk '
    BEGIN { FS="\t"; OFS=" " }
    NR==FNR { ecotype_order[$1]=NR; next }
    $2 in ecotype_order && !seen[$2]++ {
        eco_data[ecotype_order[$2]] = $2 OFS \
            $(NF-17) OFS $(NF-16) OFS $(NF-15) OFS $(NF-14) OFS \
            $(NF-13) OFS $(NF-12) OFS $(NF-11) OFS $(NF-10) OFS \
            $(NF-9)  OFS $(NF-8)  OFS $(NF-7)  OFS $(NF-6)  OFS \
            $(NF-5)  OFS $(NF-4)  OFS $(NF-3)  OFS $(NF-2)  OFS \
            $(NF-1)  OFS $NF
    }
    END {
        print "ecotype", \
              "eco1",  "eco2",  "eco3",  "eco4",  "eco5",  "eco6", \
              "eco7",  "eco8",  "eco9",  "eco10", "eco11", "eco12", \
              "eco13", "eco14", "eco15", "eco16", "eco17", "eco18"
        for (i=1; i<=length(ecotype_order); i++)
            if (eco_data[i] != "") print eco_data[i]
    }
' population_list.txt env.txt > 06.eco/eco.txt

# Stage B: strip header + ecotype column, transpose, split into 18 files
awk 'NR>1 { for(i=2; i<=NF; i++) printf "%s ", $i; print "" }' \
    06.eco/eco.txt > 06.eco/temp.txt

awk '{
    for (i=1; i<=NF; i++) a[i,NR]=$i
}
END {
    for (i=1; i<=NF; i++) {
        for (j=1; j<=NR; j++) printf "%s ", a[i,j]
        print ""
    }
}' 06.eco/temp.txt > 06.eco/transposed_eco.txt && rm 06.eco/temp.txt

awk '{ print > "06.eco/eco" NR ".txt" }' 06.eco/transposed_eco.txt

echo "[Step 11] Done."


# ============================================================
# Step 12: Verify BayPass binary
# ============================================================
# Download: https://forge.inrae.fr/mathieu.gautier/baypass_public
# Extract:  tar -xzvf baypass_public-master.tar.gz
# Compile:  follow README in sources/

echo "[Step 12] Testing BayPass binary..."

if [[ ! -x "$BAYPASS" ]]; then
    echo "ERROR: g_baypass not found or not executable at: $BAYPASS"
    exit 1
fi

"$BAYPASS" \
    -gfile  "05.geno/chr1.geno.txt" \
    -efile  "06.eco/eco1.txt"       \
    -outprefix "07.covresult/eco1.chr.1.anacovis" \
    -nthreads "$NTHREADS_COVIS"

echo "[Step 12] Test run complete."


# ============================================================
# Step 13A: Run BayPass covariate analysis (AnaCovis mode)
# ============================================================
# Processes eco variables in pairs to limit peak thread usage.
# Both jobs for a given chromosome must finish before the next
# chromosome starts (enforced by wait).
#
# Peak concurrency: 2 jobs × NTHREADS_COVIS threads.

echo "[Step 13A] Running BayPass AnaCovis analysis..."

for i in $(seq 1 2 "$N_CHR"); do
    for j in $(seq 1 "$N_CHR"); do

        echo "[$(date '+%H:%M:%S')] Starting eco${i} chr${j}..."
        nohup "$BAYPASS" \
            -gfile     "05.geno/chr${j}.geno.txt"            \
            -efile     "06.eco/eco${i}.txt"                  \
            -outprefix "07.covresult/eco${i}.chr.${j}.anacovis" \
            -nthreads  "$NTHREADS_COVIS"                     \
            > "nohup_eco${i}_chr${j}.out" 2>&1 &

        next=$(( i + 1 ))
        if [[ "$next" -le "$N_CHR" ]]; then
            echo "[$(date '+%H:%M:%S')] Starting eco${next} chr${j}..."
            nohup "$BAYPASS" \
                -gfile     "05.geno/chr${j}.geno.txt"               \
                -efile     "06.eco/eco${next}.txt"                   \
                -outprefix "07.covresult/eco${next}.chr.${j}.anacovis" \
                -nthreads  "$NTHREADS_COVIS"                         \
                > "nohup_eco${next}_chr${j}.out" 2>&1 &
        fi

        wait
        echo "[$(date '+%H:%M:%S')] eco${i}${next:+/eco${next}} chr${j} complete."

    done
done

echo "[Step 13A] AnaCovis runs complete."


# ============================================================
# Step 13B: Run BayPass final analysis (with omega matrix)
# ============================================================
# Uses the omega covariance matrix estimated in Step 13A as a
# prior, which corrects for population structure.
# Runs sequentially (no background jobs) to avoid overloading
# the system at NTHREADS_FINAL threads per job.
#
# NOTE: omegafile must exist from Step 13A before this runs.

echo "[Step 13B] Running BayPass final analysis..."

for i in $(seq 1 "$N_CHR"); do
    for j in $(seq 1 "$N_CHR"); do
        echo "[$(date '+%H:%M:%S')] Final run eco${i} chr${j}..."
        "$BAYPASS" \
            -gfile     "05.geno/chr${j}.geno.txt"                          \
            -efile     "06.eco/eco${i}.txt"                                \
            -omegafile "07.covresult/eco${i}.chr.${j}.anacovis_mat_omega.out" \
            -outprefix "09.final.result/eco${i}.chr.${j}.final"            \
            -nthreads  "$NTHREADS_FINAL"
    done
done

echo "[Step 13B] Final analysis complete."


# ============================================================
# Step 14: Initialise merged result files from chr1
# ============================================================
# Creates eco{i}_merged_file.out for each eco variable by
# copying the chr1 result (including the header line).

echo "[Step 14] Initialising merged result files from chr1..."

for i in $(seq 1 "$N_CHR"); do
    cp "09.final.result/eco${i}.chr.1.final_summary_betai_reg.out" \
       "09.final.result/eco${i}_merged_file.out"
done

echo "[Step 14] Done."


# ============================================================
# Step 15: Append chr2–18 results to merged files
# ============================================================
# sed '1d' strips the header from chromosomes 2–18 before
# appending, so the merged file has exactly one header row.

echo "[Step 15] Appending chr2–${N_CHR} results to merged files..."

for i in $(seq 1 "$N_CHR"); do
    for j in $(seq 2 "$N_CHR"); do
        file="09.final.result/eco${i}.chr.${j}.final_summary_betai_reg.out"
        if [[ -f "$file" ]]; then
            sed '1d' "$file" >> "09.final.result/eco${i}_merged_file.out"
        else
            echo "WARNING: $file not found — skipping."
        fi
    done
done

echo "[Step 15] Done."


# ============================================================
# Step 16: Annotate merged results with SNP positions
# ============================================================
# Extracts CHR, POS, and SNP ID from the BIM file, pastes them
# alongside the merged BayPass output, then selects the columns
# needed for downstream visualisation: CHR, POS, SNP, BF(dB).
#
# BIM column order: CHR(1) SNPID(2) CM(3) BP(4) A1(5) A2(6)
# awk '{print $1,$4,$2}' → CHR, BP, SNPID
#
# After paste, eco{i}_merged_file.final.txt has columns:
#   CHR(1) POS(2) SNPID(3) | BayPass cols...
# Column 10 in the BayPass output = BF(dB); verify against
# your BayPass version before running.

echo "[Step 16] Annotating results with SNP positions..."

awk '{print $1, $4, $2}' chr919qcldbeagle.bim > 919.positions.txt

for i in $(seq 1 "$N_CHR"); do
    paste \
        919.positions.txt \
        "09.final.result/eco${i}_merged_file.out" \
        > "09.final.result/eco${i}_merged_file.final.txt"

    awk '{print $1, $2, $3, $10}' \
        "09.final.result/eco${i}_merged_file.final.txt" \
        > "09.final.result/eco${i}.txt"
done

echo "[Step 16] Done."


echo "========================================================"
echo "BayPass pipeline complete: $(date)"
echo "========================================================"