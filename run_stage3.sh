#!/usr/bin/env bash
# run_stage3.sh
#
# LINGUA stage 3 (alignment extraction), steps 1-17, for one focal species of a Progressive Cactus alignment.
# Corrected pipeline (not the Brassicaceae paper reproduction):
#   step 3   dominant collinear chain per species and exon (default), keeping the orthologous blocks and dropping
#            blocks on other chromosomes/strands, more than JUMP bp apart, out of order or overlapping;
#            STEP3=strict uses the original filter, which removes the species from the whole exon instead
#   step 4   no genome lookup (LINGUA 911c653)
#   step 8   Python port; with ONE_ROW=1 (default) each species keeps only its row with the most aligned bases
#   step 9   per-exon padding removal (LINGUA 69ca427), which fixes the exon-junction errors
#
# Usage (on Dayhoff):
#   bash run_stage3.sh <workdir> <focal_species> <bed> <cds_fasta> <maf> <tree_newick>
#   <focal_species> must match the species prefix of the MAF sequence names (e.g. Homo_sapiens).
#   <bed>: tab-separated, header line, columns chrom, start (0-based), end, name, score, strand, transcript ID.
#   <cds_fasta>: CDS nucleotide sequences named by the BED transcript IDs.
# Environment: ENV_BIN (python env bin, default busco_env), THREADS (8), JOBS (parallel step 1 chunks, default THREADS), STEP3 (dominant|strict), JUMP (50000),
#   ONE_ROW (1|0), STEP1_SRC (folder of existing step 1 files to copy instead of running mafExtractor),
#   FLOOR_GB (50, refuse to start below this free space)
# Run detached and low-priority, e.g.
#   setsid nohup nice -n 19 ionice -c3 bash run_stage3.sh ... > <workdir>.log 2>&1 < /dev/null &
set -euo pipefail

WD="${1:?usage: run_stage3.sh <workdir> <focal_species> <bed> <cds_fasta> <maf> <tree_newick>}"
SP="${2:?focal species}"
BED_SRC=$(readlink -f "${3:?bed}")
CDS_SRC=$(readlink -f "${4:?cds fasta}")
MAF=$(readlink -f "${5:?maf}")
export TREE_FILE=$(readlink -f "${6:?tree}")
REPO=$(cd "$(dirname "$0")" && pwd)
export PATH=${ENV_BIN:-/opt/miniconda3/envs/busco_env/bin}:/opt/cactus/cactus/bin:$PATH   # Biopython, pandas, mafExtractor
S=$REPO/scripts
T=${THREADS:-8}
JOBS=${JOBS:-$T}
STEP3=${STEP3:-dominant}
JUMP=${JUMP:-50000}
ONE_ROW=${ONE_ROW:-1}
FLOOR_GB=${FLOOR_GB:-50}
B=$SP.bed

for f in "$BED_SRC" "$CDS_SRC" "$MAF" "$TREE_FILE"; do [ -e "$f" ] || { echo "missing input $f"; exit 2; }; done
[ -z "${STEP1_SRC:-}" ] || [ -d "$STEP1_SRC" ] || { echo "missing STEP1_SRC $STEP1_SRC"; exit 2; }
case $STEP3 in dominant|strict) ;; *) echo "STEP3 must be dominant or strict"; exit 2;; esac
mkdir -p "$WD" && cd "$WD"
free=$(df -BG --output=avail . | tail -1 | tr -dc 0-9)
[ "$free" -ge "$FLOOR_GB" ] || { echo "only $free GB free (< $FLOOR_GB GB)"; exit 2; }

step() { echo "== $1 $(date -Is)"; shift; "$@"; }
{
    echo "run_stage3.sh $(date -Is)"
    echo "repo: $REPO ($(git -C "$REPO" rev-parse --short HEAD 2>/dev/null || echo 'no git'))"
    echo "focal species: $SP; STEP3=$STEP3 JUMP=$JUMP ONE_ROW=$ONE_ROW THREADS=$T JOBS=$JOBS"
    for f in "$BED_SRC" "$CDS_SRC" "$TREE_FILE"; do echo "input: $f md5 $(md5sum < "$f" | cut -c1-32)"; done
    echo "input: $MAF ($(stat -c %s "$MAF") bytes, $(stat -c %y "$MAF"))"
    [ -n "${STEP1_SRC:-}" ] && echo "step 1 copied from: $STEP1_SRC"
} > PROVENANCE.txt

tr -d '\r' < "$BED_SRC" | awk 'BEGIN{FS=OFS="\t"} !/^#/ && NR>1 {sub(/\+$/,"",$1)} {print}' > "$B"
step prepare python3 "$S/prepare_species_inputs.py" "$B" "$CDS_SRC" coords.tsv N_EXONS ref_raw.fna

mkdir -p step1
if [ -n "${STEP1_SRC:-}" ]; then
    tail -n +2 coords.tsv | while IFS=$'\t' read -r c s e; do
        f="$SP.${c}_${s}_${e}.txt"
        if [ -f "$STEP1_SRC/$f" ]; then cp -p "$STEP1_SRC/$f" step1/; else echo "$f"; fi
    done > step1_missing.txt
else
    mkdir -p chunks
    tail -n +2 coords.tsv | split -n l/"$JOBS" -d -a 3 - chunks/c_
    pids=()
    for c in chunks/c_???; do
        { printf 'seq\tstart\tstop\n'; cat "$c"; } > "$c.tsv"
        bash "$S/step1_mafextractor.sh" "$SP" "$c.tsv" "$MAF" step1 > "$c.log" 2>&1 &
        pids+=($!)
    done
    for p in "${pids[@]}"; do wait "$p" || { echo "a step 1 chunk failed (chunks/*.log)"; exit 1; }; done
    : > step1_missing.txt
fi
n=$(find step1 -name '*.txt' | wc -l)
echo "step 1: $n files for $(cat N_EXONS) exons; $(wc -l < step1_missing.txt) missing (step1_missing.txt)"
[ "$n" -gt 0 ] || { echo "no step 1 files"; exit 1; }

step clean_headers python3 "$S/clean_nucleotide_fasta_headers.py" -f ref_raw.fna -b "$B" -o ref.fna
step step2 bash "$S/step2_missingbp_summaries.sh" step1 step2
if [ "$STEP3" = dominant ]; then
    step step3 python3 "$S/step3_dominant_chain.py" -i step1 -o step3 -r "$SP" -s step3_summary.csv -j "$JUMP" -t "$T"
else
    step step3 python3 "$S/step3_filtering_species.py" -i step1 -o step3 -r "$SP" -s step3_summary.csv -t "$T"
fi
step step4 bash "$S/step4_reverse_complement_add_missingbp.sh" -m step3 -o step4 -b "$B" -t "$T"
step step5 python3 "$S/step5_validate_reverse_complement.py" step4 step1 "$B" "$SP" "$T"
step step6 python3 "$S/step6_update_summary_files.py" step4/extraction_summary.tsv "$B" "$SP"
step step7 python3 "$S/step7_first_stage_validation.py" step1 step3 step4 "$SP"
one=(); [ "$ONE_ROW" = 1 ] && one=(--one-row-per-species)
step step8 python3 "$S/step8_stitching_exons.py" step4 "$B" step8 "$SP" "${one[@]}"
step step9 python3 "$S/step9_trimming_duplicate_bases.py" step8 step9 ref.fna "$B" "$T"
step step10 python3 "$S/step10_second_stage_validation.py" step9 ref.fna "$SP"
step step11 python3 "$S/step11_update_stitched_exons_summary.py" step8/master_summary.tsv step9 "$B" ref.fna "$SP"
V=step8/master_summary_final_validation.tsv
step step12 bash "$S/step12_identify_start_codons.sh" step9 step12 "$B" "$SP" "$V"
step step13 bash "$S/step13_identify_stop_codons.sh" step9 step13 "$B" "$SP" "$V"
step step14 bash "$S/step14_identify_indels.sh" step9 step14 "$B" "$SP" "$V"
step step15 bash "$S/step15_identify_first_and_last_three_bases.sh" step9 step15 "$SP" "$V"
step step16 python3 "$S/step16_extract_start_and_stops.py" step9 step16_starts_stops.tsv
step step17 python3 "$S/step17_combine_all_species_info.py" step9 step3 "$B" "$SP" step17_all_species.tsv
echo "== done $(date -Is)"
