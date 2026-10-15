================================================================================
SCRIPT: run_stage3.sh
================================================================================
Description: Runs stage 3 steps 1-17 for one focal (reference) species in one working folder: input preparation
(exon list and header-cleaned reference CDS from the BED), step 1 (mafExtractor, or a copy of existing step 1
files), steps 2-17, and with ONE_ROW=1 the step 14 dropped-row filter. Stops at the first failed step.

Usage:
  bash run_stage3.sh <workdir> <focal_species> <bed> <cds_fasta> <maf> <tree_newick>
  Run detached and low-priority, e.g.
  setsid nohup nice -n 19 ionice -c3 env THREADS=8 bash run_stage3.sh ... > <workdir>.log 2>&1 < /dev/null &

Inputs:
  - <workdir>: A new folder (created if needed); it must not already contain step outputs unless START_STEP is used.
  - <focal_species>: The species prefix of the MAF sequence names (e.g. Arabidopsis_thaliana, Homo_sapiens).
  - <bed>: Tab-separated: chrom, start (0-based), end, name, score, strand, transcript ID. A header line and '#'
    lines are skipped; a trailing '+' on chromosome names is removed.
  - <cds_fasta>: Reference CDS nucleotide sequences named by the BED transcript IDs (a '.p' suffix is accepted).
    The run stops if a BED transcript has no CDS.
  - <maf>: The reference-anchored MAF. Give '-' instead when STEP1_SRC is set
    or START_STEP > 1; PROVENANCE.txt then records "MAF not given (not used: ...)". A path, when given, must exist
    and is recorded with its size and date (not checksummed).
  - <tree_newick>: Newick tree with the MAF species names as tips (step 8).

Environment (all optional):
  - ENV_BIN: bin folder of a Python env with Biopython and pandas (default /opt/miniconda3/envs/busco_env/bin).
  - MAFEXTRACTOR_BIN: folder holding mafExtractor (default /opt/cactus/cactus/bin).
  - GNU parallel must already be on PATH (step 4).
  - THREADS (8); JOBS (parallel step 1 chunks, default THREADS; each chunk reads the whole MAF).
  - STEP3=dominant|strict (default dominant: step3_dominant_chain.py; strict: step3_filtering_species.py); JUMP (50000).
  - ONE_ROW=1|0 (default 1: step 8 --one-row-per-species plus the step 14 dropped-row filter).
  - STEP1_SRC: folder of existing step 1 files (<species>.<chrom>_<start>_<end>.txt) to copy instead of running
    mafExtractor; mafExtractor is then not needed.
  - ALLOW_MISSING=1: continue when some exons have no step 1 file (default: stop). Their transcripts fail step 9.
  - START_STEP=1-17 (default 1): resume in an existing workdir (see below).
  - FLOOR_GB (50): refuse to start with less free space in the workdir.

Checks before anything runs (the run stops, without writing STATUS, if one fails):
  inputs exist; mafExtractor on PATH if step 1 runs without STEP1_SRC; GNU parallel if step 4 runs; Biopython and
  pandas; the tree can be read (warning if the focal species is not a tip); free space; no outputs of steps
  >= START_STEP already in the workdir; outputs of steps < START_STEP present.

Step 1 checks (after extraction or copying): every exon must have a non-empty file containing reference lines
(problems listed in step1_problems.tsv, run stops); exons without a file are listed in step1_missing.txt and stop
the run unless ALLOW_MISSING=1. Extraction failures are listed in step1/step1_failed.tsv.

Resuming (START_STEP=N): input preparation reruns (cheap and deterministic); steps before N are skipped; outputs of
step N and later must be moved out of the workdir first (the script never deletes anything).
Outputs per step: 1 step1/; 2 step2/; 3 step3/, step3_summary.csv; 4 step4/; 6 step4/extraction_summary_fixed.tsv;
8 step8/, step8_one_row_dropped/; 9 step9/; 11 step8/master_summary_final_validation.tsv; 12 step12/; 13 step13/;
14 step14/ (and step14_raw/, step14_removed_dropped_row_indels.tsv with ONE_ROW=1); 15 step15/;
16 step16_starts_stops.tsv; 17 step17_all_species.tsv. Steps 5, 7 and 10 only report (in the log).

Status and records in the workdir:
  - STATUS: "RUNNING since <date>", then "SUCCEEDED <date>" or "FAILED at <step> (exit N) <date>"
    (read by tools/run_monitor.py).
  - PROVENANCE.txt: appended on every run: inputs with md5s, options, repo commit, tool paths.
  - coords.tsv, N_EXONS, ref_raw.fna, ref.fna, <focal_species>.bed: prepared inputs.

What to check after a run: STATUS; step1_missing.txt and step1_problems.tsv empty; refused files in
step9/step9_trimming_report.tsv; step14_removed_dropped_row_indels.tsv (ONE_ROW=1); audits/audit_stage3_endtoend.py
and, with STEP3=dominant, audits/audit_optionA.py.
