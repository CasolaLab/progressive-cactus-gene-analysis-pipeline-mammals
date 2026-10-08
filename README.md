# 3_alignment_extraction_mammals_repo

LINGUA stage 3 (alignment extraction) for the T2T mammal Progressive Cactus alignment, focal
species *Homo sapiens*. Its output feeds stage 4, which is adapted locally in
`~/Dropbox/de-novo_genes/Ne-DNGs/zoonomia-PC-datasets/95_species/t2t_mammals/lingua_t2t_mammals/`.

Assembled 2026-10-08 from the corrected Brassicaceae stage 3 scripts, following
`t2t_mammals/HANDOFF_stage3_from_brassicaceae_2026-10-08.md`. `PROVENANCE.tsv` gives the source,
version and md5 of every file. `README_LINGUA_stage3.md` and `docs/` are the LINGUA stage 3
documentation, unchanged.

**Not the paper-reproduction repo.** `/mnt/D/projects/Brassicales/softmasked/LINGUA_stage3_paper/`
reproduces the Brassicaceae paper outputs with their bugs; do not use it for mammals.

## Inputs on Dayhoff

- HAL: `/mnt/ExData/t2t_mammals.hal` (14 species, tree identical to `trimmed_t2t_tree.nwk`).
- Human-referenced MAF: `/mnt/ExData/t2t_mammals_human.maf` (2026-08-04).
- Still needed: the human CDS BED (must be named `Homo_sapiens.bed`, LF line endings), the human
  CDS nucleotide FASTA for step 9, and the Newick tree for step 8 (`TREE_FILE`).

## Versions chosen, and why

| Step | Script | Choice |
|---|---|---|
| 3 | `step3_filtering_species.py` | repo strict filter (default). `step3_dominant_chain.py` is the audited alternative; choose before running |
| 4 | `step4_reverse_complement_add_missingbp.sh` | LINGUA 911c653: no genome lookup, so no genome FASTAs or samtools needed |
| 8 | `step8_stitching_exons.py` | Python port, run with `--one-row-per-species` (split rows come from step 8 keying rows by species.chromosome) |
| 9 | `step9_trimming_duplicate_bases.py` | LINGUA **69ca427**, not main: main removes the real base when the padding base equals its neighbor |
| 12-14 | repo 911c653 | reference length read from the step 11 `actual_nt` column; step 14 `cds_length` is now each species' own length |
| others | repo main 911c653 | unchanged |

## Checks before trusting output

1. Step 9 junction errors are invisible in the reference row; check non-reference rows at a locus
   whose padding base equals its neighbor.
2. Run `audits/audit_stage3_endtoend.py` (every base at its expected column, positive control);
   `audits/audit_optionA.py` if the dominant-chain step 3 is used.
3. Steps 12-14 have not yet run on RefSeq names (`NC_060925.1`, `NP_…`); confirm Chromosome is
   not NA in their outputs.
4. Steps 15 and 17 write rows in nondeterministic order; compare content, not bytes.

## Running on Dayhoff

- `export PATH=/opt/miniconda3/envs/busco_env/bin:/opt/cactus/cactus/bin:$PATH`
- Step 1 (mafExtractor) costs ~31 s per exon, scaled by MAF size; the human MAF is ~88 GB.
- Run with `nice -n 19 ionice -c3`, few threads, under `tools/run_monitor.py` (thermal and disk guard).
