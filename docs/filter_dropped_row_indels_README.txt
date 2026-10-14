================================================================================
SCRIPT: filter_dropped_row_indels.py  (run after step 14 when step 8 used --one-row-per-species)
================================================================================
Description: Removes the step 14 indels that exist only because step 8 kept one row per species: the exons of a species' dropped rows (aligned on another chromosome) become gaps in its kept row and look like deletions. For every gene in <step8_dir>_one_row_dropped/, the reference positions covered by each dropped row are found by applying step 9's trimming to the step 8 columns (reference base columns minus the extra-base columns listed in step9_trimming_report.tsv); the script stops if this does not reproduce step 9's reference row. A step 14 row is removed when its species is the dropped row's species and its interval overlaps those positions.

Usage:
  python3 filter_dropped_row_indels.py <step8_dir> <step9_dir> <step14_raw_dir> <step14_out_dir> <removed_tsv> <reference_species>

Inputs:
  - <step8_dir>: Step 8 output (its sibling <step8_dir>_one_row_dropped/ holds the dropped rows).
  - <step9_dir>: Step 9 output, including step9_trimming_report.tsv.
  - <step14_raw_dir>: Step 14 output to filter (per-file tables and master_indels.tsv).
  - <step14_out_dir>: New folder for the filtered tables (must not exist).
  - <removed_tsv>: Every removed row, with source file and whether it also extends past the dropped positions (a real deletion next to a dropped exon is removed with it).
  - <reference_species>: The reference species name.
