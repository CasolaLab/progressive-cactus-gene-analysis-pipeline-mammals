================================================================================
SCRIPT: step14_identify_indels.sh
================================================================================
Description: For every non-reference species in each step 9 alignment, finds runs of columns where the species has a gap and the reference a base (deletion) or the species a base and the reference a gap (insertion), and reports a run if its length is not a multiple of 3 (frameshift), or if it is a 3-bp insertion whose bases are a stop codon. Transcript ID, chromosome and strand come from the step 11 summary; reference_length is the step 11 actual_nt column; cds_length is that species' own length with gaps removed. Writes one <file>_indels.tsv per alignment and master_indels.tsv; alignment files are read in sorted order.

Usage:
  ./step14_identify_indels.sh <step9_trimmed_dir> <output_dir> <bed_file> <reference_species> <step11_summary_tsv>

Inputs:
  - <step9_trimmed_dir>: Trimmed alignments from step 9.
  - <output_dir>: Destination for the indel tables.
  - <bed_file>: Reference BED file (still a required argument).
  - <reference_species>: The reference species name prefix.
  - <step11_summary_tsv>: The step 11 *_final_validation.tsv (must contain actual_nt).

Notes:
  - With step 8 --one-row-per-species, exons of a species' dropped rows become gaps and would be reported as deletions; run_stage3.sh writes step 14 to step14_raw/ and runs filter_dropped_row_indels.py to remove them (see its README).
