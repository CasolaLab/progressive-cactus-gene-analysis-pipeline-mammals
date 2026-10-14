#!/usr/bin/env python3
"""
filter_dropped_row_indels.py  (run after step 14 when step 8 used --one-row-per-species)

With --one-row-per-species, step 8 keeps each species' row with the most aligned bases; the exons of its other
rows become gaps in the kept row, and step 14 reports those gaps as deletions although they are exons aligned on
another chromosome, not deletions. This removes those indels.

For every gene with dropped rows (<step8_dir>_one_row_dropped/<gene>.aln, written by step 8), the reference
positions covered by each dropped row are found by applying step 9's column trimming to the step 8 columns:
kept columns = reference base columns minus the extra-base columns step 9 lists in step9_trimming_report.tsv.
The trimmed step 8 reference row must equal step 9's reference row, or the script stops. A step 14 row is
removed when its species is the dropped row's species and its reference interval overlaps one of these positions;
removed rows are written to <removed_tsv> with whether they also extend past the dropped positions.

Usage: python3 filter_dropped_row_indels.py <step8_dir> <step9_dir> <step14_raw_dir> <step14_out_dir>
                                            <removed_tsv> <reference_species>
"""
import csv
import os
import sys


def read_aln(path):
    recs = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                recs.append([line[1:].split()[0], []])
            elif recs:
                recs[-1][1].append(line.strip())
    return [(h, "".join(s)) for h, s in recs]


def ranges(positions):
    out = []
    for p in sorted(positions):
        if out and p == out[-1][1] + 1:
            out[-1][1] = p
        else:
            out.append([p, p])
    return [tuple(r) for r in out]


def main():
    if len(sys.argv) != 7:
        sys.exit(__doc__)
    s8, s9, raw, out, removed_tsv, ref = sys.argv[1:7]
    dropped_dir = s8.rstrip("/") + "_one_row_dropped"
    report = os.path.join(s9, "step9_trimming_report.tsv")
    for p in (s8, s9, raw, report):
        if not os.path.exists(p):
            sys.exit(f"Error: not found: {p}")
    if os.path.exists(out):
        sys.exit(f"Error: {out} exists; refusing to overwrite")
    with open(report) as f:
        rep = {r["file"]: r for r in csv.DictReader((l for l in f if not l.startswith("#")), delimiter="\t")}

    masked = {}
    genes = sorted(f for f in os.listdir(dropped_dir) if f.endswith(".aln")) if os.path.isdir(dropped_dir) else []
    not_in_step9 = 0
    for f in genes:
        r = rep.get(f)
        if not r or r["status"] != "OK" or not os.path.exists(os.path.join(s9, f)):
            not_in_step9 += 1
            continue
        main_rows = read_aln(os.path.join(s8, f))
        ref_rows = [s for h, s in main_rows if h.split(".")[0] == ref]
        if len(ref_rows) != 1:
            sys.exit(f"Error: {f}: expected one {ref} row in step 8, found {len(ref_rows)}")
        ref_s8 = ref_rows[0]
        extra = {int(c) - 1 for c in r["extra_base_columns_removed_1based"].split(",") if c not in ("", "NA")}
        keep = [c for c, ch in enumerate(ref_s8) if ch != "-" and c not in extra]
        ref_s9 = [s for h, s in read_aln(os.path.join(s9, f)) if h.split(".")[0] == ref][0]
        if "".join(ref_s8[c] for c in keep) != ref_s9:
            sys.exit(f"Error: {f}: step 9 column mapping does not reproduce its reference row; stopping")
        pos_of = {c: k + 1 for k, c in enumerate(keep)}
        for h, s in read_aln(os.path.join(dropped_dir, f)):
            if len(s) != len(ref_s8):
                sys.exit(f"Error: {f}: dropped row {h} has {len(s)} columns, step 8 alignment {len(ref_s8)}")
            pos = {pos_of[c] for c, ch in enumerate(s) if ch != "-" and c in pos_of}
            if pos:
                masked.setdefault((f[:-4], h.split(".")[0]), []).extend(ranges(pos))

    os.makedirs(out)
    files = sorted(f for f in os.listdir(raw) if f.endswith(".tsv"))
    n_in = n_out = n_rm = n_past = 0
    with open(removed_tsv, "w") as rem:
        wrote = False
        for fn in files:
            with open(os.path.join(raw, fn)) as fh:
                lines = fh.read().splitlines()
            if not lines:
                open(os.path.join(out, fn), "w").close()
                continue
            hdr = lines[0].split("\t")
            fi, si = hdr.index("file_name"), hdr.index("Species")
            ai, bi = hdr.index("IndelStart_RefPos"), hdr.index("IndelEnd_RefPos")
            if not wrote:
                rem.write("source_file\textends_past_dropped_positions\t" + lines[0] + "\n")
                wrote = True
            kept = [lines[0]]
            for line in lines[1:]:
                p = line.split("\t")
                a, b = int(p[ai]), int(p[bi])
                hits = [(x, y) for x, y in masked.get((p[fi], p[si]), []) if a <= y and b >= x]
                if hits:
                    inside = any(x <= a and b <= y for x, y in hits)
                    rem.write(f"{fn}\t{'no' if inside else 'yes'}\t{line}\n")
                    n_rm += 1
                    n_past += not inside
                else:
                    kept.append(line)
            n_in += len(lines) - 1
            n_out += len(kept) - 1
            with open(os.path.join(out, fn), "w") as o:
                o.write("\n".join(kept) + "\n")
    assert n_in == n_out + n_rm
    print(f"filter_dropped_row_indels: {len(genes)} genes with dropped rows ({not_in_step9} not written by step 9); "
          f"{len(masked)} gene/species with dropped positions; step 14 rows in {n_in}, kept {n_out}, "
          f"removed {n_rm} ({n_past} also extend past the dropped positions) -> {out}")


if __name__ == "__main__":
    main()
