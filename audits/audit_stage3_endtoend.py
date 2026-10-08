"""
audit_stage3_endtoend.py

Runs ON Dayhoff. Base-by-base check that the final stage 3 alignments (step 9) contain exactly the step 3 blocks
they should, at the right positions, for one reference species. Independent of the pipeline code.

For every final alignment, every non-reference row (species.chromosome) and every exon of the gene:
  expected = bases of that species' step-3 blocks on the row's chromosome, at reference positions inside the CDS
             (complemented on minus-strand genes, since rows are in transcript orientation)
  CHECK 1  every expected base is present at its column                         (missing -> error)
  CHECK 2  every non-gap base in the row is an expected base                     (unexplained base -> error)
  CHECK 3  blocks of that species on OTHER chromosomes do not appear in this row  (counted inside CHECK 2)
  CHECK 4  the reference row equals the CDS from ref.fna
  CHECK 5  no species is in more than one row (Option A) / report how many (strict)
Optionally, with a second run folder (the strict run), rows whose species and exons were untouched by the new
step 3 and by step 8b must be byte-identical between the two runs (CHECK 6, regression).

Usage: python3 audit_stage3_endtoend.py <species_workdir> <reference_species> <bed_prefix_species> [strict_workdir]
"""
import collections
import csv
import os
import sys

wd, ref, sp_bed = sys.argv[1:4]
strict = sys.argv[4] if len(sys.argv) > 4 else None
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

tx = collections.defaultdict(list)
for line in open(f"{wd}/{sp_bed}.bed"):
    if line.startswith("#") or not line.strip():
        continue
    f = line.rstrip("\n").split("\t")
    tx[f[6]].append((int(f[1]), int(f[2]), f[5], f"{sp_bed}.{f[0]}_{f[1]}_{f[2]}.txt"))


def fasta(path):
    d, k = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            k = line[1:].split()[0]
            d[k] = []
        elif k:
            d[k].append(line)
    return {k: "".join(v).upper() for k, v in d.items()}


cds = fasta(f"{wd}/ref.fna")


def aln(path):
    recs, k = [], None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            recs.append([line[1:], ""])
        elif recs:
            recs[-1][1] += line
    return recs


def blocks(path):
    """species -> list of (chrom, [(ref_pos, base)])"""
    out = collections.defaultdict(list)
    ref_line = None
    rows = []
    for line in open(path):
        if line.startswith("a"):
            if ref_line:
                for sp, ch, seq in rows:
                    cols, pos = [], ref_line[0]
                    for r, q in zip(ref_line[1], seq):
                        if r != "-":
                            if q != "-":
                                cols.append((pos, q.upper()))
                            pos += 1
                    out[sp].append((ch, cols))
            ref_line, rows = None, []
        elif line.startswith("s"):
            p = line.split()
            sp, ch = p[1].split(".", 1)
            if sp == ref:
                ref_line = (int(p[2]), p[6])
            else:
                rows.append((sp, ch, p[6]))
    if ref_line:
        for sp, ch, seq in rows:
            cols, pos = [], ref_line[0]
            for r, q in zip(ref_line[1], seq):
                if r != "-":
                    if q != "-":
                        cols.append((pos, q.upper()))
                    pos += 1
            out[sp].append((ch, cols))
    return out


err = collections.Counter()
ok = collections.Counter()
examples = collections.defaultdict(list)
step9 = f"{wd}/step9"
B = {}
for f in sorted(os.listdir(step9)):
    if not f.endswith(".aln"):
        continue
    t = next((t for t in tx if f"_{t}_" in f), None)
    if t is None:
        err["no transcript match"] += 1
        continue
    exons = sorted(tx[t], reverse=(tx[t][0][2] == "-"))
    strand = exons[0][2]
    idx, k = {}, 0
    for s, e, _, _ in exons:
        for p in (range(s, e) if strand == "+" else range(e - 1, s - 1, -1)):
            idx[p] = k
            k += 1
    recs = aln(f"{step9}/{f}")
    refrow = recs[0][1]
    if refrow.replace("-", "").upper() != cds.get(t, ""):
        err["CHECK4 reference row != CDS"] += 1
        examples["4"].append(f)
    else:
        ok["CHECK4"] += 1
    colmap = [i for i, c in enumerate(refrow) if c != "-"]
    exp_by_sp = collections.defaultdict(dict)          # species -> chrom -> {col: base}
    for s, e, _, fn in exons:
        if fn not in B:
            B[fn] = blocks(f"{wd}/step3/{fn}")
        for sp, bl in B[fn].items():
            for ch, cols in bl:
                d = exp_by_sp[sp].setdefault(ch, {})
                for pos, q in cols:
                    if pos in idx:
                        d[colmap[idx[pos]]] = q.translate(COMP) if strand == "-" else q
    seen = collections.Counter(h.split(".")[0] for h, _ in recs[1:])
    err["CHECK5 alignments with a species in >1 row"] += any(v > 1 for v in seen.values())
    for h, row in recs[1:]:
        sp, ch = h.split(".", 1)
        exp = exp_by_sp.get(sp, {}).get(ch, {})
        other = {c for c2, d in exp_by_sp.get(sp, {}).items() if c2 != ch for c in d}
        miss = sum(1 for c, q in exp.items() if row[c].upper() != q and not (row[c] == "-" and False))
        miss_gap = sum(1 for c, q in exp.items() if row[c] == "-")
        wrong = sum(1 for c, q in exp.items() if row[c] not in "-" and row[c].upper() != q)
        extra = [c for c, x in enumerate(row) if x != "-" and c not in exp]
        extra_other = sum(1 for c in extra if c in other)
        ok["rows"] += 1
        ok["expected bases"] += len(exp)
        if miss_gap:
            err["CHECK1 expected base missing (gap)"] += miss_gap
            examples["1"].append((f, h))
        if wrong:
            err["CHECK1 expected base wrong"] += wrong
            examples["1w"].append((f, h))
        if extra:
            err["CHECK2 unexplained base in row"] += len(extra)
            if extra_other:
                err["CHECK3 base from another-chromosome block"] += extra_other
            examples["2"].append((f, h, len(extra)))

reg = collections.Counter()
if strict:
    trimmed = set()
    p3 = f"{wd}/step3_summary.csv"
    for r in csv.DictReader(open(p3)):
        if r.get("status") == "trimmed":
            trimmed.add((r["file_name"], r["species"]))
    merged = set()
    if os.path.exists(f"{wd}/step8b_log.tsv"):
        for r in csv.DictReader(open(f"{wd}/step8b_log.tsv"), delimiter="\t"):
            merged.add((r["file"], r["species"]))
    for f in sorted(os.listdir(step9)):
        if not f.endswith(".aln") or not os.path.exists(f"{strict}/step9/{f}"):
            continue
        t = next(t for t in tx if f"_{t}_" in f)
        a = dict(aln(f"{step9}/{f}"))
        b = dict(aln(f"{strict}/step9/{f}"))
        for h, row in a.items():
            sp = h.split(".")[0]
            touched = any((fn, sp) in trimmed for _, _, _, fn in tx[t]) or (f, sp) in merged
            if touched or h not in b:
                continue
            reg["untouched rows compared"] += 1
            if row != b[h]:
                reg["untouched rows DIFFERENT from strict"] += 1
                examples["6"].append((f, h))

print(f"== end-to-end audit of {wd}")
print("passed:", dict(ok))
print("errors:", {k: v for k, v in err.items() if v} or "none")
if strict:
    print("regression vs strict:", dict(reg))
for k, v in examples.items():
    print(f"  example {k}: {v[:4]}")
