"""
audit_optionA.py

Runs ON Dayhoff. Independent checks of the Option A scripts (step3_dominant_chain.py, step8b_one_row_per_species.py)
on their real outputs, written separately from the scripts themselves:
  step 3  A. every step-3 s-line exists verbatim in the step-1 file (nothing invented or altered)
          B. reference lines: step 3 keeps every step-1 reference block
          C. per (exon, species): one chromosome and strand, consecutive MAF-order gaps <= jump
          D. the kept (chromosome, strand) group has the most aligned bases, and the kept chain is the
             largest chain of that group (recomputed independently)
          E. no species with aligned bases in step 1 is lost from an exon
  step 8b F. no alignment has a species in more than one row after step 8b
          G. for every logged merge, the kept row has >= bases than each dropped row
          H. step 8 master_summary 'aligned_species' vs rows actually in the step 9 file (stale counts)
  step 9  I. step 9 report statuses
Usage: python3 audit_optionA.py <species_workdir> <reference_species> [jump]
"""
import collections
import csv
import os
import sys

wd, ref = sys.argv[1:3]
JUMP = int(sys.argv[3]) if len(sys.argv) > 3 else 50000


def slines(path):
    return [l.rstrip("\n") for l in open(path) if l.startswith("s")]


def blocks(path):
    out, cur = [], None
    for l in open(path):
        if l.startswith("a"):
            if cur is not None:
                out.append(cur)
            cur = []
        elif l.startswith("s") and cur is not None:
            p = l.split()
            sp, ch = p[1].split(".", 1) if "." in p[1] else (p[1], "unknown")
            cur.append((sp, ch, int(p[2]), int(p[3]), p[4], sum(c != "-" for c in p[6])))
    if cur is not None:
        out.append(cur)
    return out


bad = collections.Counter()
ex = collections.defaultdict(list)
n_exons = 0
for f in sorted(os.listdir(f"{wd}/step1")):
    n_exons += 1
    s1, s3 = slines(f"{wd}/step1/{f}"), slines(f"{wd}/step3/{f}")
    if set(s3) - set(s1):
        bad["A step3 line not in step1"] += 1; ex["A"].append(f)
    r1 = [l for l in s1 if l.split()[1].split(".")[0] == ref]
    if set(r1) - set(s3):
        bad["B reference block lost"] += 1; ex["B"].append(f)
    b1, b3 = blocks(f"{wd}/step1/{f}"), blocks(f"{wd}/step3/{f}")
    occ1, occ3 = collections.defaultdict(list), collections.defaultdict(list)
    for b in b1:
        for o in b:
            if o[0] != ref:
                occ1[o[0]].append(o)
    for b in b3:
        for o in b:
            if o[0] != ref:
                occ3[o[0]].append(o)
    for sp, lst in occ1.items():
        kept = occ3.get(sp, [])
        if sum(o[5] for o in lst) > 0 and not kept:
            bad["E species lost from exon"] += 1; ex["E"].append((f, sp))
            continue
        if not kept:
            continue
        if len({(o[1], o[4]) for o in kept}) > 1:
            bad["C >1 chromosome/strand kept"] += 1; ex["C"].append((f, sp))
        last = None
        for o in kept:
            if last is not None and abs(o[2] - last) > JUMP:
                bad["C kept jump > threshold"] += 1; ex["C2"].append((f, sp)); break
            last = o[2] + o[3]
        # collinearity of the kept blocks (MAF order): species start never before the previous end, gap <= JUMP
        for a, b in zip(kept, kept[1:]):
            if b[2] < a[2] + a[3]:
                bad["C kept blocks out of order or overlapping"] += 1; ex["Co"].append((f, sp)); break
        # optimality, recomputed independently: memoized recursion over blocks sorted by species start,
        # a successor must also come later in MAF order
        best_total = 0
        groups = collections.defaultdict(list)
        for k, o in enumerate(lst):
            groups[(o[1], o[4])].append((k, o))
        for g in groups.values():
            g = sorted(g, key=lambda x: (x[1][2], x[0]))
            memo = {}

            def best_from(i):
                if i in memo:
                    return memo[i]
                ki, oi = g[i]
                end = oi[2] + oi[3]
                tail = 0
                for j in range(i + 1, len(g)):
                    kj, oj = g[j]
                    if kj > ki and 0 <= oj[2] - end <= JUMP:
                        tail = max(tail, best_from(j))
                memo[i] = oi[5] + tail
                return memo[i]
            sys.setrecursionlimit(10000)
            best_total = max([best_total] + [best_from(i) for i in range(len(g))])
        if sum(o[5] for o in kept) != best_total:
            bad["D kept chain not optimal"] += 1; ex["D"].append((f, sp, sum(o[5] for o in kept), best_total))

dup = 0
for f in os.listdir(f"{wd}/step9"):
    if f.endswith(".aln"):
        c = collections.Counter(l[1:].split(".")[0] for l in open(f"{wd}/step9/{f}") if l.startswith(">"))
        dup += any(v > 1 for v in c.values())
if dup:
    bad["F species in >1 row after step 8b"] += dup
g_bad = g_n = 0
# the drop log moved into step 8 at commit ac4ce53 (step8/one_row_per_species_log.tsv); older runs
# wrote step8b_log.tsv. Say which one was read, and warn instead of skipping silently.
g_log = next((p for p in (f"{wd}/step8/one_row_per_species_log.tsv", f"{wd}/step8b_log.tsv") if os.path.exists(p)), None)
if g_log:
    for r in csv.DictReader(open(g_log), delimiter="\t"):
        g_n += 1
        dropped = [int(x) for x in str(r["dropped_bases"]).replace(",", ";").split(";") if x.strip()]
        if any(int(r["kept_bases"]) < d for d in dropped):
            g_bad += 1
else:
    print("WARNING: no step 8 drop log found (step8/one_row_per_species_log.tsv or step8b_log.tsv); check G not run")
if g_bad:
    bad["G kept row smaller than dropped"] += g_bad

stale = total = 0
ms = f"{wd}/step8b/master_summary.tsv" if os.path.exists(f"{wd}/step8b/master_summary.tsv") else f"{wd}/step8/master_summary.tsv"
if os.path.exists(ms):
    for r in csv.DictReader(open(ms), delimiter="\t"):
        p = f"{wd}/step9/{r['file_name']}.aln"
        if os.path.exists(p):
            total += 1
            n = sum(1 for l in open(p) if l.startswith(">"))
            stale += int(r["aligned_species"]) != n
st9 = collections.Counter()
for l in open(f"{wd}/step9/step9_trimming_report.tsv"):
    if l.startswith(("#", "file")):
        continue
    st9[l.split("\t")[5]] += 1

print(f"== {os.path.basename(wd)}: {n_exons} exons")
print("violations:", dict(bad) if bad else "none")
for k, v in ex.items():
    print(f"  example {k}: {v[:3]}")
print(f"H master_summary aligned_species differs from step-9 rows: {stale} of {total} alignments")
print(f"I step 9 statuses: {dict(st9)}")
print(f"step 8 drop log: {g_log or none}; entries {g_n}; kept row smaller than a dropped row: {g_bad}")
