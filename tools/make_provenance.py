#!/usr/bin/env python3
"""
make_provenance.py

Rebuild PROVENANCE.tsv: one row per script, audit and tool in the repo, with its source, the md5 of the
version used for the 2026-10-08 T2T human run (Dayhoff repo commit 6416c38), its current md5 and whether it
is unchanged, changed or new since then. Source, version and note are carried over from the previous
PROVENANCE.tsv; files with no earlier row get the note given with --new-note.

Usage: python3 tools/make_provenance.py OLD_PROVENANCE.tsv DAYHOFF_6416c38_MD5.tsv OUT.tsv [--changed-note TEXT]
       [--new-note TEXT]
       DAYHOFF_6416c38_MD5.tsv: two columns, md5 and path, for the files of commit 6416c38.
"""
import argparse
import csv
import glob
import hashlib
import os
import sys

PATTERNS = ["scripts/*", "audits/*", "tools/*.py", "run_stage3.sh"]


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("old")
    ap.add_argument("dayhoff_md5")
    ap.add_argument("out")
    ap.add_argument("--changed-note", default="modified after 6416c38")
    ap.add_argument("--new-note", default="added after 6416c38")
    a = ap.parse_args()
    for p in (a.old, a.dayhoff_md5):
        if not os.path.exists(p):
            sys.exit(f"ERROR: not found: {p}")
    old = {r["file"]: r for r in csv.DictReader(open(a.old), delimiter="\t")}
    ref = {}
    for line in open(a.dayhoff_md5):
        m, f = line.rstrip("\n").split("\t", 1)
        ref[f] = m
    files = sorted({f for p in PATTERNS for f in glob.glob(p) if os.path.isfile(f) and "__pycache__" not in f})
    rows = []
    for f in files:
        now = md5(f)
        before = ref.get(f, "")
        status = "new" if not before else ("unchanged" if before == now else "changed")
        o = old.get(f, {})
        note = o.get("note", "")
        if status == "changed":
            note = f"{a.changed_note}; earlier: {note}" if note else a.changed_note
        elif status == "new":
            note = a.new_note
        rows.append([f, o.get("source", "new"), o.get("version", ""), before, now, status, note])
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["file", "source", "version", "md5_dayhoff_6416c38", "md5_current", "status", "note"])
        w.writerows(rows)
    counts = {s: sum(r[5] == s for r in rows) for s in ("unchanged", "changed", "new")}
    print(f"{len(rows)} files: {counts}", file=sys.stderr)


if __name__ == "__main__":
    main()
