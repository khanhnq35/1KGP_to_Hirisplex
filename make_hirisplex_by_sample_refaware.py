#!/usr/bin/env python3
"""
make_hirisplex_by_sample_refaware.py

Purpose
-------
Convert a 1KGP genotype matrix (TSV) into a HiRisPlex-S by-sample CSV with
allele-aware counting (REF/ALT-aware + strand complements + simple 1-base INDELs).

Inputs
------
1) Matrix (TSV, tab-delimited):
   ~/Documents/2025.1/Bio_Labs/1kgp_out/rs41_GT_matrix_with_refalt.tsv

   Header example:
     CHROM  POS  ID  REF  ALT  HG00096  HG00097  ...
   Rows example:
     5  33951693  rs16891982  C  A  1|1  1|0  ...

   - Each sample column is raw %GT from VCF (e.g., 0|1, 1/1, 0/0, ./., 0, 1).

2) Panel (CSV, comma-delimited, UTF-8):
   ~/Documents/2025.1/Bio_Labs/hirisplexs.csv

   Required columns (case-insensitive names accepted):
     - SNP  (or rsID or ID): e.g., rs312262906
     - Allele: e.g., A

   Panel order defines output column order. If duplicate rsIDs occur, the
   first occurrence is kept.

Output
------
CSV UTF-8:
  ~/Documents/2025.1/Bio_Labs/1kgp_out/hirisplex_by_sample_refaware.csv

Columns:
  sampleid, rs312262906_A, rs11547464_A, ..., rs8051733_C
Values per cell: 0, 1, 2, or NA

Counting Logic (allele-aware)
-----------------------------
Let target allele be 'A':
- If target matches ALT (or complement(ALT)): count number of '1's in GT.
- If target matches REF (or complement(REF)): count number of '0's in GT.
- Simple 1-base INDELs:
    * If ALT == REF + allele  -> treat as ALT target (inserted base)
    * If REF == ALT + allele  -> treat as REF target (deleted base)
- If nothing matches and GT is present      -> 0
- If GT is missing (./. or .)               -> NA
- Haploid genotypes (single allele like '0' or '1'):
    * When counting ALT: 0->0, 1->1, .->NA
    * When counting REF: 0->1, 1->0, .->NA

Notes
-----
- No pandas/numpy used. Only csv/os/sys.
- Robust to whitespace and mixed separators ('|' or '/').
- Logs informative messages with [+], ⚠️, ✅.

Run
---
python3 make_hirisplex_by_sample_refaware.py
"""

import os
import sys
import csv

BASE = os.path.expanduser("~/Documents/2025.1/Bio_Labs")
OUTDIR = os.path.join(BASE, "1kgp_out")
MATRIX_TSV = os.path.join(OUTDIR, "rs41_GT_matrix_with_refalt.tsv")
PANEL_CSV = os.path.join(BASE, "hirisplexs.csv")
OUT_CSV = os.path.join(OUTDIR, "hirisplex_by_sample_refaware.csv")


def die(msg, code=1):
    sys.stderr.write("❌ " + msg + "\n")
    sys.exit(code)


def exists_nonempty(path):
    return os.path.isfile(path) and os.path.getsize(path) > 0


def read_panel(path):
    """Read panel CSV -> ordered list of (rsid, allele), keeping first dup only."""
    if not exists_nonempty(path):
        die("Thiếu file panel: %s" % path)

    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        # Sniff delimiter if needed (but panel is specified as CSV)
        reader = csv.reader(f)
        rows = list(reader)

    if not rows:
        die("Panel rỗng: %s" % path)

    header = [c.strip() for c in rows[0]]
    # Find SNP column (accept SNP/rsID/ID) and Allele column
    def find_col(cands):
        for cand in cands:
            for i, h in enumerate(header):
                if h.lower() == cand.lower():
                    return i
        return None

    snp_idx = find_col(["SNP", "rsID", "ID"])
    allele_idx = find_col(["Allele", "ALLELE", "allele"])
    if snp_idx is None or allele_idx is None:
        die("Không tìm thấy cột 'SNP/rsID/ID' hoặc 'Allele' trong panel: %s" % path)

    seen = set()
    ordered = []
    for row in rows[1:]:
        if not row or max(snp_idx, allele_idx) >= len(row):
            continue
        rsid = row[snp_idx].strip()
        allele = row[allele_idx].strip()
        if not rsid or not allele:
            continue
        # Normalize rsid to lower 'rs' prefix + digits if possible
        rsid = rsid.strip()
        allele = allele.upper()
        if rsid not in seen:
            ordered.append((rsid, allele))
            seen.add(rsid)
    return ordered


def dna_complement(base):
    """Return complement of single base A<->T, C<->G. Otherwise return None."""
    base = base.upper()
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return comp.get(base)


def gt_count_for_alt(gt):
    """Count number of ALT alleles from GT string."""
    if gt == "." or gt == "./.":
        return "NA"
    # tolerate '|' or '/' or no separator
    if "|" in gt:
        a, b = gt.split("|", 1)
    elif "/" in gt:
        a, b = gt.split("/", 1)
    else:
        # haploid
        if gt == ".":
            return "NA"
        return "NA" if gt not in ("0", "1") else ("1" if gt == "1" else "0")
    if a == "." or b == ".":
        return "NA"
    cnt = 0
    if a == "1":
        cnt += 1
    if b == "1":
        cnt += 1
    return str(cnt)


def gt_count_for_ref(gt):
    """Count number of target REF alleles (i.e., count zeros)."""
    if gt == "." or gt == "./.":
        return "NA"
    if "|" in gt:
        a, b = gt.split("|", 1)
    elif "/" in gt:
        a, b = gt.split("/", 1)
    else:
        # haploid: '0' means REF present (count 1), '1' means ALT (count 0)
        if gt == ".":
            return "NA"
        return "NA" if gt not in ("0", "1") else ("1" if gt == "0" else "0")
    if a == "." or b == ".":
        return "NA"
    cnt = 0
    if a == "0":
        cnt += 1
    if b == "0":
        cnt += 1
    return str(cnt)


def decide_target_side(ref, alt, allele):
    """
    Decide whether target allele corresponds to REF or ALT.
    Returns "REF", "ALT", or "NONE".

    Also handle complement equivalence and simple 1-base INDEL patterns:
      - ALT == REF + allele  -> ALT
      - REF == ALT + allele  -> REF
    Only applied when len(allele)==1.
    """
    r = ref.upper()
    a = alt.upper()
    t = allele.upper()

    # Direct match (SNP)
    if t == a:
        return "ALT"
    if t == r:
        return "REF"

    # Complement match
    tc = dna_complement(t)  # not used directly
    rc = dna_complement(r) if len(r) == 1 else None
    ac = dna_complement(a) if len(a) == 1 else None
    if rc is not None and t == rc:
        return "REF"
    if ac is not None and t == ac:
        return "ALT"

    # Simple 1-base INDEL heuristics
    if len(t) == 1:
        # insertion: ALT = REF + t
        if a == (r + t):
            return "ALT"
        # deletion: REF = ALT + t
        if r == (a + t):
            return "REF"

    return "NONE"


def read_matrix(path):
    """Read matrix TSV into:
       - sample_names (list)
       - data dict: rsid -> {"REF": ref, "ALT": alt, "GT": [gt per sample index]}
    """
    if not exists_nonempty(path):
        die("Thiếu file matrix: %s" % path)

    with open(path, "r", encoding="utf-8") as f:
        header_line = f.readline()
        if not header_line:
            die("Matrix rỗng: %s" % path)
        header = [h.strip() for h in header_line.rstrip("\n").split("\t")]
        # Expect at least 6 columns
        if len(header) < 6:
            die("Header matrix không hợp lệ.")
        # positions: CHROM POS ID REF ALT samples...
        sample_names = header[5:]

        data = {}
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 5:
                continue
            chrom, pos, rsid, ref, alt = fields[:5]
            gts = fields[5:]
            # Align length; if missing fill with "."
            if len(gts) < len(sample_names):
                gts += ["./."] * (len(sample_names) - len(gts))
            rsid = rsid.strip()
            data[rsid] = {"REF": ref.strip(), "ALT": alt.strip(), "GT": gts}
    return sample_names, data


def main():
    if not os.path.isdir(OUTDIR):
        die("Thiếu thư mục output: %s" % OUTDIR)

    print("[+] Đọc panel:", PANEL_CSV)
    panel = read_panel(PANEL_CSV)
    if not panel:
        die("Panel không có bản ghi hợp lệ.")

    print("[+] Đọc matrix:", MATRIX_TSV)
    samples, matrix = read_matrix(MATRIX_TSV)
    print("[+] Số sample:", len(samples))
    print("[+] Số rsID trong matrix:", len(matrix))
    print("[+] Số rsID trong panel:", len(panel))

    # Prepare output columns in order
    out_cols = ["sampleid"] + [f"{rs}_{allele}" for (rs, allele) in panel]

    # Build per-sample rows
    # Initialize with 'NA' for all target columns
    rows = []
    for s in samples:
        rows.append([s] + ["NA"] * len(panel))

    # For each target in panel, compute counts
    all_na_flags = [True] * len(panel)  # to warn columns that are entirely NA

    for idx, (rsid, allele) in enumerate(panel):
        col_idx = idx + 1  # because first column is sampleid
        info = matrix.get(rsid)
        if info is None:
            # Not found -> all NA; warn later
            continue
        ref = info["REF"]
        alt = info["ALT"]
        side = decide_target_side(ref, alt, allele)

        # Fill counts
        gts = info["GT"]
        if side == "ALT":
            for si, gt in enumerate(gts):
                val = gt_count_for_alt(gt)
                rows[si][col_idx] = val
                if val != "NA":
                    all_na_flags[idx] = False
        elif side == "REF":
            for si, gt in enumerate(gts):
                val = gt_count_for_ref(gt)
                rows[si][col_idx] = val
                if val != "NA":
                    all_na_flags[idx] = False
        else:
            # Unmatched allele -> if GT present then 0 else NA
            for si, gt in enumerate(gts):
                if gt in (".", "./."):
                    val = "NA"
                else:
                    # Normalize haploid/diploid presence check
                    parts = []
                    if "|" in gt:
                        parts = gt.split("|", 1)
                    elif "/" in gt:
                        parts = gt.split("/", 1)
                    else:
                        parts = [gt]
                    if any(p == "." for p in parts):
                        val = "NA"
                    else:
                        val = "0"
                rows[si][col_idx] = val
                if val != "NA":
                    all_na_flags[idx] = False

    # Write CSV
    os.makedirs(OUTDIR, exist_ok=True)
    with open(OUT_CSV, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(out_cols)
        for r in rows:
            writer.writerow(r)

    print("✅ Đã tạo:", OUT_CSV)
    print("Tổng: {} mẫu × {} SNP_Allele".format(len(samples), len(panel)))

    # Warnings for all-NA columns
    warned = 0
    for idx, flag in enumerate(all_na_flags):
        if flag:
            rsid, allele = panel[idx]
            print("⚠️  {}_{} toàn NA (không khớp hoặc không có trong matrix)".format(rsid, allele))
            warned += 1
    if warned == 0:
        print("[+] Không có cột nào toàn NA.")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        die("Bị huỷ bởi người dùng.", 130)