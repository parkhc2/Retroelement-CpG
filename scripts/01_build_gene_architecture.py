#!/usr/bin/env python
"""
01_build_gene_architecture.py
=============================
Classify every GENCODE gene by promoter CpG-island status and the
upstream/downstream retroelement (Alu, L1, LTR) context of its TSS.

Inputs (UCSC hg38 / GENCODE, see 00_download_raw.sh)
  --gtf        gencode.v26.annotation.gtf.gz   (GTEx v8 annotation)
  --cpg        cpgIslandExt.txt.gz             (UCSC CpG island track)
  --rmsk       rmsk.txt.gz                     (UCSC RepeatMasker track)
  --gtf-tags   gencode.v49.annotation.gtf.gz   (only for --tss mane; supplies
               MANE_Select / Ensembl_canonical transcript tags, which do not
               exist in v26)

TSS definition (--tss)
  longest  : longest annotated transcript per gene (main analysis)
  mane     : MANE_Select > Ensembl_canonical (from --gtf-tags) > longest
             (alternative TSS definition)

Distance conventions
  * All distances are strand-oriented gap distances in bp between the TSS
    and the nearest element of each class on the upstream / downstream side.
    An element overlapping the TSS has distance 0.  NA = no element of that
    class on that side of the chromosome.
  * up_class / dn_class = the closest of Alu, L1, LTR on that side if it lies
    within 10 kb of the TSS, otherwise "none".
  * config: Alu-Alu (Alu closest on both sides), Alu-mixed (one side),
    LTR/L1-hybrid (neither), incomplete_flank (a side has no element <=10 kb).
  * group12: CpG-island genes are further split by Alu distance
    (nearest Alu on either side, i.e. min(Alu_up, Alu_dn)); thresholds were
    fixed from the empirical Alu-TSS distance distribution before any
    expression analysis (Fig. S1 of the manuscript):
        Alu-Alu   : <=2 kb | 2-4 kb | up<=7 kb or dn<=8 kb | far
        Alu-mixed : <=2 kb | up<=5 kb or dn<=4 kb          | far
        LTR/L1    : <=2 kb (min of up/dn nearest element)  | far
    Non-CpG-island genes are not split by distance.

Output: TSV with one row per gene (ENSG without version).
"""
import argparse
import gzip
import re
import sys

import numpy as np
import pandas as pd

RE_ATTR = re.compile(r'(\S+) "([^"]*)"')


def parse_gtf_transcripts(path, want_tags=False):
    """Return DataFrame of transcripts: gene_id, ENSG, ENST, chrom, strand, start, end, gene_type, tags."""
    rows = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] != "transcript":
                continue
            attrs = {}
            tags = []
            for k, v in RE_ATTR.findall(f[8]):
                if k == "tag":
                    tags.append(v)
                else:
                    attrs[k] = v
            rows.append((attrs["gene_id"], attrs["transcript_id"], f[0], f[6], int(f[3]), int(f[4]),
                         attrs.get("gene_type", ""), attrs.get("gene_name", ""),
                         ";".join(tags) if want_tags else ""))
    tx = pd.DataFrame(rows, columns=["gene_id", "transcript_id", "chrom", "strand", "start", "end",
                                     "gene_type", "gene_name", "tags"])
    # drop the chrY pseudo-autosomal copies (gene_id suffix _PAR_Y) so that each ENSG is unique
    tx = tx[~tx.gene_id.str.endswith("_PAR_Y")].copy()
    tx["ENSG"] = tx.gene_id.str.split(".").str[0]
    tx["ENST"] = tx.transcript_id.str.split(".").str[0]
    tx["length"] = tx.end - tx.start + 1
    tx["tss"] = np.where(tx.strand == "+", tx.start, tx.end)
    return tx


def pick_longest(tx):
    """Longest transcript per gene; ties broken by transcript ID (deterministic)."""
    s = tx.sort_values(["ENSG", "length", "ENST"], ascending=[True, False, True])
    out = s.drop_duplicates("ENSG").copy()
    out["tss_source"] = "longest"
    return out


def pick_mane(tx26, tx_tags):
    """MANE_Select > Ensembl_canonical transcript (tags from a newer GENCODE) mapped onto v26 transcripts;
    genes without a tagged transcript present in v26 fall back to the longest v26 transcript."""
    tagged = tx_tags[["ENSG", "ENST", "tags"]].copy()
    tagged["is_mane"] = tagged.tags.str.contains("MANE_Select")
    tagged["is_canon"] = tagged.tags.str.contains("Ensembl_canonical")
    tagged = tagged[tagged.is_mane | tagged.is_canon]
    # priority: MANE_Select first
    tagged["prio"] = np.where(tagged.is_mane, 0, 1)
    tagged = tagged.sort_values(["ENSG", "prio"]).drop_duplicates("ENSG")
    m = tx26.merge(tagged[["ENSG", "ENST", "prio"]], on=["ENSG", "ENST"], how="inner")
    m["tss_source"] = np.where(m.prio == 0, "MANE_Select", "Ensembl_canonical")
    longest = pick_longest(tx26)
    rest = longest[~longest.ENSG.isin(m.ENSG)].copy()
    rest["tss_source"] = "longest_fallback"
    cols = ["gene_id", "ENSG", "ENST", "chrom", "strand", "tss", "gene_type", "gene_name", "tss_source"]
    return pd.concat([m[cols], rest[cols]], ignore_index=True)


def load_cpg(path):
    cpg = pd.read_csv(path, sep="\t", header=None, usecols=[1, 2, 3], names=["chrom", "start", "end"])
    return cpg


def load_rmsk(path):
    """Return dict class -> DataFrame(chrom,start,end) for Alu (SINE/Alu), L1 (LINE/L1), LTR (class LTR)."""
    rk = pd.read_csv(path, sep="\t", header=None, usecols=[5, 6, 7, 11, 12],
                     names=["chrom", "start", "end", "repClass", "repFamily"])
    return {
        "Alu": rk[rk.repFamily == "Alu"][["chrom", "start", "end"]],
        "L1": rk[rk.repFamily == "L1"][["chrom", "start", "end"]],
        "LTR": rk[rk.repClass == "LTR"][["chrom", "start", "end"]],
    }


def build_index(df):
    """Per-chromosome sorted starts, ends, and running max of ends (for overlap detection)."""
    idx = {}
    for ch, sub in df.groupby("chrom"):
        s = sub.sort_values("start")
        st = s.start.to_numpy()
        en = s.end.to_numpy()
        idx[ch] = (st, en, np.maximum.accumulate(en))
    return idx


def dist_left_right(idx, chrom, tss):
    """Gap distance from each TSS to the nearest element on the genomic left and right.
    0 if an element overlaps the TSS (start <= tss < end, UCSC 0-based half-open vs 1-based TSS
    as used in the original analysis); NaN if no element on that side."""
    L = np.full(len(tss), np.nan)
    R = np.full(len(tss), np.nan)
    for ch in np.unique(chrom):
        m = chrom == ch
        if ch not in idx:
            continue
        st, en, cummax_en = idx[ch]
        t = tss[m]
        i = np.searchsorted(st, t, side="right")          # number of elements starting at/before t
        left = np.full(len(t), np.nan)
        has = i > 0
        ce = cummax_en[i[has] - 1]
        left[has] = np.where(ce > t[has], 0, t[has] - ce)
        right = np.full(len(t), np.nan)
        hr = i < len(st)
        right[hr] = st[i[hr]] - t[hr]
        right = np.where(left == 0, 0, right)
        L[m] = left
        R[m] = right
    return L, R


def classify(genes, cpg, retro, window=10000):
    """genes: DataFrame with ENSG, chrom, strand, tss. Returns architecture table."""
    out = genes.copy()
    chrom = out.chrom.to_numpy()
    tss = out.tss.to_numpy().astype(np.int64)
    plus = out.strand.to_numpy() == "+"

    L, _ = dist_left_right(build_index(cpg), chrom, tss)
    out["cpg_island"] = L == 0

    classes = list(retro)
    for k in classes:
        L, R = dist_left_right(build_index(retro[k]), chrom, tss)
        out[f"{k}_up"] = np.where(plus, L, R)
        out[f"{k}_dn"] = np.where(plus, R, L)

    for side in ("up", "dn"):
        M = out[[f"{k}_{side}" for k in classes]].to_numpy()
        M2 = np.where(np.isnan(M), np.inf, M)
        j = M2.argmin(axis=1)
        d = M2.min(axis=1)
        out[f"{side}_class"] = np.where(d <= window, np.array(classes)[j], "none")
        out[f"{side}_dist"] = np.where(np.isinf(d), np.nan, d)

    n_alu = (out.up_class == "Alu").astype(int) + (out.dn_class == "Alu").astype(int)
    incomplete = (out.up_class == "none") | (out.dn_class == "none")
    out["config"] = np.select([incomplete, n_alu == 2, n_alu == 1],
                              ["incomplete_flank", "Alu-Alu", "Alu-mixed"], "LTR/L1-hybrid")

    au = out.Alu_up.fillna(np.inf)
    ad = out.Alu_dn.fillna(np.inf)
    alu_min = np.minimum(au, ad)
    any_min = np.minimum(out.up_dist.fillna(np.inf), out.dn_dist.fillna(np.inf))
    cg = out.cpg_island.to_numpy()
    cf = out.config
    out["group12"] = np.select(
        [
            ~cg & (cf == "incomplete_flank"), ~cg,
            cf == "incomplete_flank",
            (cf == "Alu-Alu") & (alu_min <= 2000),
            (cf == "Alu-Alu") & (alu_min <= 4000),
            (cf == "Alu-Alu") & ((au <= 7000) | (ad <= 8000)),
            cf == "Alu-Alu",
            (cf == "Alu-mixed") & (alu_min <= 2000),
            (cf == "Alu-mixed") & ((au <= 5000) | (ad <= 4000)),
            cf == "Alu-mixed",
            any_min <= 2000,
        ],
        [
            "nonCGI_incomplete", "nonCGI_" + cf.astype(str),
            "CGI_incomplete",
            "CGI_AluAlu_le2k", "CGI_AluAlu_2-4k", "CGI_AluAlu_4-7/8k", "CGI_AluAlu_far",
            "CGI_AluMixed_le2k", "CGI_AluMixed_2-5/4k", "CGI_AluMixed_far",
            "CGI_hybrid_le2k",
        ],
        "CGI_hybrid_far",
    )
    # continuous / symmetric-bin helpers used by robustness analyses
    out["Alu_min_dist"] = np.where(np.isinf(alu_min), np.nan, alu_min)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--cpg", required=True)
    ap.add_argument("--rmsk", required=True)
    ap.add_argument("--tss", choices=["longest", "mane"], default="longest")
    ap.add_argument("--gtf-tags", help="newer GENCODE GTF carrying MANE_Select/Ensembl_canonical tags (for --tss mane)")
    ap.add_argument("--chroms", default="chr1-22,chrX,chrY", help="'chr1-22,chrX,chrY' (default) or 'all'")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    print("[1/4] parsing GTF", file=sys.stderr)
    tx = parse_gtf_transcripts(a.gtf)
    if a.tss == "longest":
        genes = pick_longest(tx)
    else:
        if not a.gtf_tags:
            sys.exit("--tss mane requires --gtf-tags")
        tx_tags = parse_gtf_transcripts(a.gtf_tags, want_tags=True)
        genes = pick_mane(tx, tx_tags)
    if a.chroms != "all":
        keep = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
        genes = genes[genes.chrom.isin(keep)]
    genes = genes.sort_values(["chrom", "tss"]).reset_index(drop=True)

    print("[2/4] loading CpG islands and RepeatMasker", file=sys.stderr)
    cpg = load_cpg(a.cpg)
    retro = load_rmsk(a.rmsk)

    print("[3/4] classifying %d genes" % len(genes), file=sys.stderr)
    res = classify(genes, cpg, retro)

    cols = ["ENSG", "gene_id", "gene_name", "chrom", "strand", "tss", "tss_source", "ENST", "gene_type",
            "cpg_island", "Alu_up", "Alu_dn", "L1_up", "L1_dn", "LTR_up", "LTR_dn",
            "up_class", "up_dist", "dn_class", "dn_dist", "Alu_min_dist", "config", "group12"]
    res = res[cols].rename(columns={"gene_type": "biotype"})
    print("[4/4] writing", a.out, file=sys.stderr)
    res.to_csv(a.out, sep="\t", index=False, float_format="%.0f")
    print(res.group12.value_counts().to_string(), file=sys.stderr)


if __name__ == "__main__":
    main()
