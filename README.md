# Retroelement-CpG

Genome-wide classification of human genes by **promoter CpG-island status** and the
**retroelement context of the transcription start site (TSS)** — the identity and distance
of the nearest Alu, L1 and LTR elements upstream and downstream of each TSS.

This table is the gene-annotation layer used in

> Park H-C, Choi H, Kim S, Choi J, Oh J-H, Rhyu M-G, Hong S-J.
> *Promoter-proximal Alu–CpG island architecture stratifies transcriptional stability during hematopoietic differentiation.*

## Contents

```
scripts/00_download_raw.sh              fetch GENCODE and UCSC inputs into data/raw/
scripts/01_build_gene_architecture.py   build the classification table
results/gene_architecture_hg38_v26_longestTSS.tsv.gz   TSS = longest transcript (main)
results/gene_architecture_hg38_v26_maneTSS.tsv.gz      TSS = MANE Select / Ensembl canonical
```

## Inputs

| File | Source | Content |
|---|---|---|
| `gencode.v26.annotation.gtf.gz` | GENCODE release 26 (GRCh38, Ensembl 88) | gene and transcript models |
| `gencode.v49.annotation.gtf.gz` | GENCODE release 49 | `MANE_Select` / `Ensembl_canonical` transcript tags (only for `--tss mane`) |
| `cpgIslandExt.txt.gz` | UCSC hg38 | CpG islands |
| `rmsk.txt.gz` | UCSC hg38 RepeatMasker | Alu (`repFamily = Alu`), L1 (`repFamily = L1`), LTR (`repClass = LTR`) |

```bash
pip install numpy pandas            # Python >= 3.10
bash scripts/00_download_raw.sh     # ~290 MB

python scripts/01_build_gene_architecture.py \
    --gtf data/raw/gencode.v26.annotation.gtf.gz \
    --cpg data/raw/cpgIslandExt.txt.gz --rmsk data/raw/rmsk.txt.gz \
    --tss longest --out results/gene_architecture_hg38_v26_longestTSS.tsv

# alternative TSS definition
python scripts/01_build_gene_architecture.py \
    --gtf data/raw/gencode.v26.annotation.gtf.gz --gtf-tags data/raw/gencode.v49.annotation.gtf.gz \
    --cpg data/raw/cpgIslandExt.txt.gz --rmsk data/raw/rmsk.txt.gz \
    --tss mane --out results/gene_architecture_hg38_v26_maneTSS.tsv
```

Run time ≈ 20 s (longest) / 40 s (mane); memory < 2 GB.

## Processing steps

| Step | Rule |
|---|---|
| Gene universe | All GENCODE v26 genes on chr1–22, X, Y. chrM genes and the chrY pseudo-autosomal copies (`gene_id` suffix `_PAR_Y`) are removed → **58,137 genes**, one row per Ensembl gene ID (version stripped). |
| TSS | `--tss longest`: 5′ end of the longest annotated transcript of the gene (ties broken by transcript ID). `--tss mane`: MANE Select transcript, else Ensembl canonical transcript (tags from `--gtf-tags` mapped onto v26 transcript IDs), else the longest transcript (`tss_source` records which rule applied). |
| CpG-island gene | The TSS lies inside a `cpgIslandExt` interval (`cpg_island = True`). |
| Retroelement distances | For each class (Alu, L1, LTR) and each side of the TSS — upstream / downstream, oriented by gene strand — the gap in bp between the TSS and the nearest element. 0 if the element overlaps the TSS; empty if no element of that class exists on that side of the chromosome. |
| Nearest element per side | `up_class` / `dn_class` = the closest of Alu, L1, LTR on that side **if within 10 kb**, otherwise `none`; `up_dist` / `dn_dist` = its distance. |
| Configuration (`config`) | `Alu-Alu` — Alu is the nearest element on both sides; `Alu-mixed` — on one side only; `LTR/L1-hybrid` — on neither side; `incomplete_flank` — one side has no element within 10 kb. |
| `Alu_min_dist` | Nearest Alu on either side, min(`Alu_up`, `Alu_dn`) — continuous Alu-distance variable. |
| Architecture group (`group12`) | CpG-island genes are split by Alu distance as in the table below; non-CpG-island genes are not. Thresholds were fixed from the empirical Alu–TSS distance distribution before any expression data were examined; the upstream/downstream asymmetry reflects the asymmetric retroelement density around TSSs. |

### Architecture groups

| `group12` | CpG island | `config` | Distance rule (bp) |
|---|---|---|---|
| `CGI_AluAlu_le2k` | yes | Alu-Alu | `Alu_min_dist` ≤ 2,000 |
| `CGI_AluAlu_2-4k` | yes | Alu-Alu | 2,000 < `Alu_min_dist` ≤ 4,000 |
| `CGI_AluAlu_4-7/8k` | yes | Alu-Alu | otherwise `Alu_up` ≤ 7,000 or `Alu_dn` ≤ 8,000 |
| `CGI_AluAlu_far` | yes | Alu-Alu | beyond |
| `CGI_AluMixed_le2k` | yes | Alu-mixed | `Alu_min_dist` ≤ 2,000 |
| `CGI_AluMixed_2-5/4k` | yes | Alu-mixed | otherwise `Alu_up` ≤ 5,000 or `Alu_dn` ≤ 4,000 |
| `CGI_AluMixed_far` | yes | Alu-mixed | beyond |
| `CGI_hybrid_le2k` | yes | LTR/L1-hybrid | min(`up_dist`, `dn_dist`) ≤ 2,000 |
| `CGI_hybrid_far` | yes | LTR/L1-hybrid | beyond |
| `nonCGI_Alu-Alu` | no | Alu-Alu | — |
| `nonCGI_Alu-mixed` | no | Alu-mixed | — |
| `nonCGI_LTR/L1-hybrid` | no | LTR/L1-hybrid | — |
| `CGI_incomplete`, `nonCGI_incomplete` | — | incomplete_flank | — |

## Output columns

| Column | Description |
|---|---|
| `ENSG` | Ensembl gene ID without version |
| `gene_id`, `gene_name`, `biotype` | GENCODE v26 gene ID (versioned), symbol, gene type |
| `chrom`, `strand`, `tss` | TSS coordinate (1-based, GTF convention) |
| `tss_source` | `longest` · `MANE_Select` · `Ensembl_canonical` · `longest_fallback` |
| `ENST` | transcript whose 5′ end defines the TSS |
| `cpg_island` | `True` if the TSS overlaps a CpG island |
| `Alu_up`, `Alu_dn`, `L1_up`, `L1_dn`, `LTR_up`, `LTR_dn` | distance (bp) to the nearest element of each class upstream / downstream; empty = none on that side |
| `up_class`, `up_dist`, `dn_class`, `dn_dist` | nearest element class within 10 kb on each side and its distance |
| `Alu_min_dist` | min(`Alu_up`, `Alu_dn`) |
| `config` | `Alu-Alu` · `Alu-mixed` · `LTR/L1-hybrid` · `incomplete_flank` |
| `group12` | architecture group (see above) |

### Group sizes (longest-transcript TSS, 58,137 genes)

| group | n | | group | n |
|---|---|---|---|---|
| CGI_AluAlu_le2k | 5,886 | | nonCGI_Alu-Alu | 13,519 |
| CGI_AluAlu_2-4k | 869 | | nonCGI_Alu-mixed | 14,367 |
| CGI_AluAlu_4-7/8k | 175 | | nonCGI_LTR/L1-hybrid | 14,536 |
| CGI_AluAlu_far | 14 | | nonCGI_incomplete | 1,421 |
| CGI_AluMixed_le2k | 3,271 | | | |
| CGI_AluMixed_2-5/4k | 1,293 | | | |
| CGI_AluMixed_far | 300 | | | |
| CGI_hybrid_le2k | 1,265 | | | |
| CGI_hybrid_far | 371 | | | |
| CGI_incomplete | 850 | | | |

## Notes

* UCSC tracks are 0-based half-open and GTF coordinates are 1-based; distances are computed directly between the 1-based TSS and the UCSC intervals, scoring an element with `start ≤ TSS < end` as overlapping (distance 0).
* The longest-transcript and MANE/canonical definitions select the same transcript for 87 % of genes (67 % of protein-coding genes) and give the same `group12` for 95 % (88 %).

## License

MIT — see `LICENSE`.
