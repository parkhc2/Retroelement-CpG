#!/usr/bin/env bash
# Download the public reference tracks used to build the gene architecture table.
# Usage: bash scripts/00_download_raw.sh [data/raw]
set -euo pipefail
RAW=${1:-data/raw}
mkdir -p "$RAW"
cd "$RAW"

# GENCODE v26 gene annotation (GRCh38) -- gene models and TSS
wget -nc https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz
# GENCODE v49 -- only needed for --tss mane (MANE_Select / Ensembl_canonical transcript tags,
# which are absent from v26)
wget -nc https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz
# UCSC hg38 CpG island and RepeatMasker tracks
wget -nc https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz
wget -nc https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz

ls -la
