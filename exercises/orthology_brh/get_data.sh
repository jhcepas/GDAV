#!/usr/bin/env bash
# Build 4proteomes.faa (human, African elephant, zebrafish, Ciona) from Ensembl release 116:
# one protein per gene (the Ensembl canonical one), primary assembly only.
#
# Usage:  bash get_data.sh            (needs curl or wget, python3, md5sum; ~70 MB download)
# Output: 4proteomes.faa  (~83,000 proteins, ~44 MB), headers >SPECIES|gene_symbol|protein_id
#
# Sources (all from ftp.ensembl.org/pub/release-116/):
#   fasta/<species>/pep/<Species>.<assembly>.pep.all.fa.gz       all translations (all isoforms)
#   mysql/<species>_core_116_<v>/gene.txt.gz, transcript.txt.gz  core-database tables; the gene table
#       stores each gene's canonical transcript (column canonical_transcript_id)
# build_proteome.py keeps the protein of the canonical transcript of each gene and drops patches and
# alternative haplotypes. Every input file is checked against the md5 recorded below (computed on
# 2026-11-17), and the final 4proteomes.faa against its expected md5.
set -euo pipefail

RELEASE=116
BASE=https://ftp.ensembl.org/pub/release-${RELEASE}
HERE=$(cd "$(dirname "$0")" && pwd)
DL=ensembl_${RELEASE}
mkdir -p "$DL"

# code  ftp-species-dir     file prefix                  core database
SPECIES="HUMAN homo_sapiens Homo_sapiens.GRCh38 homo_sapiens_core_116_38
LOXAF loxodonta_africana Loxodonta_africana.loxAfr3 loxodonta_africana_core_116_3
DANRE danio_rerio Danio_rerio.GRCz11 danio_rerio_core_116_11
CIOIN ciona_intestinalis Ciona_intestinalis.KH ciona_intestinalis_core_116_3"

# md5 of the Ensembl 116 files used to build the course version
cat > "$DL/MD5SUMS" <<'EOF'
2b76acfd14109da0dc22bdbb93a4e70a  Homo_sapiens.GRCh38.pep.all.fa.gz
b319bb0cdba8e104cb3eba4b9f354077  Loxodonta_africana.loxAfr3.pep.all.fa.gz
bb9339902b58bc0758ae4c959482dc9e  Danio_rerio.GRCz11.pep.all.fa.gz
f32015ef77ea70815476b61ab5add638  Ciona_intestinalis.KH.pep.all.fa.gz
5c26129eb7a784680d2442d65974b64d  homo_sapiens_core_116_38.gene.txt.gz
5ee4afaf6ea7efaf5b94c5b87de888a3  homo_sapiens_core_116_38.transcript.txt.gz
592d02a2055d15f9e54096405615849b  loxodonta_africana_core_116_3.gene.txt.gz
6ecd3cb28630dc6fd9da4856aaa7a814  loxodonta_africana_core_116_3.transcript.txt.gz
ec79305991029d021e671c11d3869e52  danio_rerio_core_116_11.gene.txt.gz
624637e820cdaeebecd455f965427241  danio_rerio_core_116_11.transcript.txt.gz
83a59435aba6e3acc0552d714b049421  ciona_intestinalis_core_116_3.gene.txt.gz
83a49d0560030c06d246764df905ab6d  ciona_intestinalis_core_116_3.transcript.txt.gz
EOF
EXPECTED_MD5_4PROTEOMES=f66e8193226babaa42efe2e713842768

fetch() {  # url output
    [ -s "$2" ] && return 0
    echo "  downloading $1"
    if command -v curl >/dev/null; then curl -sSfL -o "$2.part" "$1"; else wget -q -O "$2.part" "$1"; fi
    mv "$2.part" "$2"
}

echo "Downloading Ensembl release $RELEASE files into $DL/"
while read -r code dir prefix core; do
    fetch "$BASE/fasta/$dir/pep/$prefix.pep.all.fa.gz" "$DL/$prefix.pep.all.fa.gz"
    fetch "$BASE/mysql/$core/gene.txt.gz" "$DL/$core.gene.txt.gz"
    fetch "$BASE/mysql/$core/transcript.txt.gz" "$DL/$core.transcript.txt.gz"
done <<< "$SPECIES"

echo "Checking md5 of the downloaded files"
(cd "$DL" && md5sum -c --quiet MD5SUMS)

echo "Building one-protein-per-gene proteomes"
rm -f 4proteomes.faa.part
while read -r code dir prefix core; do
    python3 "$HERE/build_proteome.py" "$code" "$DL/$prefix.pep.all.fa.gz" \
        "$DL/$core.gene.txt.gz" "$DL/$core.transcript.txt.gz" >> 4proteomes.faa.part
done <<< "$SPECIES"
mv 4proteomes.faa.part 4proteomes.faa

md5=$(md5sum 4proteomes.faa | cut -d' ' -f1)
if [ "$md5" = "$EXPECTED_MD5_4PROTEOMES" ]; then
    echo "OK: 4proteomes.faa ($(grep -c '>' 4proteomes.faa) proteins, md5 $md5) is identical to the course version"
else
    echo "WARNING: md5 of 4proteomes.faa ($md5) differs from the course version ($EXPECTED_MD5_4PROTEOMES)" >&2
    exit 1
fi
echo "You can now delete $DL/ to save space."
