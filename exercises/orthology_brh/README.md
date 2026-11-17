# Orthology by reciprocal BLAST (Comparative Genomics, Exercise 1, after session 2)

Material for the wiki page
[Comparative Genomics Exercise 1: Orthology by reciprocal BLAST](https://github.com/jhcepas/GDAV/wiki/Comparative-Genomics-Exercise-1:-Orthology-by-reciprocal-BLAST).

| File | Content |
|---|---|
| `queries/LOXAF_G3T0S8.faa` | elephant protein G3T0S8 (UniProt) = Ensembl ENSLAFP00000006601, gene ENSLAFG00000007862 (Goal 1) |
| `queries/DANRE_tph1a.faa` | zebrafish tph1a, ENSDARP00000067158 (Goal 2) |
| `queries/CIOIN_TPH.faa` | *Ciona* tryptophan hydroxylase, ENSCINP00000007100, gene ENSCING00000003462 (NCBI Gene 100049066) (Goal 3) |
| `queries/HUMAN_TPH2.faa` | human TPH2, ENSP00000329093 (Goal 3) |
| `aaah_family.faa` | the 18 aromatic amino acid hydroxylases (TPH1, TPH2, TH, PAH and relatives) of the four species + 6 unrelated proteins (ACTB, GAPDH, HBB, INS). Enough to do the exercise without downloading the proteomes (only E-values change) |
| `data/4proteomes.faa.gz` | **the four proteomes, ready to use** (82,853 proteins, 24.6 MB compressed, `gzip -9 -n`) |
| `get_data.sh` | *optional*: how the data were built (downloads Ensembl release 116 and rebuilds `4proteomes.faa`) |
| `build_proteome.py` | keeps one protein per gene (Ensembl canonical), used by `get_data.sh` |
| `coverage.py` | adds query/subject coverage to BLAST tabular output; `--best` keeps the best hit per species |
| `MD5SUMS` | checksums of the files of this folder, of `data/4proteomes.faa.gz` and of the uncompressed `4proteomes.faa` |

All scripts use bash and the Python 3 standard library only.

## The proteome file

`4proteomes.faa`: 82,853 proteins, headers `>SPECIES|gene_symbol|protein_id` (Ensembl gene ID when the gene
has no symbol), one sequence per line.

| Code | Species | Ensembl 116 assembly | Proteins in pep.all | Kept (genes) |
|---|---|---|---|---|
| HUMAN | *Homo sapiens* | GRCh38 | 382,428 | 20,557 |
| LOXAF | *Loxodonta africana* | loxAfr3 (Loxafr3.0, 2009; gene set 2011) | 25,635 | 20,033 |
| DANRE | *Danio rerio* | GRCz11 | 52,089 | 25,592 |
| CIOIN | *Ciona intestinalis* (KH assembly; the species now named *Ciona robusta*) | KH | 17,302 | 16,671 |

How "one protein per gene" is chosen: Ensembl stores a **canonical transcript** for every gene (for human
protein-coding genes it is MANE Select when available). `build_proteome.py` reads it from the core-database
dump of the same release (`gene.txt.gz`, column `canonical_transcript_id`, joined to `transcript.txt.gz`
for the stable ID) and keeps the `pep.all` translation of that transcript. Translations on patches and
alternative haplotypes (human `*_PATCH`, `HSCHR*`; zebrafish `ALT_*`) are removed. The pep headers of
Ensembl do not carry the canonical flag, and the GTF files that do are much larger (141 MB for human).

## Use

```bash
gunzip -c data/4proteomes.faa.gz > 4proteomes.faa
makeblastdb -dbtype prot -in 4proteomes.faa -out 4proteomes.blastdb
```

## Reproduce (optional: how the data were built)

```bash
bash get_data.sh       # downloads ~70 MB into ensembl_116/, checks md5, builds 4proteomes.faa, checks its md5
md5sum -c --ignore-missing MD5SUMS   # checks the files of this folder (and 4proteomes.faa if present)
```

Expected md5 of `4proteomes.faa`: `f66e8193226babaa42efe2e713842768`. Ensembl keeps old releases on its FTP
(`pub/release-116/`), so the script should keep giving the same file.

Results on the wiki were obtained with NCBI BLAST+ 2.17.0 (`makeblastdb -dbtype prot`, `blastp -task
blastp -evalue 0.001`, default parameters otherwise).
