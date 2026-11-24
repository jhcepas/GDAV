# Orthology from gene trees (Comparative Genomics, Exercise 2, after sessions 3 and 4)

Material for the wiki page
[Comparative Genomics Exercise 2: Orthology from gene trees](https://github.com/jhcepas/GDAV/wiki/Comparative-Genomics-Exercise-2:-Orthology-from-gene-trees).

The proteomes are the same as in Exercise 1 and are **not copied here**: use
[`../orthology_brh/data/4proteomes.faa.gz`](../orthology_brh/data/4proteomes.faa.gz) (Ensembl release 116,
human, African elephant, zebrafish and *Ciona*, one canonical protein per gene; see the README of
`orthology_brh/` for how it was built).

| File | Content |
|---|---|
| `queries/HUMAN_TPH2.faa` | human TPH2, Ensembl ENSP00000329093 (gene ENSG00000139287), 490 aa (Goal 1). Same file as `../orthology_brh/queries/HUMAN_TPH2.faa` |
| `queries/HUMAN_TP53.faa` | human TP53, Ensembl ENSP00000269305 (gene ENSG00000141510), 393 aa, identical to UniProt P04637 (Goal 2) |
| `extract_hits.py` | BLAST tabular hits → FASTA of the hit sequences, taken from the proteome (plain or `.gz`) |
| `root_tree.py` | roots a Newick tree (midpoint, or outgroup given by name fragments) and draws it as text with the support values; can save the rooted tree |
| `backup/` | the alignments and IQ-TREE trees of the wiki page (`TPH2_homologs.aln`, `TP53_homologs.aln` and their `.treefile`), for whoever cannot run MAFFT or IQ-TREE |
| `MD5SUMS` | checksums of the files of this folder |

Both scripts use the **Python 3 standard library only** (no ete3, no Biopython). `root_tree.py` does what
`ete3` did in previous years (`set_outgroup`, `get_midpoint_outgroup`, `ete3 view --text`): ete3 3.1.3 no
longer imports on Python ≥ 3.13 without extra packages, so the exercise does not depend on it. Midpoint
rooting places the root exactly halfway along the longest leaf-to-leaf path; the topology is identical to
ete3's midpoint rooting (checked on both trees of the exercise).

## Use

```bash
conda install -c conda-forge -c bioconda blast mafft "iqtree>=3"     # into your own environment
gunzip -c ../orthology_brh/data/4proteomes.faa.gz > 4proteomes.faa
makeblastdb -dbtype prot -in 4proteomes.faa -out 4proteomes.blastdb
blastp -task blastp -query queries/HUMAN_TPH2.faa -db 4proteomes.blastdb -evalue 0.001 \
       -outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore" > TPH2_homologs.tsv
python3 extract_hits.py TPH2_homologs.tsv 4proteomes.faa > TPH2_homologs.faa
mafft --auto TPH2_homologs.faa > TPH2_homologs.aln
iqtree3 -s TPH2_homologs.aln -m MFP -B 1000 --seed 12345 -T 1
python3 root_tree.py TPH2_homologs.aln.treefile                                   # midpoint
python3 root_tree.py TPH2_homologs.aln.treefile --outgroup 'TH,th,PAH,pah,ENSCING00000006533'
```

## Versions used for the results on the wiki (2026-11-17)

| Tool | Version | Source |
|---|---|---|
| NCBI BLAST+ | 2.17.0 | official Linux x64 binaries |
| MAFFT | 7.526 (2024/Apr/26) | official Linux portable package (mafft.cbrc.jp) |
| IQ-TREE | 3.1.4 | official Linux release (github.com/iqtree/iqtree3) |
| Python | 3.14 (scripts work with any Python ≥ 3.6) | |

Re-run on 2026-11-24 with the conda packages (bioconda `blast` 2.17.0 and `iqtree` 3.1.4, conda-forge `mafft`
7.526): identical alignments and `.treefile`s to those in `backup/`.

Real outputs (same seed, same thread count → identical trees; seeds 1, 2, 3 gave the same topology and
UFBoot values within ±4):

| Tree | Sequences | Alignment columns | Best model (BIC) | Run time (1 thread) |
|---|---|---|---|---|
| TPH (`TPH2_homologs`) | 18 | 546 | Q.INSECT+G4 | ~30 s |
| p53 (`TP53_homologs`) | 22 | 802 | JTT+F+G4 | ~1 min |

`mafft --auto` chose L-INS-i for both files.
