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
| `view_tree.py` | draws a Newick tree (or, from Python, an ete4 `Tree`) as ASCII text with the support values at the internal nodes; `--dist` adds branch lengths (needs ete4) |
| `backup/` | the alignments and IQ-TREE trees of the wiki page (`TPH2_homologs.aln`, `TP53_homologs.aln` and their `.treefile`), for whoever cannot run MAFFT or IQ-TREE |
| `MD5SUMS` | checksums of the files of this folder |

`extract_hits.py` uses the Python 3 standard library only. Trees are rooted and saved with
[ete4](https://etetoolkit.github.io/ete/) in Python (code on the wiki page) and drawn with `view_tree.py`
(same text layout as the former `ete3 view --text`). The current conda-forge
package `ete4` 4.4.0 declares `brotli` instead of `brotli-python` and fails on import (`No module named
'brotli'`): install `brotli-python` with it.

## Use

```bash
conda install -c conda-forge -c bioconda blast mafft "iqtree>=3" ete4 brotli-python   # into your own environment
gunzip -c ../orthology_brh/data/4proteomes.faa.gz > 4proteomes.faa
makeblastdb -dbtype prot -in 4proteomes.faa -out 4proteomes.blastdb
blastp -task blastp -query queries/HUMAN_TPH2.faa -db 4proteomes.blastdb -evalue 0.001 \
       -outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore" > TPH2_homologs.tsv
python3 extract_hits.py TPH2_homologs.tsv 4proteomes.faa > TPH2_homologs.faa
mafft --auto TPH2_homologs.faa > TPH2_homologs.aln
iqtree3 -s TPH2_homologs.aln -m MFP -B 1000 --seed 12345 -T 1
python3 view_tree.py TPH2_homologs.aln.treefile          # rooting with ete4: see the wiki page
```

## Versions used for the results on the wiki (2026-11-17)

| Tool | Version | Source |
|---|---|---|
| NCBI BLAST+ | 2.17.0 | official Linux x64 binaries |
| MAFFT | 7.526 (2024/Apr/26) | official Linux portable package (mafft.cbrc.jp) |
| IQ-TREE | 3.1.4 | official Linux release (github.com/iqtree/iqtree3) |
| Python | 3.14 (scripts work with any Python ≥ 3.6) | |
| ete4 | 4.4.0 | conda-forge (+ `brotli-python`) |

Re-run on 2026-11-24 with the conda packages (bioconda `blast` 2.17.0 and `iqtree` 3.1.4, conda-forge `mafft`
7.526): identical alignments and `.treefile`s to those in `backup/`. ete4 rooting gives the same
topologies and support values as the former `root_tree.py`.

Real outputs (same seed, same thread count → identical trees; seeds 1, 2, 3 gave the same topology and
UFBoot values within ±4):

| Tree | Sequences | Alignment columns | Best model (BIC) | Run time (1 thread) |
|---|---|---|---|---|
| TPH (`TPH2_homologs`) | 18 | 546 | Q.INSECT+G4 | ~30 s |
| p53 (`TP53_homologs`) | 22 | 802 | JTT+F+G4 | ~1 min |

`mafft --auto` chose L-INS-i for both files.
