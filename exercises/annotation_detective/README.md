# Annotation detective (Comparative Genomics, after session 1)

Material for the wiki page
[Comparative Genomics Challenge: Annotation detective](https://github.com/jhcepas/GDAV/wiki/Comparative-Genomics-Challenge:-Annotation-detective).

| File | Content |
|---|---|
| `entries.tsv` | the five UniProtKB entries of the challenge (cases A-E) |
| `annotation_evidence.py` | prints, for each entry, where its name, EC number, function and catalytic activity come from (evidence code, meaning, source) plus its InterPro families. Python 3 standard library only |
| `snapshots/` | the five entries as JSON (UniProt release 2026_04) and the meaning of their evidence codes, so the challenge works offline and gives the same answers after UniProt changes |

```bash
python3 annotation_evidence.py --file entries.tsv            # live data from rest.uniprot.org
python3 annotation_evidence.py --offline --file entries.tsv  # saved snapshots (UniProt 2026_04)
python3 annotation_evidence.py P12345                        # any other UniProtKB entry
```
