#!/usr/bin/env python3
"""Keep one protein per gene (the Ensembl canonical one) from an Ensembl pep.all FASTA file.

Usage:
    python3 build_proteome.py SPECIES_CODE pep.all.fa.gz gene.txt.gz transcript.txt.gz > out.faa

- gene.txt.gz and transcript.txt.gz are the tables of the Ensembl core database dump of the same
  release (ftp.ensembl.org/pub/release-N/mysql/<species>_core_N_x/). Column 12 of `gene` is
  `canonical_transcript_id` (internal ID); column 14 of `transcript` is its stable ID (ENST...).
  Every gene in Ensembl has exactly one canonical transcript; for protein-coding genes it is the
  "Ensembl canonical" transcript shown in the browser (MANE Select in human).
- Proteins whose canonical transcript is not in the list are dropped (alternative isoforms).
- Proteins on patches and alternative haplotypes (human/zebrafish) are dropped: they are extra
  copies of genes of the primary assembly.
- Output headers: >SPECIES|gene_symbol|protein_id  (gene stable ID when there is no symbol).
  Each sequence is written on ONE line, so `grep -A 1 'HUMAN|TPH1|' 4proteomes.faa` extracts it.

Only the Python standard library is used.
"""
import gzip, re, sys


def canonical_transcripts(gene_file, transcript_file):
    canon_internal = set()
    with gzip.open(gene_file, 'rt', encoding='latin-1') as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if f[10] == '1':                      # is_current
                canon_internal.add(f[11])         # canonical_transcript_id
    stable = set()
    with gzip.open(transcript_file, 'rt', encoding='latin-1') as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if f[0] in canon_internal:
                stable.add(f[13])                 # stable_id (ENST...)
    return stable


NON_PRIMARY = re.compile(r'PATCH|^HSCHR|^ALT_')


def main():
    code, pep, gene_file, transcript_file = sys.argv[1:5]
    canon = canonical_transcripts(gene_file, transcript_file)
    keep, seen_genes, n_in, n_out = False, set(), 0, 0
    out = sys.stdout
    seq = []

    def flush():
        if seq:
            out.write(''.join(seq) + '\n')
            seq.clear()
    with gzip.open(pep, 'rt') as fh:
        for line in fh:
            if line.startswith('>'):
                flush()
                n_in += 1
                f = line[1:].split()
                info = dict(x.split(':', 1) for x in f[1:] if ':' in x and not x.startswith(('chromosome', 'scaffold', 'primary_assembly', 'supercontig')))
                region = f[2]                     # e.g. chromosome:GRCh38:11:...
                tid = info['transcript'].split('.')[0]
                gid = info['gene'].split('.')[0]
                keep = (tid in canon and gid not in seen_genes and not NON_PRIMARY.search(region.split(':')[2]))
                if keep:
                    seen_genes.add(gid)
                    n_out += 1
                    symbol = info.get('gene_symbol', gid).replace('|', '_')
                    out.write('>%s|%s|%s\n' % (code, symbol, f[0].split('.')[0]))
            elif keep:
                seq.append(line.strip())
        flush()
    sys.stderr.write('%s: %d proteins in pep.all, %d kept (one per gene, canonical, primary assembly)\n'
                     % (code, n_in, n_out))


if __name__ == '__main__':
    main()
