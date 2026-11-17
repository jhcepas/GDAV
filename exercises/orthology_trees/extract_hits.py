#!/usr/bin/env python3
"""Extract the sequences of the BLAST hits from a FASTA file (e.g. the proteome you searched).

Usage:
    python3 extract_hits.py hits.tsv 4proteomes.faa > homologs.faa
    python3 extract_hits.py hits.tsv 4proteomes.faa --evalue 1e-10 > homologs.faa

- hits.tsv: BLAST tabular output (-outfmt 6 ...). Column 2 must be the subject ID (sseqid) and, if you
  use --evalue, the E-value must be the second-to-last column (true for `-outfmt 6` and for
  `-outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore"`).
- 4proteomes.faa: FASTA file (plain or .gz). The ID of a sequence is the first word of its header.

Every hit is written once (BLAST can report several alignments for the same subject), in the order
of the BLAST output. The query is included if it is in the FASTA file (it is its own best hit).
A summary is printed to stderr.

Only the Python standard library is used.
"""
import gzip
import sys


def read_hits(path, max_evalue=None):
    ids = []
    seen = set()
    fh = sys.stdin if path == '-' else open(path)
    for line in fh:
        if not line.strip() or line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if max_evalue is not None and float(f[-2]) > max_evalue:
            continue
        if f[1] not in seen:
            seen.add(f[1])
            ids.append(f[1])
    return ids


def read_fasta(path, wanted):
    seqs = {}
    opener = gzip.open if path.endswith('.gz') else open
    name = None
    with opener(path, 'rt') as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                name = line[1:].split()[0]
                if name in wanted:
                    seqs[name] = []
                else:
                    name = None
            elif name is not None:
                seqs[name].append(line)
    return {k: ''.join(v) for k, v in seqs.items()}


def main():
    args = sys.argv[1:]
    max_evalue = None
    if '--evalue' in args:
        i = args.index('--evalue')
        max_evalue = float(args[i + 1])
        del args[i:i + 2]
    if len(args) != 2:
        sys.exit(__doc__)
    ids = read_hits(args[0], max_evalue)
    seqs = read_fasta(args[1], set(ids))
    missing = [i for i in ids if i not in seqs]
    for i in ids:
        if i in seqs:
            print('>' + i)
            print(seqs[i])
    sys.stderr.write('%d hits, %d sequences written\n' % (len(ids), len(ids) - len(missing)))
    if missing:
        sys.stderr.write('WARNING: not found in %s: %s\n' % (args[1], ', '.join(missing)))


if __name__ == '__main__':
    main()
