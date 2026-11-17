#!/usr/bin/env python3
"""Add coverage columns to BLAST tabular output and (optionally) keep the best hit per species.

Expects BLAST output made with:
    -outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore"
and sequence names like SPECIES|gene|protein (e.g. HUMAN|TPH1|ENSP00000508368).

Usage:
    python3 coverage.py hits.tsv                 # all hits + qcov and scov
    python3 coverage.py hits.tsv --best          # best hit (highest bit score) per species, self-hit removed
    blastp ... | python3 coverage.py - --best    # read from a pipe

qcov = alignment length / query length (%), scov = alignment length / subject length (%).
(Approximate: the alignment length includes gaps. BLAST can also report coverage itself with the
`qcovhsp` column.)

Only the Python standard library is used.
"""
import sys


def main():
    args = sys.argv[1:]
    best = '--best' in args
    files = [a for a in args if a != '--best'] or ['-']
    fh = sys.stdin if files[0] == '-' else open(files[0])
    rows = []
    for line in fh:
        if not line.strip() or line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) != 8:
            sys.exit('Expected 8 columns: qseqid sseqid pident length qlen slen evalue bitscore')
        q, s, pident, length, qlen, slen, evalue, bits = f
        rows.append((q, s, float(pident), int(length), int(qlen), int(slen), evalue, float(bits)))
    if best:
        top = {}
        for r in rows:
            if r[0] == r[1]:
                continue                                  # self-hit
            key = (r[0], r[1].split('|')[0])              # (query, subject species)
            if key not in top or r[7] > top[key][7]:
                top[key] = r
        rows = sorted(top.values(), key=lambda r: (r[0], -r[7]))
    print('\t'.join(['query', 'subject', 'pident', 'length', 'qlen', 'slen', 'qcov%', 'scov%', 'evalue', 'bitscore']))
    for q, s, pident, length, qlen, slen, evalue, bits in rows:
        print('\t'.join([q, s, '%.1f' % pident, str(length), str(qlen), str(slen),
                         '%.0f' % (100.0 * length / qlen), '%.0f' % (100.0 * length / slen), evalue, '%g' % bits]))


if __name__ == '__main__':
    main()
