#!/usr/bin/env python3
"""Where does a UniProt annotation come from?

Prints, for each UniProtKB entry, its review status, annotation score, protein existence level, the
protein name and EC numbers WITH the evidence behind each one (ECO code, its meaning and its source),
the function comment (if any) and the protein domains/families matched (InterPro).

Usage:
    python3 annotation_evidence.py P62593 A0A086MRT6          # live data from rest.uniprot.org
    python3 annotation_evidence.py --offline P62593           # use the saved snapshots/ (UniProt 2026_04)
    python3 annotation_evidence.py --file entries.tsv         # accessions from the first column

Only the Python standard library is used.
"""
import json, os, sys, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
_eco_cache = {}


def fetch_json(url):
    with urllib.request.urlopen(urllib.request.Request(url, headers={'Accept': 'application/json'}), timeout=60) as r:
        return json.load(r)


def get_entry(acc, offline=False):
    snap = os.path.join(HERE, 'snapshots', acc + '.json')
    if offline:
        with open(snap) as fh:
            return json.load(fh)
    return fetch_json('https://rest.uniprot.org/uniprotkb/%s.json' % acc)


def eco_label(code, offline=False):
    """Name of an Evidence & Conclusion Ontology term (from EBI QuickGO; cached in snapshots/eco_terms.json)."""
    cache_file = os.path.join(HERE, 'snapshots', 'eco_terms.json')
    if not _eco_cache and os.path.exists(cache_file):
        with open(cache_file) as fh:
            _eco_cache.update(json.load(fh))
    if code not in _eco_cache and not offline:
        try:
            d = fetch_json('https://www.ebi.ac.uk/QuickGO/services/ontology/eco/terms/' + code)
            _eco_cache[code] = d['results'][0]['name']
        except Exception:
            _eco_cache[code] = '?'
    return _eco_cache.get(code, '?')


def evidences(obj, offline=False):
    out = []
    for e in obj.get('evidences', []):
        src = e.get('source', '')
        sid = e.get('id', '')
        out.append('%s (%s)%s' % (e['evidenceCode'], eco_label(e['evidenceCode'], offline),
                                  ' from %s:%s' % (src, sid) if src else ''))
    return out or ['no evidence attached']


def report(acc, offline=False):
    d = get_entry(acc, offline)
    desc = d.get('proteinDescription', {})
    name_block = desc.get('recommendedName') or (desc.get('submissionNames') or [{}])[0]
    print('=' * 100)
    print('%s  %s' % (d['primaryAccession'], d.get('uniProtkbId', '')))
    print('  Database section : %s' % d.get('entryType'))
    print('  Organism         : %s' % d.get('organism', {}).get('scientificName'))
    print('  Length           : %s aa' % d.get('sequence', {}).get('length'))
    print('  Annotation score : %s (1 = lowest, 5 = best annotated)' % d.get('annotationScore'))
    print('  Protein existence: %s' % d.get('proteinExistence'))
    full = name_block.get('fullName', {})
    print('  Protein name     : %s' % full.get('value'))
    for ev in evidences(full, offline):
        print('      evidence     : %s' % ev)
    for ec in name_block.get('ecNumbers', []):
        print('  EC number        : %s' % ec['value'])
        for ev in evidences(ec, offline):
            print('      evidence     : %s' % ev)
    for c in d.get('comments', []):
        if c.get('commentType') == 'FUNCTION':
            for t in c.get('texts', []):
                print('  Function         : %s' % t['value'][:300])
                for ev in evidences(t, offline):
                    print('      evidence     : %s' % ev)
        if c.get('commentType') == 'CATALYTIC ACTIVITY':
            r = c.get('reaction', {})
            print('  Catalytic act.   : %s' % r.get('name', '')[:200])
            for ev in evidences(r, offline):
                print('      evidence     : %s' % ev)
    refs = [x for x in d.get('uniProtKBCrossReferences', []) if x['database'] == 'InterPro']
    for x in refs:
        props = {p['key']: p['value'] for p in x.get('properties', [])}
        print('  InterPro         : %s %s' % (x['id'], props.get('EntryName', '')))
    pubs = [r for r in d.get('references', []) if r.get('citation', {}).get('citationType') == 'journal article']
    print('  Journal articles : %d' % len(pubs))


if __name__ == '__main__':
    args = sys.argv[1:]
    offline = '--offline' in args
    args = [a for a in args if a != '--offline']
    if args and args[0] == '--file':
        with open(args[1]) as fh:
            args = [l.split('\t')[0] for l in fh if l.strip() and not l.startswith('#')]
    if not args:
        sys.exit(__doc__)
    for acc in args:
        report(acc, offline)
