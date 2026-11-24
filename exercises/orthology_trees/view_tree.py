#!/usr/bin/env python3
"""Draw a tree as ASCII text in the terminal, with the support values at the internal nodes.

From the terminal:
    python3 view_tree.py TPH2_homologs.aln.treefile          # draw a Newick file
    python3 view_tree.py TPH_outgroup.nwk --dist             # also show branch lengths

From Python, after editing a tree with ete4:
    import sys; sys.path.append('GDAV/exercises/orthology_trees')
    from view_tree import show
    show(t)                     # t is an ete4 Tree
    show(t, dist=True)

The tree is not modified: children are only drawn smaller clade first (ladderized).
"""
import argparse

from ete4 import Tree


def _label(node, dist):
    if node.is_leaf:
        text = node.name or ''
    else:
        text = '' if node.support is None else '%g' % node.support
    if dist and node.dist is not None and not node.is_root:
        text += ' (%.3f)' % node.dist
    return text


def _draw(node, dist, char='-'):
    """Return (lines, row of the stem) for the subtree of `node` (same layout as ete3's get_ascii)."""
    label = _label(node, dist)
    if node.is_leaf:
        return ['%s-%s' % (char, label)], 0
    width = max(3, len(label) + 2)
    pad = ' ' * width
    lines, stems = [], []
    children = sorted(node.children, key=lambda c: len(c))
    for i, child in enumerate(children):
        c = '/' if i == 0 else ('\\' if i == len(children) - 1 else '-')
        sub, stem = _draw(child, dist, c)
        stems.append(stem + len(lines))
        lines.extend(sub)
        lines.append('')
    lines.pop()
    lo, hi = stems[0], stems[-1]
    prefix = [pad] * (lo + 1) + [pad[:-1] + '|'] * (hi - lo - 1) + [pad] * (len(lines) - hi)
    mid = (lo + hi) // 2
    stem = char + '-' * (width - 2) + prefix[mid][-1]
    prefix[mid] = stem[0] + label + stem[len(label) + 1:]
    return [p + line for p, line in zip(prefix, lines)], mid


def show(tree, dist=False):
    """Print an ete4 tree as ASCII (support values at the internal nodes; dist=True adds branch lengths)."""
    print('\n'.join(_draw(tree, dist)[0]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('tree', help='Newick file (e.g. the .treefile of IQ-TREE)')
    ap.add_argument('--dist', action='store_true', help='show branch lengths')
    args = ap.parse_args()
    show(Tree(open(args.tree)), dist=args.dist)


if __name__ == '__main__':
    main()
