#!/usr/bin/env python3
"""Root a Newick tree (midpoint or outgroup) and draw it as text in the terminal.

Usage:
    python3 root_tree.py tree.treefile                          # midpoint rooting, text drawing
    python3 root_tree.py tree.treefile --outgroup 'PAH,TH'      # outgroup rooting
    python3 root_tree.py tree.treefile --unrooted               # draw the tree as it is in the file
    python3 root_tree.py tree.treefile -o rooted.nwk            # also save the rooted tree (Newick)
    python3 root_tree.py tree.treefile --dist                   # show branch lengths in the drawing

--outgroup takes a comma-separated list of text fragments: every leaf whose name contains any of them
(case-sensitive) is in the outgroup. The outgroup must form one clade in the unrooted tree; the root is
placed in the middle of the branch that separates it from the rest.

Midpoint rooting places the root halfway along the longest leaf-to-leaf path of the tree.

The drawing shows the support values of IQ-TREE (UFBoot, or SH-aLRT/UFBoot) at the internal nodes.
Children are sorted so that the smaller clade is drawn first (ladderized).

Only the Python standard library is used (it does the same as `ete3`'s set_outgroup /
get_midpoint_outgroup and `ete3 view --text`, without installing anything).
"""
import argparse
import sys


# ---------------------------------------------------------------- Newick parsing (into an unrooted graph)
class Node:
    __slots__ = ('name', 'children', 'dist', 'support')

    def __init__(self, name='', dist=0.0, support=''):
        self.name, self.children, self.dist, self.support = name, [], dist, support


def parse_newick(text):
    """Return the rooted Node structure exactly as written in the file."""
    text = text.strip()
    if not text.endswith(';'):
        sys.exit('ERROR: the Newick tree must end with ";"')
    pos = 0

    def read_label():
        nonlocal pos
        if text[pos] == "'":
            end = text.index("'", pos + 1)
            label, pos = text[pos + 1:end], end + 1
            return label
        start = pos
        while text[pos] not in ',():;':
            pos += 1
        return text[start:pos].strip()

    def read_node():
        nonlocal pos
        node = Node()
        if text[pos] == '(':
            pos += 1
            while True:
                node.children.append(read_node())
                if text[pos] == ',':
                    pos += 1
                elif text[pos] == ')':
                    pos += 1
                    break
                else:
                    sys.exit('ERROR: unexpected character %r in the Newick tree' % text[pos])
        label = read_label()
        if node.children:
            node.support = label          # internal labels are support values (IQ-TREE)
        else:
            node.name = label
        if text[pos] == ':':
            pos += 1
            start = pos
            while text[pos] not in ',();':
                pos += 1
            node.dist = float(text[start:pos])
        return node

    return read_node()


def to_graph(root):
    """Unrooted graph: adjacency {node: {neighbour: (length, support)}}. Degree-2 roots are removed."""
    adj = {}

    def link(a, b, length, support):
        adj.setdefault(a, {})[b] = (length, support)
        adj.setdefault(b, {})[a] = (length, support)

    def walk(node):
        adj.setdefault(node, {})
        for c in node.children:
            link(node, c, c.dist, c.support if c.children else '')
            walk(c)
    walk(root)
    if len(root.children) == 2:          # rooted input: merge the two branches below the root
        a, b = root.children
        (la, sa), (lb, sb) = adj[root][a], adj[root][b]
        del adj[a][root], adj[b][root], adj[root]
        link(a, b, la + lb, sa or sb)
    return adj


def leaves(adj):
    return [n for n in adj if not n.children]


# ---------------------------------------------------------------- rooting
def rooted_from_edge(adj, a, b, frac):
    """Build a rooted tree with the root on edge a-b, at fraction `frac` of its length from a."""
    length, support = adj[a][b]
    root = Node()

    def build(node, parent, dist, sup):
        new = Node(node.name, dist, sup)
        for nb, (l, s) in adj[node].items():
            if nb is not parent:
                new.children.append(build(nb, node, l, s))
        return new
    root.children = [build(a, b, length * frac, support), build(b, a, length * (1 - frac), support)]
    return root


def paths_from(adj, start):
    """Distance and parent pointers from `start` to every node."""
    dist, parent, stack = {start: 0.0}, {start: None}, [start]
    while stack:
        n = stack.pop()
        for nb, (l, _) in adj[n].items():
            if nb not in dist:
                dist[nb], parent[nb] = dist[n] + l, n
                stack.append(nb)
    return dist, parent


def midpoint(adj):
    lv = leaves(adj)
    best = (-1, None, None)
    for a in lv:
        d, _ = paths_from(adj, a)
        for b in lv:
            if d[b] > best[0]:
                best = (d[b], a, b)
    total, a, b = best
    d, parent = paths_from(adj, a)
    half = total / 2.0
    node = b                                  # walk back from b towards a
    while parent[node] is not None:
        p = parent[node]
        if d[p] <= half <= d[node]:
            length = adj[p][node][0]
            frac = (half - d[p]) / length if length > 0 else 0.5
            return rooted_from_edge(adj, p, node, frac), (a.name, b.name, total)
        node = p
    sys.exit('ERROR: could not find the midpoint')


def leaf_names(node):
    if not node.children:
        return {node.name}
    return set().union(*(leaf_names(c) for c in node.children))


def outgroup_root(adj, fragments):
    lv = leaves(adj)
    out = {n.name for n in lv if any(f in n.name for f in fragments)}
    if not out:
        sys.exit('ERROR: no leaf name contains any of: %s' % ', '.join(fragments))
    if len(out) == len(lv):
        sys.exit('ERROR: all leaves match the outgroup')
    start = next(n for n in lv if n.name not in out)
    tmp = rooted_from_edge(adj, start, next(iter(adj[start])), 0.0)
    best = None                               # smallest clade containing the whole outgroup

    def find(node):
        nonlocal best
        names = leaf_names(node)
        if out <= names and (best is None or len(names) < len(best[1])):
            best = (node, names)
        for c in node.children:
            find(c)
    find(tmp)
    if best[1] != out:
        extra = sorted(best[1] - out)
        sys.exit('ERROR: the outgroup is not a clade in this tree. The smallest clade that contains it '
                 'also contains:\n  ' + '\n  '.join(extra))
    # locate the corresponding edge in the graph: the outgroup clade hangs from one graph edge
    for a in adj:
        for b in adj[a]:
            if a.children or b.children:
                side = side_leaves(adj, b, a)
                if side == out:
                    return rooted_from_edge(adj, a, b, 0.5)
    sys.exit('ERROR: outgroup edge not found')


def side_leaves(adj, node, parent):
    """Leaf names reachable from `node` without going through `parent`."""
    names, stack = set(), [(node, parent)]
    while stack:
        n, p = stack.pop()
        if not n.children:
            names.add(n.name)
        for nb in adj[n]:
            if nb is not p:
                stack.append((nb, n))
    return names


# ---------------------------------------------------------------- output
def size(node):
    return 1 if not node.children else sum(size(c) for c in node.children)


def ladderize(node):
    for c in node.children:
        ladderize(c)
    node.children.sort(key=size)


def to_newick(node, root=True):
    if node.children:
        s = '(' + ','.join(to_newick(c, False) for c in node.children) + ')' + node.support
    else:
        s = node.name
    if not root:
        s += ':%.6g' % node.dist
    return s + (';' if root else '')


def ascii_art(node, show_dist, char1='-'):
    """Text drawing (same layout as ete3's get_ascii), with support values on internal nodes."""
    if node.children:
        label = node.support
    else:
        label = node.name
    if show_dist and char1 != '=':
        label += ' (%.3f)' % node.dist
    if not node.children:
        return ['%s-%s' % (char1 if char1 != '=' else '-', label)], 0
    LEN = max(3, len(label) + 2)
    PAD, PA = ' ' * LEN, ' ' * (LEN - 1)
    result, mids = [], []
    for i, c in enumerate(node.children):
        char2 = '/' if i == 0 else ('\\' if i == len(node.children) - 1 else '-')
        lines, mid = ascii_art(c, show_dist, char2)
        mids.append(mid + len(result))
        result.extend(lines)
        result.append('')
    result.pop()
    lo, hi, end = mids[0], mids[-1], len(result)
    prefixes = [PAD] * (lo + 1) + [PA + '|'] * (hi - lo - 1) + [PAD] * (end - hi)
    mid = (lo + hi) // 2
    first = '-' if char1 == '=' else char1
    prefixes[mid] = first + '-' * (LEN - 2) + prefixes[mid][-1]
    if label:
        stem = prefixes[mid]
        prefixes[mid] = stem[0] + label + stem[len(label) + 1:]
    return [p + l for p, l in zip(prefixes, result)], mid


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('tree', help='Newick file (e.g. the .treefile of IQ-TREE); "-" reads stdin')
    g = ap.add_mutually_exclusive_group()
    g.add_argument('--outgroup', help='comma-separated name fragments defining the outgroup')
    g.add_argument('--unrooted', action='store_true', help='do not re-root: draw the tree as in the file')
    ap.add_argument('-o', '--output', help='write the (re)rooted tree in Newick to this file')
    ap.add_argument('--dist', action='store_true', help='show branch lengths in the drawing')
    args = ap.parse_args()

    text = sys.stdin.read() if args.tree == '-' else open(args.tree).read()
    tree = parse_newick(text)
    if args.unrooted:
        info = 'tree drawn as written in the file (an unrooted tree drawn from an arbitrary point)'
    else:
        adj = to_graph(tree)
        if args.outgroup:
            frags = [f for f in args.outgroup.split(',') if f]
            tree = outgroup_root(adj, frags)
            info = 'rooted with the outgroup: ' + ', '.join(frags)
        else:
            tree, (a, b, total) = midpoint(adj)
            info = 'midpoint rooting: longest path %.3f substitutions/site, between\n  %s and\n  %s' % (total, a, b)
    ladderize(tree)
    print('# ' + info)
    print('\n'.join(ascii_art(tree, args.dist, '=')[0]))
    if args.output:
        with open(args.output, 'w') as fh:
            fh.write(to_newick(tree) + '\n')
        print('# rooted tree saved in %s' % args.output)


if __name__ == '__main__':
    main()
