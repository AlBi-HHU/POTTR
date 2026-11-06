"""
Draw trajectories with Graphviz (one vector PDF page per graph).

Nodes of the same cluster ('-?-' relations, unordered alterations) are merged into one grey box listing all
cluster members, the graph is drawn as transitive reduction with root '0', as in draw_trajectory.py.
Alterations in `highlight` are printed in bold red and their boxes are filled orange.
"""
import html
import os
import subprocess
from collections import defaultdict

import networkx as nx

CLUSTER_FILL = '#e6e6e6'
HIGHLIGHT_FILL = '#fdd9a0'
HIGHLIGHT_FONT = '#b30000'


def condense_clusters(graph: nx.DiGraph):
    """Merge each cluster into one node (attribute 'members') and return the transitive reduction."""
    parent = {n: n for n in graph.nodes}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for n, cluster in graph.nodes(data='cluster_nodes'):
        for m in cluster or ():
            if m in parent:
                parent[find(m)] = find(n)

    members = defaultdict(list)
    for n in graph.nodes:
        members[find(n)].append(n)

    condensed = nx.DiGraph()
    for rep, nodes in members.items():
        condensed.add_node(rep, members=sorted(nodes))
    for u, v in graph.edges:
        if find(u) != find(v):
            condensed.add_edge(find(u), find(v))

    if nx.is_directed_acyclic_graph(condensed):
        reduced = nx.transitive_reduction(condensed)
        reduced.add_nodes_from(condensed.nodes(data=True))
        condensed = reduced
    return condensed


def _title_label(lines: list):
    """Left-aligned multi-line graph label."""
    escaped = [line.replace('\\', '\\\\').replace('"', '\\"') for line in lines]
    return '"' + '\\l'.join(escaped) + '\\l"'


def _node_label(members: list, highlight: set, per_line: int):
    parts = []
    for m in members:
        text = html.escape(m, quote=False)
        if m in highlight:
            text = f'<B><FONT COLOR="{HIGHLIGHT_FONT}">{text}</FONT></B>'
        parts.append(text)
    rows = [', '.join(parts[i:i + per_line]) for i in range(0, len(parts), per_line)]
    return '<' + '<BR/>'.join(rows) + '>'


def to_dot(graph: nx.DiGraph, title_lines: list, highlight: set = frozenset(), per_line: int = 4,
           title_fontsize: int = 14):
    """Graphviz source of a graph with root '0' (as returned by create_graphs.get_graph_from_line)."""
    condensed = condense_clusters(graph)
    lines = ['digraph G {',
             f'  graph [label={_title_label(title_lines)}, labelloc=t, labeljust=l, fontsize={title_fontsize}, '
             f'fontname="Helvetica", margin=0.3, nodesep=0.3, ranksep=0.4];',
             '  node [shape=box, style="rounded,filled", fillcolor="white", fontname="Helvetica", fontsize=11];',
             '  edge [arrowsize=0.7];']
    ids = {}
    for i, node in enumerate(condensed.nodes):
        ids[node] = f'n{i}'
        node_members = condensed.nodes[node]['members']
        if node == '0':
            lines.append(f'  n{i} [label="0", shape=circle, fillcolor="#bdbdbd", fontsize=9, width=0.3];')
            continue
        if any(m in highlight for m in node_members):
            fill = HIGHLIGHT_FILL
        elif len(node_members) > 1:
            fill = CLUSTER_FILL
        else:
            fill = 'white'
        lines.append(f'  n{i} [label={_node_label(node_members, highlight, per_line)}, fillcolor="{fill}"];')
    for u, v in condensed.edges:
        lines.append(f'  {ids[u]} -> {ids[v]};')
    lines.append('}')
    return '\n'.join(lines) + '\n'


def message_dot(title_lines: list, title_fontsize: int = 14):
    """Graphviz source of a page that only shows text, e.g. for an input tree that could not be found."""
    return ('digraph G {\n'
            f'  graph [label={_title_label(title_lines)}, labelloc=t, labeljust=l, fontsize={title_fontsize}, '
            f'fontname="Helvetica", margin=0.3];\n'
            '  node [shape=plaintext, label=""];\n  empty;\n}\n')


def render_pdf(dot_source: str, out_pdf: str, dot_binary: str = 'dot'):
    """Write the .dot file next to the pdf and render it with Graphviz."""
    os.makedirs(os.path.dirname(out_pdf), exist_ok=True)
    dot_file = os.path.splitext(out_pdf)[0] + '.dot'
    with open(dot_file, 'w') as f:
        f.write(dot_source)
    subprocess.run([dot_binary, '-Tpdf', '-o', out_pdf, dot_file], check=True)
    return out_pdf
