"""
Read the output of get_output_sizes.py (one run_POTTR.py output folder per k).

Expected layout of the get_output_sizes.py output directory (<output-path><timestamp>/):
    out_k2/converted_graphs.txt
    out_k2/significance_output.txt      (only if run_POTTR.py ran the significance test)
    out_k3/...
    graph_orders.csv
Folders created with a solution pool are named out_k_pool<N>_<k>; both namings are recognized.
"""
import ast
import os
import re
import sys
from dataclasses import dataclass

import networkx as nx

# reuse the parsers of the POTTR code base (Code/ is the parent folder of this script's folder)
CODE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if CODE_DIR not in sys.path:
    sys.path.insert(0, CODE_DIR)

from create_graphs import get_graph_from_line  # noqa: E402


# matches out_k2, out_k_pool10_2, tracerx_k7, lung_k12, ...
K_FOLDER_PATTERN = re.compile(r'(?:^|_)k(?:_pool\d+_)?(\d+)')

# run_POTTR.py only runs the significance test if the trajectory has fewer than 13 nodes (root included)
NO_SIGNIFICANCE_NOTE = ('not computed: no significance_output.txt (run_POTTR.py only runs the test '
                        'for trajectories with fewer than 13 nodes incl. root)')


@dataclass
class Trajectory:
    k: int
    index: int                      # 1-based position in converted_graphs.txt
    relations: list                 # relations as written in converted_graphs.txt, e.g. 'A->-B', 'A-?-B'
    support: int                    # number of supporting patients
    supporting_graphs: list         # names of the supporting input trees, e.g. 'CRUK0063-0'
    graph: nx.DiGraph               # transitively closed trajectory incl. root '0', clusters as 'cluster_nodes'
    pval_ind: float = None          # MASTRO test, independent model (column pval_ind)
    pval_pottr: float = None        # POTTR test (column pval_POTTR)
    evidence: float = None          # evidence score for cluster resolution (column resolution_evidence)
    significance_note: str = ''
    number: int = 0                 # 1-based position among the unique trajectories of the summary
    found_for_k: tuple = ()         # all k for which this trajectory was reported

    @property
    def alterations(self):
        return set(self.graph.nodes) - {'0'}

    @property
    def size(self):
        return len(self.alterations)


def find_k_folders(results_dir: str):
    """Return [(k, folder)] sorted by k for all run_POTTR.py output folders in results_dir."""
    results_dir = os.path.abspath(results_dir)
    if os.path.isfile(os.path.join(results_dir, 'converted_graphs.txt')):
        # a single run_POTTR.py output folder was given
        match = K_FOLDER_PATTERN.search(os.path.basename(results_dir))
        return [(int(match.group(1)) if match else 0, results_dir)]

    folders = []
    for entry in os.listdir(results_dir):
        path = os.path.join(results_dir, entry)
        match = K_FOLDER_PATTERN.search(entry)
        if match and os.path.isdir(path) and os.path.isfile(os.path.join(path, 'converted_graphs.txt')):
            folders.append((int(match.group(1)), path))
    return sorted(folders)


def read_trajectories(converted_file: str, k: int):
    """Parse converted_graphs.txt: line 2i holds 'relations (support)', line 2i+1 the supporting trees."""
    with open(converted_file, 'r') as f:
        lines = [line.rstrip('\n') for line in f]
    while lines and not lines[-1].strip():
        lines.pop()

    trajectories = []
    for i in range(0, len(lines) - 1, 2):
        traj_line = lines[i].strip()
        pos = traj_line.rfind('(')
        relations = traj_line[:pos].split()
        support = int(traj_line[pos + 1:].replace(')', '').strip())
        supporting_graphs = lines[i + 1].split()
        index = len(trajectories) + 1
        graph = get_graph_from_line(' '.join(relations), f'k{k}_trajectory_{index}')
        trajectories.append(Trajectory(k=k, index=index, relations=relations, support=support,
                                       supporting_graphs=supporting_graphs, graph=graph))
    return trajectories


def read_significance(significance_file: str):
    """Rows of significance_output.txt as dicts (';'-separated, every row ends with a trailing ';')."""
    with open(significance_file, 'r') as f:
        header = f.readline().strip().split(';')
        rows = []
        for line in f:
            if line.strip():
                rows.append(dict(zip(header, line.rstrip('\n').split(';'))))
    return header, rows


def _row_graphs(row):
    try:
        return {str(g) for g in ast.literal_eval(row.get('traj_occ_list', '[]'))}
    except (ValueError, SyntaxError):
        return set()


def _row_edges(row):
    edges = row.get('edges_traj', '').strip().strip('[]')
    return {tuple(e.split('->-', 1)) for e in edges.split(',') if '->-' in e}


def _to_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def assign_pvalues(trajectories: list, significance_file: str):
    """
    Add pval_ind (MASTRO, independent), pval_POTTR and the resolution evidence to the trajectories. compute_significance.py writes one row
    per trajectory in the order of converted_graphs.txt; the order is checked with the supporting trees and, if it
    does not match, the row is searched by supporting trees and edges.
    """
    if not os.path.isfile(significance_file):
        for traj in trajectories:
            traj.significance_note = NO_SIGNIFICANCE_NOTE
        return

    header, rows = read_significance(significance_file)
    for traj in trajectories:
        graphs = set(traj.supporting_graphs)
        row = None
        i = traj.index - 1
        if i < len(rows) and _row_graphs(rows[i]) == graphs:
            row = rows[i]
        else:
            candidates = [r for r in rows if _row_graphs(r) == graphs]
            exact = [r for r in candidates if _row_edges(r) == set(traj.graph.edges)]
            if exact:
                row = exact[0]
            elif len(candidates) == 1:
                row = candidates[0]

        if row is None:
            traj.significance_note = 'no matching row in significance_output.txt'
            continue

        traj.pval_ind = _to_float(row.get('pval_ind'))
        traj.pval_pottr = _to_float(row.get('pval_POTTR'))
        traj.evidence = _to_float(row.get('resolution_evidence'))
        if 'pval_POTTR' not in header:
            traj.significance_note = 'significance_output.txt has no pval_POTTR column (older output)'


def load_results(results_dir: str, max_per_k: int = None):
    """Return {k: [Trajectory]} for all k folders, with p-values assigned."""
    results = {}
    for k, folder in find_k_folders(results_dir):
        trajectories = read_trajectories(os.path.join(folder, 'converted_graphs.txt'), k)
        if max_per_k:
            trajectories = trajectories[:max_per_k]
        assign_pvalues(trajectories, os.path.join(folder, 'significance_output.txt'))
        results[k] = trajectories
    return results


def trajectory_key(t: Trajectory):
    """
    Identity of a trajectory independent of k: alterations, edges of the transitive closure and cluster pairs.
    The graph is transitively closed, so the same trajectory written with or without transitive relations matches.
    """
    if t.graph is None:
        return 'relations', frozenset(t.relations)
    clusters = frozenset(frozenset((n, m)) for n, cluster in t.graph.nodes(data='cluster_nodes')
                         for m in (cluster or ()) if m != n)
    return frozenset(t.graph.nodes), frozenset(t.graph.edges), clusters


def unique_trajectories(results: dict):
    """
    Merge trajectories that were reported for several k. Each trajectory is kept once, represented by the
    occurrence with the highest support (then most supporting trees, then one with p-values, then smallest k).
    The result is sorted by support and size (both descending) and numbered from 1.
    """
    groups = {}
    for k, trajectories in sorted(results.items()):
        for t in trajectories:
            groups.setdefault(trajectory_key(t), []).append(t)

    unique = []
    for occurrences in groups.values():
        best = max(occurrences, key=lambda t: (t.support, len(t.supporting_graphs),
                                               t.pval_ind is not None or t.pval_pottr is not None,
                                               -t.k, -t.index))
        best.found_for_k = tuple(sorted({t.k for t in occurrences}))
        unique.append(best)

    unique.sort(key=lambda t: (-t.support, -t.size, t.k, t.index))
    for number, t in enumerate(unique, start=1):
        t.number = number
    return unique
