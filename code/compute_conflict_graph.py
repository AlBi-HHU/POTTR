from os import cpu_count
from itertools import combinations
from networkx import Graph, DiGraph, MultiGraph
from numpy import sum
from collections import defaultdict
from multiprocessing import Pool


def process_graph_relations(graph: DiGraph) -> tuple[str, frozenset, set, set, set, set, set]:
    incomparable, cluster = set(), set()
    for node1, node2 in combinations(graph.nodes, 2):
        if graph.has_edge(node1, node2) or graph.has_edge(node2, node1):
            continue
        # cluster attributes are symmetric, checking for node1 is sufficient
        if node1 in graph.nodes[node2].get('cluster_nodes', ()):
            cluster.add(frozenset((node1, node2)))
        # node pair must be incomparable
        else:
            incomparable.add(frozenset((node1, node2)))

    directed_edges = set(graph.edges)
    reversed_edges = {(v, u) for u, v in directed_edges}
    undirected_edges = {frozenset(e) for e in graph.edges}
    return graph.name, frozenset(graph.nodes), directed_edges, reversed_edges, undirected_edges, incomparable, cluster


def get_conflict_graph_for_pair(graph_pair: list):
    name1, N1, D1, R1, U1, I1, C1 = graph_pair[0]
    name2, N2, D2, R2, U2, I2, C2 = graph_pair[1]
    shared = N1 & N2
    # {frozenset(e) for e in D1 & R2}: add conflicts of elements in opposite order
    # | (U1 & I2) | (I1 & U2): add conflicts of ordered vs incomparable elements
    # | (I1 & C2) | (C1 & I2): add conflicts of incomparable vs clustered elements
    conflicts = {frozenset(e) for e in D1 & R2} | (U1 & I2) | (I1 & U2) | (I1 & C2) | (C1 & I2)

    conflict_graph = Graph()
    conflict_graph.name = ':'.join(sorted([name1, name2]))
    conflict_graph.add_nodes_from(shared)
    conflict_graph.add_edges_from(conflicts)

    potential_conflicts = defaultdict()
    for node1, node2 in (U1 & C2):
        edge = (node1, node2) if (node1, node2) in D1 else (node2, node1)
        potential_conflicts[edge] = (conflict_graph.name, name1)
    for node1, node2 in (C1 & U2):
        edge = (node1, node2) if (node1, node2) in D2 else (node2, node1)
        potential_conflicts[edge] = (conflict_graph.name, name2)

    return conflict_graph, potential_conflicts


def collect_potential_conflicts(conflicts: dict, dict_to_update: dict):
    for pcp in conflicts:
        if pcp not in dict_to_update:
            dict_to_update[pcp] = {'labels': set(), 'edge_graph_names': set()}
        dict_to_update[pcp]['labels'].add(conflicts[pcp][0])
        dict_to_update[pcp]['edge_graph_names'].add(conflicts[pcp][1])

    return dict_to_update


def process_split(split: list):
    conflict_graphs = []
    potential_conflicts_dict = dict()
    for graph_pair in split:
        graph, potential_conflicts = get_conflict_graph_for_pair(graph_pair)
        conflict_graphs.append(graph)
        potential_conflicts_dict = collect_potential_conflicts(potential_conflicts, potential_conflicts_dict)

    return conflict_graphs, potential_conflicts_dict


def get_conflict_graphs_parallel(graph_pairs: list, verbose: bool=False, num_workers: int=cpu_count()):
    if verbose:
        print('Create pairwise conflict graphs')

    # list of pairs
    num_workers = num_workers or cpu_count()
    split_size = max(1, len(graph_pairs) // num_workers)

    # list of pairs is split into chunks
    splits = [graph_pairs[i:i + split_size] for i in range(0, len(graph_pairs), split_size)]

    with Pool(processes=num_workers) as pool:
        results = pool.map(process_split, splits)
        pool.close()  # No more tasks will be submitted to the pool
        pool.join()  # Wait for worker processes to finish

    if verbose:
        print('Done computing pairwise conflict graphs')
        print('Start creating union conflict graph')

    flattened_results = []
    flattened_potential_conflicts_dict = defaultdict(lambda: {'labels': set(), 'edge_graph_names': set()})
    for result in results:
        conflict_graphs, potential_conflicts_dict = result
        for conflict_graph in conflict_graphs:
            flattened_results.append(conflict_graph)
        for pcp in potential_conflicts_dict:
            flattened_potential_conflicts_dict[pcp]['labels'].update(potential_conflicts_dict[pcp]['labels'])
            flattened_potential_conflicts_dict[pcp]['edge_graph_names'].update(potential_conflicts_dict[pcp]['edge_graph_names'])

    return flattened_results, flattened_potential_conflicts_dict


def get_conflict_graphs_single_thread(graph_pairs: list, verbose: bool=False):
    if verbose:
        print('Create union conflict graph')

    union_graph = MultiGraph()
    union_edges_dict = dict()

    potential_conflicts_dict = dict()
    for graph_pair in graph_pairs:
        conflict_graph, potential_conflicts = get_conflict_graph_for_pair(graph_pair)
        potential_conflicts_dict = collect_potential_conflicts(potential_conflicts, potential_conflicts_dict)
        key = frozenset((graph_pair[0][0], graph_pair[1][0]))
        union_edges_dict[key] = conflict_graph.edges
        union_graph.add_edges_from(conflict_graph.edges, label=':'.join([graph_pair[0][0], graph_pair[1][0]]))
        union_graph.add_nodes_from(conflict_graph.nodes)

    return union_graph, potential_conflicts_dict


def get_union_conflict_graph(pairwise_conflict_graphs: list, verbose: bool=False):
    if verbose:
        print('Create union conflict graph parallel')
    union_graph = MultiGraph()
    union_edges_dict = dict()

    for conflict_graph in pairwise_conflict_graphs:
        key = frozenset(set((conflict_graph.name).split(':')))
        union_edges_dict[key] = conflict_graph.edges
        union_graph.add_edges_from(conflict_graph.edges, label=conflict_graph.name)
        union_graph.add_nodes_from(conflict_graph.nodes)

    return union_graph


def add_low_frequency_edges_union_graph(union_graph: MultiGraph, potential_conflicts: dict):
    keys = list(potential_conflicts.keys())
    for key in keys:
        node1, node2 = key
        reversed_key = (node2, node1)
        if key in potential_conflicts and reversed_key in potential_conflicts:
            if len(potential_conflicts[key]['edge_graph_names']) > len(
                    potential_conflicts[reversed_key]['edge_graph_names']):
                for label in potential_conflicts[reversed_key]['labels']:
                    union_graph.add_edge(*reversed_key, label=label)
                del potential_conflicts[reversed_key]
            elif len(potential_conflicts[reversed_key]['edge_graph_names']) > len(
                    potential_conflicts[key]['edge_graph_names']):
                for label in potential_conflicts[key]['labels']:
                    union_graph.add_edge(*key, label=label)
                del potential_conflicts[key]
            else:
                print('Same frequency of edges: ', key, len(potential_conflicts[key]['edge_graph_names']), reversed_key,
                      len(potential_conflicts[reversed_key]['edge_graph_names']))

    return


def add_resolution_threshold_edges(union_graph: MultiGraph, potential_conflicts: dict, threshold: int):
    keys = list(potential_conflicts.keys())
    for key in keys:
        if len(potential_conflicts[key]['edge_graph_names']) < threshold:
            for label in potential_conflicts[key]['labels']:
                union_graph.add_edge(*key, label=label)
            del potential_conflicts[key]

    return





