import os
import itertools
import networkx as nx
import numpy as np
from collections import defaultdict
from scipy.special import comb
from multiprocessing import Pool


def get_conflict_graph_for_pair(graph_pair: list):
    # add nodes and edges from first graph
    nodes = set(graph_pair[0].nodes)
    union_edges = set(graph_pair[0].edges)
    cluster_node_dict = defaultdict(set)
    incomparable_dict = defaultdict(set)
    for graph in graph_pair:
        # if nodes or edges differ from first graph, they are removed from intersection
        nodes.intersection_update(set(graph.nodes))
        union_edges.update(set(graph.edges))

        for node1, node2 in itertools.combinations(graph.nodes, 2):
            if (node1, node2) not in graph.edges and (node2, node1) not in graph.edges:
                # check if node1 and node2 are in a cluster
                if 'cluster_nodes' in graph.nodes[node1]:
                    if node2 in graph.nodes[node1]['cluster_nodes']:
                        cluster_node_dict[node1].add(node2)
                        cluster_node_dict[node2].add(node1)
                    else:
                        incomparable_dict[node1].add(node2)
                        incomparable_dict[node2].add(node1)
                else:
                    incomparable_dict[node1].add(node2)
                    incomparable_dict[node2].add(node1)

    conflict_graph = nx.Graph()
    conflict_graph.name = ':'.join(sorted([graph.name for graph in graph_pair]))
    conflict_graph.add_nodes_from(nodes)
    potential_conflicts = defaultdict()
    for node_pair in itertools.combinations(nodes, 2):
        node1, node2 = node_pair
        conflict_row = [0, 0, 0, 0]

        # if node1 precedes node2 in any graph set row entry 0 to 1
        if (node1, node2) in union_edges:
            conflict_row[0] = 1
        # if node2 precedes node1 in any graph set row entry 1 to 1
        if (node2, node1) in union_edges:
            conflict_row[1] = 1
        # if node1 and node2 are incomparable set row entry 2 to 1
        if node2 in incomparable_dict[node1]:
            conflict_row[2] = 1
        # if node1 and node2 are in the same cluster set row entry 3 to 1
        if node2 in cluster_node_dict[node1]:
            conflict_row[3] = 1

        # add edge if conflict
        summed_3cols = np.sum(conflict_row[:3])
        summed_last2cols = np.sum(conflict_row[-2:])
        if summed_3cols > 1 or summed_last2cols > 1:
            conflict_graph.add_edge(node1, node2)

        # without resolution threshold, we can resolve clusters when observing a directed edge for a node pair;
        # however, if we only want to resolve clusters when we observe a required amount of directed edges, we
        # keep the cluster and edge pairs as a potential conflict and add an edge between the two nodes if the
        # threshold is not met; otherwise, no edge is added to the conflict graph
        if conflict_row[0] and conflict_row[3]:
            edge_graph_name = graph_pair[0].name if (node1, node2) in graph_pair[0].edges else graph_pair[1].name
            potential_conflicts[(node1, node2)] = (conflict_graph.name, edge_graph_name)
        if conflict_row[1] and conflict_row[3]:
            edge_graph_name = graph_pair[0].name if (node2, node1) in graph_pair[0].edges else graph_pair[1].name
            potential_conflicts[(node2, node1)] = (conflict_graph.name, edge_graph_name)

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


def get_conflict_graphs_parallel(graph_pairs: list, verbose: bool=False, num_workers=os.cpu_count()):
    if verbose:
        print('Create pairwise conflict graphs')

    # list of pairs
    if num_workers:
        split_size = max(1, len(graph_pairs) // num_workers)
    else:
        split_size = max(1, len(graph_pairs) // os.cpu_count())

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

    union_graph = nx.MultiGraph()
    union_edges_dict = dict()

    potential_conflicts_dict = dict()
    for graph_pair in graph_pairs:
        conflict_graph, potential_conflicts = get_conflict_graph_for_pair(graph_pair)
        potential_conflicts_dict = collect_potential_conflicts(potential_conflicts, potential_conflicts_dict)
        key = frozenset(set(graph.name for graph in graph_pair))
        union_edges_dict[key] = conflict_graph.edges
        union_graph.add_edges_from(conflict_graph.edges, label=':'.join(graph.name for graph in graph_pair))
        union_graph.add_nodes_from(conflict_graph.nodes)

    return union_graph, potential_conflicts_dict


def get_union_conflict_graph(pairwise_conflict_graphs: list, verbose: bool=False):
    if verbose:
        print('Create union conflict graph parallel')
    union_graph = nx.MultiGraph()
    union_edges_dict = dict()

    for conflict_graph in pairwise_conflict_graphs:
        key = frozenset(set((conflict_graph.name).split(':')))
        union_edges_dict[key] = conflict_graph.edges
        union_graph.add_edges_from(conflict_graph.edges, label=conflict_graph.name)
        union_graph.add_nodes_from(conflict_graph.nodes)

    return union_graph


def add_low_frequency_edges_union_graph(union_graph: nx.MultiGraph, potential_conflicts: dict):
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


def add_resolution_threshold_edges(union_graph: nx.MultiGraph, potential_conflicts: dict, threshold: int):
    keys = list(potential_conflicts.keys())
    for key in keys:
        if len(potential_conflicts[key]['edge_graph_names']) < threshold:
            for label in potential_conflicts[key]['labels']:
                union_graph.add_edge(*key, label=label)
            del potential_conflicts[key]

    return





