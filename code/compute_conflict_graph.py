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

    node_tuples = []
    conflict_matrix = np.zeros((int(comb(len(nodes), 2)), 4), dtype=int)
    for i, node_pair in enumerate(itertools.combinations(nodes, 2)):
        node1, node2 = node_pair
        node_tuples.append((node1, node2))

        # if node1 precedes node2 in any graph set row entry 0 to 1
        if (node1, node2) in union_edges:
            conflict_matrix[i][0] = 1
        # if node2 precedes node1 in any graph set row entry 1 to 1
        if (node2, node1) in union_edges:
            conflict_matrix[i][1] = 1
        # if node1 and node2 are incomparable set row entry 2 to 1
        if node2 in incomparable_dict[node1]:
            conflict_matrix[i][2] = 1
        # if node1 and node2 are in the same cluster set row entry 3 to 1
        if node2 in cluster_node_dict[node1]:
            conflict_matrix[i][3] = 1

    conflict_graph = nx.Graph()
    conflict_graph.name = ':'.join(sorted([graph.name for graph in graph_pair]))
    conflict_graph.add_nodes_from(nodes)
    summed_3cols = np.sum(conflict_matrix[:, :3], axis=1)
    summed_last2cols = np.sum(conflict_matrix[:, -2:], axis=1)
    for i, val in enumerate(summed_3cols):
        # if summ of columns larger than 1 we have a conflict and need an edge in the conflict graph
        if val > 1 or summed_last2cols[i] > 1:
            conflict_graph.add_edge(*node_tuples[i])

    return conflict_graph, conflict_matrix, node_tuples, union_edges


def get_union_conflict_graph(pairwise_conflict_graphs: list, verbose: bool=False):
    if verbose:
        print('Create union conflict graph parallel')
    union_graph = nx.MultiGraph()
    union_edges_dict = dict()

    for conflict_graph in pairwise_conflict_graphs:
        # TODO check if this works with new implementation of selection
        key = frozenset(set((conflict_graph.name).split(':')))
        union_edges_dict[key] = conflict_graph.edges
        union_graph.add_edges_from(conflict_graph.edges, label=conflict_graph.name)
        union_graph.add_nodes_from(conflict_graph.nodes)

    return union_graph


def process_split(split: list):
    conflict_graphs = []
    for graph_pair in split:
        graph, _, _, _ = get_conflict_graph_for_pair(graph_pair)
        conflict_graphs.append(graph)

    return conflict_graphs


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
    for result in results:
        for conflict_graph in result:
            flattened_results.append(conflict_graph)

    return flattened_results


def get_conflict_graphs_single_thread(graph_pairs: list, verbose: bool=False):
    if verbose:
        print('Create union conflict graph')

    union_graph = nx.MultiGraph()
    union_edges_dict = dict()

    for graph_pair in graph_pairs:
        conflict_graph, _, _, _ = get_conflict_graph_for_pair(graph_pair)
        key = frozenset(set(graph.name for graph in graph_pair))
        union_edges_dict[key] = conflict_graph.edges
        union_graph.add_edges_from(conflict_graph.edges, label=':'.join(graph.name for graph in graph_pair))
        union_graph.add_nodes_from(conflict_graph.nodes)

    return union_graph





