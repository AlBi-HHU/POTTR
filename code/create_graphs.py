from os import cpu_count
from networkx import is_directed_acyclic_graph, transitive_closure_dag, DiGraph
from multiprocessing import Pool
from collections.abc import Iterable
from collections import defaultdict


def filter_graphs(graphs_dict: dict, k: int):
    def filter_graph(graph, filter_nodes):
        filtered = graph.copy()
        filtered.remove_nodes_from(filter_nodes)
        return filtered

    # count how often a node occurs across patients
    node_dict = defaultdict(lambda: 0)
    for p in graphs_dict:
        for n in set().union(*(g.nodes for g in graphs_dict[p])):
            node_dict[n] += 1
    # remove nodes that are less frequent than k as they would never be part of the solution
    rare_nodes = {node for node, count in node_dict.items() if count < k}

    result = {p: [filter_graph(g, rare_nodes) for g in graphs_dict[p]] for p in graphs_dict}
    return result


def add_attributes(node1, node2, graph):
    if node1 in graph and 'cluster_nodes' in graph.nodes[node1]:
        # make sure that node2 is in cluster of all nodes that are in the cluster with node1
        for n in graph.nodes[node1]['cluster_nodes']:
            if n == node1 or n == node2:
                continue
            graph.nodes[n]['cluster_nodes'].add(node2)
        graph.nodes[node1]['cluster_nodes'].add(node2)
    else:
        graph.add_node(node1, cluster_nodes={node2})


def get_graph_from_line(read_line, name):
    graph = DiGraph()
    graph.add_node('0')

    split_line = read_line.strip('\n').strip().split(',')
    if len(split_line) == 2:
        graph.name = split_line[0]
        line = split_line[1].split(' ')
    else:
        graph.name = str(name)
        line = read_line.strip('\n').strip().split(' ')

    for i, e in enumerate(line):
        # _ is special character for us
        if e.__contains__('->-'):
            a, b = e.split('->-')
            graph.add_edge(a, b)
        elif e.__contains__('-?-'):
            a, b = e.split('-?-')
            # mark nodes as nodes from the same cluster
            add_attributes(a, b, graph)
            add_attributes(b, a, graph)
        elif e.__contains__('-/-'):
            a, b = e.split('-/-')
            graph.add_nodes_from([a, b])
        else:
            # e must be a single node instead of an edge
            graph.add_node(e)

    for node in graph.nodes:
        if node == '0':
            continue
        graph.add_edge('0', node)

    try:
        graph = transitive_closure_dag(graph)
        return graph
    except Exception as e:
        print('Error message: ', e)
        print(graph.name)
        print('Is DAG?', is_directed_acyclic_graph(graph))
        return None


def process_split(split: Iterable):
    evol_id, phylo_tree, line = split
    graph = get_graph_from_line(line, phylo_tree)
    return evol_id, graph


def get_graphs_parallel(trees_to_build: list, num_workers: int=cpu_count()):
    graphs = {x[0]: [] for x in trees_to_build}
    split_size = max(1, len(graphs) // num_workers)

    with Pool(processes=num_workers) as pool:
        for evol_id, graph in pool.imap_unordered(process_split, trees_to_build, chunksize=split_size):
            if graph:
                graphs[evol_id].append(graph)

    print(len(graphs), sum([len(graphs[i]) for i in graphs]))
    return graphs


def get_graphs_single_thread(lines: list, names: list, verbose: bool=False):
    if verbose:
        print('Start single threaded computation of graphs')
    graphs = []
    for i, line in enumerate(lines):
        graph = get_graph_from_line(line, names[i])
        graphs.append(graph)

    return graphs