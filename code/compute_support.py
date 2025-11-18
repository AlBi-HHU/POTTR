import os
import argparse
import networkx as nx


def get_parser():
    """Get parser object for combinatorial vaccine design."""
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', '-i', required=True, dest='input', type=str,
                        help='Input trajectories txt file to compute support')
    parser.add_argument('--output_dir', '-o', required=True, dest='dir', type=str,
                        help='Output directory to store max trajectories and support')
    return parser


def hash_graph(G):
    tuple_of_edges = tuple(sorted(G.edges))
    tuple_of_nodes = tuple(sorted(G.nodes))
    return hash(tuple_of_edges + tuple_of_nodes)


def get_graphs_from_computation(sorted_trajectories):
    rec_traj_graphs = dict()
    for traj in sorted_trajectories:
        graph_name = set(traj.graph.get('name', 'Graph has no name').split(':'))
        graph_hash = hash_graph(traj)
        if graph_hash not in rec_traj_graphs:
            rec_traj_graphs[graph_hash] = [traj, graph_name]
        else:
            rec_traj_graphs[graph_hash][1].update(graph_name)
    return rec_traj_graphs


def compute_support(trajectories: list, input_graphs: dict, output_dir: str):
    sorted_graphs = sorted(trajectories, key=lambda n: n.number_of_nodes(), reverse=True)
    rec_traj_graphs = get_graphs_from_computation(sorted_graphs)

    # iterate over all input graphs in check for any of the trajectories
    for patient in input_graphs:
        for G in input_graphs[patient]:
            graph_ = nx.transitive_closure_dag(G)
            for key, val in rec_traj_graphs.items():
                traj, graph_names = val
                traj = nx.transitive_closure_dag(traj)
                sub_g = graph_.subgraph(traj.nodes).copy()
                if sub_g.edges == traj.edges:
                    name = G.graph.get('name', 'Graph has no name')
                    rec_traj_graphs[key][1].add(name)

    dir = os.path.join(output_dir, 'processed_graphs/')
    if not os.path.isdir(dir):
        os.mkdir(dir)

    file_name = dir + '/processed_graphs_support.csv'
    file_info = open(file_name, 'w')
    file_info.write('File Index,Support,Supporting Graphs,Edges\n')
    i = 1
    for key, val in rec_traj_graphs.items():
        G, graph_names = val
        trans_G = nx.transitive_closure_dag(G)
        graph_names = sorted(list(set(graph_names)))
        count = len(graph_names)
        node_num = len(G.nodes)
        if node_num > 1:
            edges = set()
            for node in G.nodes:
                if 'cluster_nodes' in G.nodes[node]:
                    # only nodes without adjacency in trajectory are allowed to be part of cluster
                    for cluster_node in G.nodes[node]['cluster_nodes']:
                        if trans_G.has_edge(node, cluster_node) or trans_G.has_edge(cluster_node, node):
                            continue
                        else:
                            edges.add(node + '-?-' + cluster_node)
            for edge in G.edges:
                a, b = edge
                edges.add(a + '->-' + b)

            #print(str(i) + ',' + str(count) + ',' + ' '.join(graph_names) + ',')
            file_info.write(str(i) + ',' + str(count) + ',' + ' '.join(graph_names) + ',')
            file_info.write(' '.join([f'{edge}' for edge in edges]) + '\n')
            i += 1

    return file_name