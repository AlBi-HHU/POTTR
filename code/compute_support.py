import os
import argparse
import networkx as nx
from compute_conflict_graph import process_graph_relations


def get_parser():
    """Get parser object for combinatorial vaccine design."""
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', '-i', required=True, dest='input', type=str,
                        help='Input trajectories txt file to compute support')
    parser.add_argument('--output_dir', '-o', required=True, dest='dir', type=str,
                        help='Output directory to store max trajectories and support')
    return parser


def trajectory_key(T):
    # one set per cluster; each member yields the same set, the outer frozenset removes the repeats
    clusters = frozenset(frozenset(cluster | {n})
                         for n, cluster in T.nodes(data='cluster_nodes', default=None) if cluster)
    return frozenset(T.nodes), frozenset(T.edges), clusters


def get_graphs_from_computation(sorted_trajectories):
    rec_traj_graphs = dict()
    for traj in sorted_trajectories:
        graph_name = set(traj.graph.get('name', 'Graph has no name').split(':'))
        key = trajectory_key(traj)
        if key not in rec_traj_graphs:
            rec_traj_graphs[key] = [traj, graph_name]
        else:
            rec_traj_graphs[key][1].update(graph_name)
    return rec_traj_graphs


def find_conflicts(traj_relations, graph_relations):
    """
    Same rules as in the conflict graph computation, applied to the transitively closed trajectory (1) and the
    subgraph of a candidate graph induced by the trajectory nodes (2). Ordered vs clustered pairs are no conflict.
    """
    _, _, D1, _, U1, I1, C1 = traj_relations
    _, _, _, R2, U2, I2, C2 = graph_relations
    return {frozenset(e) for e in D1 & R2} | (U1 & I2) | (I1 & U2) | (I1 & C2) | (C1 & I2)


def compute_support(trajectories: list, input_graphs: dict, output_dir: str):
    sorted_graphs = sorted(trajectories, key=lambda n: n.number_of_nodes(), reverse=True)
    rec_traj_graphs = get_graphs_from_computation(sorted_graphs)

    # relations of the transitively closed trajectories only need to be computed once
    traj_relations = {key: process_graph_relations(nx.transitive_closure_dag(traj))
                      for key, (traj, _) in rec_traj_graphs.items()}

    # iterate over all input graphs and check for any of the trajectories
    for patient in input_graphs:
        for G in input_graphs[patient]:
            graph_ = nx.transitive_closure_dag(G)
            name = G.graph.get('name', 'Graph has no name')
            for key, val in rec_traj_graphs.items():
                traj, graph_names = val
                # only graphs containing all trajectory nodes are candidates
                if not set(traj.nodes) <= set(graph_.nodes):
                    continue
                sub_relations = process_graph_relations(graph_.subgraph(traj.nodes))
                if find_conflicts(traj_relations[key], sub_relations):
                    continue
                # a trajectory cluster that is ordered in a supporting graph is no conflict, but unexpected:
                # the ILP could have selected this graph and resolved the cluster
                unresolved = traj_relations[key][6] & sub_relations[4]
                if unresolved:
                    print(f'Warning: graph {name} supports trajectory {traj.name} but orders its cluster pairs',
                          [tuple(pair) for pair in unresolved])
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
        # graphs of the same patient are not independent: support counts patients, the output still lists all graphs
        # graph names are '<patient>-<tree number>', patient ids may contain '-' themselves
        count = len({name.rsplit('-', 1)[0] for name in graph_names})
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

            file_info.write(str(i) + ',' + str(count) + ',' + ' '.join(graph_names) + ',')
            file_info.write(' '.join([f'{edge}' for edge in edges]) + '\n')
            i += 1

    return file_name