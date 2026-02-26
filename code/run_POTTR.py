#!/usr/bin/env python3
# ---------------------------------------------------------------------------- #
#                                    IMPORTS                                   #
# ---------------------------------------------------------------------------- #
import os
import argparse
import itertools
import POTTR
import compute_support
import convert_to_mastro_format
import networkx as nx
from read_input_dags import read_multiple_graphs_per_evolution
from create_graphs import add_attributes, get_graphs_parallel, get_graphs_single_thread
from compute_conflict_graph import *
from MASTRO_significance_test import compute_significance
from draw_trajectory import draw_trajectory_graph


# ---------------------------------------------------------------------------- #
#                               GLOBAL VARIABLES                               #
# ---------------------------------------------------------------------------- #
verbose: bool
parallel: bool


# ---------------------------------------------------------------------------- #
#                                   FUNCTIONS                                  #
# ---------------------------------------------------------------------------- #
def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-path', '-o', required=True, dest='path', type=str,
                        help='Path to store output files')
    parser.add_argument('--dags', '-d', required=True, dest='dags', type=str,
                        help='File or directory containing transitively closed DAGs (incomplete posets)')
    parser.add_argument('--k', '-k', required=True, dest='k', type=int,
                        help='Number k of incomplete posets to search for common trajectory')
    parser.add_argument('--resolution_threshold', '-rt', dest='resolution_threshold', type=int,
                        help='Number of edges required to resolve a cluster')
    parser.add_argument('--resolution_frequency', '-rf', action='store_true',
                        help='Only allow resolution in the direction of the the most frequent edge for a cluster node pair')
    parser.add_argument('--cores', '-c', dest='cores', type=int, default=1,
                        help='Number cores / threads Gurobi should use; default 0, Gurobi will use all available cores')
    parser.add_argument('--parallelize', '-parallel', action='store_true',
                        help='Enable parallel processing')
    parser.add_argument('--solution-pool-size', '-pool', default=0, dest='pool_size', type=int,
                        help='Solution pool size for Gurobi to retrieve multiple solutions')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Increase output verbosity')
    parser.add_argument('--draw_dots', '-dots', action='store_true',
                        help='Create trajectory png files (only recommended for small instances)')
    return parser


def log(string):
    if verbose:
        print(string)


def label_match(G1, G2):
    if set(G1.nodes) != set(G2.nodes):
        return False

    if set(G1.edges) != set(G2.edges):
        return False

    return True


def filter_duplicates(trajectories):
    """
    ILP might return the same trajectory multiple times but with a different graph selection. Filter out duplicated
    trajectories and join the selected graphs in one output.
    This might cause a combination of trajectories found in multiple trees of the same patient.
    """

    indices_to_remove = []
    for i, G in enumerate(trajectories):
        for j in range(i+1, len(trajectories)):
            # check if graph j is the same as G, if so, remove it
            if nx.is_isomorphic(G, trajectories[j]) and label_match(G, trajectories[j]):
                indices_to_remove.append(j)
                G_name = G.graph.get('name').split(':')
                traj_name = trajectories[j].graph.get('name').split(':')
                new_name = list(set(G_name) | set(traj_name))
                new_name = ':'.join(new_name)
                G.name = new_name

    print('Number of duplicates', len(indices_to_remove))
    unique_trajectories = []
    for i, tra in enumerate(trajectories):
        if i not in indices_to_remove:
            unique_trajectories.append(tra)

    return unique_trajectories


def main():
    args = get_parser().parse_args()
    global verbose
    verbose = args.verbose

    global parallel
    parallel = args.parallelize

    directory = os.path.expanduser(args.path + '/')
    os.makedirs(directory, exist_ok=True)

    log('Start process')

    # read in dags from input path and create all pairwise combinations for computing conflict graphs
    dags = args.dags
    log(f'Reading dags {dags}')

    # graphs will be transformed into a dict of graphs
    graphs_dict = read_multiple_graphs_per_evolution(path=dags, out=directory,
                                                     parallel_processes=args.cores, verbose_flag=verbose)
    # since graphs is a dictionary, the graph pairs for the conflict computation must be computed from the dict values
    graph_pairs = []
    for p1, p2 in itertools.combinations(graphs_dict.keys(), 2):
        for g1, g2 in itertools.product(graphs_dict[p1], graphs_dict[p2]):
            graph_pairs.append([g1, g2])

    k = args.k
    if len(graphs_dict) < k:
        print(f'Value of k is larger than number of patients, set k to maximum value of {len(graphs_dict)}')
        k = len(graphs_dict)

    # create union conflict graph
    if parallel:
        log('Create pairwise conflict graphs parallel')
        pairwise_conflict_graphs, potential_conflicts = get_conflict_graphs_parallel(graph_pairs=graph_pairs, verbose=verbose, num_workers=args.cores)
        log('Compute union conflict graph')
        union_graph = get_union_conflict_graph(pairwise_conflict_graphs=pairwise_conflict_graphs, verbose=verbose)
    else:
        union_graph, potential_conflicts = get_conflict_graphs_single_thread(graph_pairs=graph_pairs, verbose=verbose)

    # add edges if certain clusters should not be resolved, e.g. if confidence of resolution is too low
    if args.resolution_frequency:
        add_low_frequency_edges_union_graph(union_graph, potential_conflicts)
    if args.resolution_threshold:
        add_resolution_threshold_edges(union_graph, potential_conflicts, args.resolution_threshold)
    log('Done creating conflict graph')

    log('Start ILP')
    node_selection_list, graph_selection_list = POTTR.find_max_k_common_trajectory(union_conflict_graph=union_graph,
                                                                                   input_graphs=graphs_dict,
                                                                                   k=k,
                                                                                   cores=args.cores,
                                                                                   solution_pool_size=args.pool_size,
                                                                                   verbose=verbose)

    """
    build trajectory from selected graphs and nodes; since conflict graph has no information about original edges,
    we need to infer them from the selected graphs
    """
    trajectories = []
    for i, graph_list in enumerate(graph_selection_list):
        log(f'Solution {i+1} of {len(graph_selection_list)}')
        max_trajectory = nx.DiGraph()
        max_trajectory.add_node('0')
        names = []
        all_edges = set()
        is_new_order = set()
        for graph in graph_list:
            names.append(graph.name)
            edges = [(u, v) for u, v in graph.edges if u in node_selection_list[i] and v in node_selection_list[i]]
            max_trajectory.add_edges_from(edges)

            """
            identify resolved clusters from original edges: if graph does not have an edge that we added to our 
            trajectory, this must be due to a resolved cluster
            """
            # only do comparison after we added edges to our all_edges set
            if len(all_edges) > 0:
                symm_diff = all_edges ^ set(edges)
                if len(symm_diff) > 0:
                    is_new_order = symm_diff
            all_edges.update(set(edges))

            if args.draw_dots:
                dot = nx.drawing.nx_pydot.to_pydot(graph)
                dot.write_png(directory + 'input_graph_' + graph.name + '.png')

        if len(is_new_order) > 0:
            print('Trajectory introduces order!', names, is_new_order)
            with open(os.path.join(directory, 'resolution.txt'), 'a') as f:
                f.write(graph.name + ' new order ' + ','.join('%s %s' % x for x in list(is_new_order)))

        if len(max_trajectory.nodes) != len(node_selection_list[i]):
            print('Error: Nodes added to max trajectory that were not selected!')
            exit(-1)

        if nx.is_directed_acyclic_graph(max_trajectory):
            max_trajectory = nx.transitive_reduction(max_trajectory)
            max_trajectory.name = ':'.join(names)
        else:
            print('Max Trajectory is no DAG!')
            print(names)
            print(max_trajectory.edges)
            for graph in graph_list:
                if not nx.is_directed_acyclic_graph(graph):
                    print('Selected input graph is no DAG!', graph.name)
            exit(-1)

        for node_pair in list(itertools.combinations(max_trajectory.nodes, 2)):
            a, b = node_pair
            if a in max_trajectory.neighbors(b) or b in max_trajectory.neighbors(a):
                continue
            for g in graph_list:
                if 'cluster_nodes' in g.nodes[a]:
                    if b in g.nodes[a]['cluster_nodes']:
                        add_attributes(a, b, max_trajectory)
                        add_attributes(b, a, max_trajectory)

        trajectories.append(max_trajectory)

    # filter out duplicate results reported by ILP
    trajectories = filter_duplicates(trajectories)
    log('Compute support')
    support_file = compute_support.compute_support(trajectories, graphs_dict, directory)
    log('Convert to output format')
    converted_file = convert_to_mastro_format.convert(support_file, directory)

    gexf_directory = os.path.join(directory, 'trajectories_gexf/')
    os.makedirs(gexf_directory, exist_ok=True)
    with open(os.path.join(gexf_directory, 'traj_graphs_names.csv'), 'w') as f:
        f.write('trajectory,input graphs\n')
        for i, traj in enumerate(trajectories):
            traj_copy = nx.DiGraph()
            traj_copy.add_edges_from(traj.edges())
            traj_copy.add_nodes_from(traj.nodes())
            mapping = {n: ','.join(sorted(list(v) + [n])) for n, attrs in traj.nodes(data=True) for v in attrs.values()}
            traj_copy = nx.relabel_nodes(traj_copy, mapping)
            # must convert node attributes, which are sets, to lists to write gexf
            # nx.set_node_attributes(traj, {n: {k: list(v) for k, v in attr.items()} for n, attr in traj.nodes(data=True)})
            nx.write_gexf(traj_copy, gexf_directory + str(i) + '_trajectory.gexf')
            f.write(str(i) + ',' + traj.name + '\n')

    if args.draw_dots:
        log('Draw trajectories')
        draw_trajectory_graph(converted_file, directory)

    # run MASTRO significance test if trajectories are small enough, e.g. fewer than 10 nodes
    if len(trajectories[0].nodes) < 13:
        log('Run significance test')
        results_significance = os.path.join(directory, 'significance_output.txt')
        compute_significance.run_stat_significancce_test(support_file=converted_file, graph_file=dags,
                                                         output_file=results_significance, cores=args.cores)
    else:
        print('Trajectory size is too large to execute the significance test.')

    trajectory_size = len(list(trajectories[0].nodes))
    return trajectory_size


if __name__ == '__main__':
    main()