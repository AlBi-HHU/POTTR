import networkx as nx
from typing import Any
from math import comb, factorial
from networkx.algorithms.isomorphism.vf2userfunc import DiGraphMatcher
from MASTRO_significance_test import compute_significance
from read_input_dags import read_multiple_graphs_per_evolution
from compute_support import process_graph_relations, find_conflicts


def expand(graph: nx.DiGraph):
    expanded = set()
    expanded_graph = graph.copy()
    for node, cluster in graph.nodes(data='cluster_nodes'):
        if not cluster:
            continue
        if node in expanded:
            continue
        curr = node
        expanded.add(curr)
        cluster = list(cluster)
        for i in range(len(cluster)):
            next = cluster[i]
            expanded_graph.add_edge(curr, next)
            expanded.add(next)
            curr = next

    expanded_graph = nx.transitive_closure_dag(expanded_graph)
    return expanded_graph


def compute_traj_order_prob(trajectory: nx.DiGraph, expanded_traj: nx.DiGraph) -> float:
    # get the nodes in hierarchical order using the outdegree for sorting
    degrees = dict(trajectory.out_degree(trajectory.nodes))
    sorted_degrees = dict(sorted(degrees.items(), key=lambda item: item[1], reverse=True))

    # get clusters in according order, treat singletons as cluster
    sorted_clusters = []
    seen = set()
    for node in sorted_degrees:
        # ignore root node
        if node in seen or node == '0':
            continue
        cluster_nodes = trajectory.nodes[node].get('cluster_nodes', ())
        if cluster_nodes:
            cluster = {node}
            cluster.update(cluster_nodes)
            sorted_clusters.append(cluster)
            seen.update(cluster_nodes)
        else:
            sorted_clusters.append({node})
            seen.add(node)

    # compute number of potential label orderings along the trajectorie's topology
    cluster_sizes = [len(c) for c in sorted_clusters]
    # ignore root node
    traj_labels = len(trajectory.nodes) - 1
    possible_orders = 1
    for i in range(len(sorted_clusters)):
        n = traj_labels - sum(cluster_sizes[:i])
        k = cluster_sizes[i]
        possible_orders *= comb(n, k)

    # compute probability that labels of the trajectory are observed in exactly this order
    n_automorphisms = compute_significance.compute_num_automorph(expanded_traj)
    prob = n_automorphisms / possible_orders
    return prob


def compute_stat_significance_POTTR(trajectory: nx.DiGraph, graphs_dict: dict) -> tuple[list[Any], int, list[Any]]:
    expanded_trajectory = expand(trajectory)
    prob_label_order = compute_traj_order_prob(trajectory, expanded_trajectory)
    per_patient_probs = []
    candidate_graphs = []
    for patient in graphs_dict:
        tmp_patient_props = []
        for i, graph_ in enumerate(graphs_dict[patient]):
            if set(expanded_trajectory.nodes).issubset(set(graph_.nodes)):
                candidate_graphs.append(graph_.name)
                expanded_graph = expand(graph_)

                # get number of alterations
                t = len(list(expanded_graph.nodes)) - 1
                k = len(list(expanded_trajectory.nodes)) - 1
                set_g_nodes_nogerm = set(expanded_graph.nodes)
                set_g_nodes_nogerm.remove('0')

                # compute number of subsets of nodes of graph isomorphic to trajectory
                subsets = compute_significance.findsubsets(set_g_nodes_nogerm, k)
                num_isomorph = 0

                for subset_g in subsets:
                    subset_g = list(subset_g)
                    subset_g.append('0')
                    sub_g = expanded_graph.subgraph(subset_g).copy()
                    if nx.is_isomorphic(sub_g, expanded_trajectory):
                        num_isomorph += 1

                # compute probability to observe the trajectory anywhere in expanded_graph
                # probability for label order * probability for a certain position
                prob_position = prob_label_order / comb(t, k)
                # trajectory can only be at one of the isomorphic positions
                prob_graph = num_isomorph * prob_position
                tmp_patient_props.append(prob_graph)

        if tmp_patient_props:
            per_patient_probs.append(sum(tmp_patient_props) / len(tmp_patient_props))

    n = len(graphs_dict.keys())
    return per_patient_probs, n, candidate_graphs


def compute_support_probs(trajectory: nx.DiGraph, graphs_dict: dict):
    def get_clusters(graph: nx.DiGraph):
        cluster_graph = nx.Graph()
        cluster_graph.add_nodes_from(graph.nodes)
        for node1, cluster in graph.nodes(data='cluster_nodes'):
            if cluster:
                for node2 in cluster:
                    cluster_graph.add_edge(node1, node2)

        clusters = []
        for nodes in nx.connected_components(cluster_graph):
            clusters.append(frozenset(nodes))

        clusters.remove({'0'})
        return clusters


    def compute_probs_cluster_fit(trajectory, graph):
        clusters_graph = get_clusters(graph)
        clusters_trajectory = get_clusters(trajectory)
        traj_nodes = set(trajectory.nodes)

        # We are computing the inverted multinomial coefficient here as we want the probability for one possible partition
        # A multinomial coefficient counts the number of ways to partition a set of n distinct objects into k ordered
        # subsets of sizes n_1, n_2, ..., n_k
        p_G = 1
        for cl_graph in clusters_graph:
            b_i = 1
            for cl_traj in clusters_trajectory:
                b_i *= factorial(len(cl_graph & cl_traj))
            p_G *= (b_i / factorial(len(cl_graph & traj_nodes)))

        return p_G


    traj_relations = process_graph_relations(nx.transitive_closure_dag(trajectory))
    traj_edges = {e for e in traj_relations[4] if '0' not in e}
    patient_probs = []
    edges_supp = []
    for patient in graphs_dict:
        probs_true_fit = []
        patient_edge_supp = []
        for G in graphs_dict[patient]:
            # only graphs containing all trajectory nodes are candidates
            if not set(trajectory.nodes) <= set(G.nodes):
                continue
            # check for any conflicts: a supporting graph is not supposed to have any
            graph_ = nx.transitive_closure_dag(G)
            subgraph = graph_.subgraph(trajectory.nodes)
            sub_relations = process_graph_relations(subgraph)
            sub_edges = {e for e in sub_relations[4] if '0' not in e}
            if find_conflicts(traj_relations, sub_relations):
                continue
            # compute direct edge support
            edge_support = (sum(1 for e in traj_edges if e in sub_edges) / len(traj_edges) if traj_edges else 1.0)
            patient_edge_supp.append(edge_support)
            # compute probability that the order of Gs clusters fits trajectory
            prob_supporting_T = compute_probs_cluster_fit(trajectory, subgraph)
            probs_true_fit.append(prob_supporting_T)

        if len(probs_true_fit) > 0:
            # average over all patient trees, given that all trees are equally likely
            avg_prob = sum(probs_true_fit) / len(probs_true_fit)
            patient_probs.append(avg_prob)

            avg_supp_patient = sum(patient_edge_supp) / len(patient_edge_supp)
            edges_supp.append(avg_supp_patient)

    avg_supp_total = sum(edges_supp) / len(edges_supp)

    return patient_probs, len(patient_probs), avg_supp_total


