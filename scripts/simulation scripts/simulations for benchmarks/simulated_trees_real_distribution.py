import argparse
import sys
import os.path
import networkx as nx
from random import normalvariate, choices, seed
from statistics import stdev
from collections import defaultdict

# make POTTR's modules in code/ importable regardless of the working directory
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'code'))
from read_input_dags import read_multiple_graphs_per_evolution


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--number', '-n', required=True, dest='num', type=int,
                        help='Number of graphs to generate')
    parser.add_argument('--input_data', '-i', required=True, dest='input', type=str,
                        help='Input graphs used as reference to estimate edge frequencies')
    parser.add_argument('--output', '-o', required=True, dest='output', type=str,
                        help='Path to store generated graphs')
    parser.add_argument('--threads', '-t', dest='threads', type=int, default=2,
                        help='Number of parallel threads (default 2)')
    parser.add_argument('--edge_mean', '-mu', dest='edge_mean', type=int, default=None,
                        help='Mean number of edges in graphs (default takes mean from reference data)')
    parser.add_argument('--edge_dev', '-sigma', dest='edge_dev', type=int, default=None,
                        help='Standard deviation of number of edges in graphs (default takes sigma from reference data)')
    parser.add_argument('--min_edges', '-min', dest='min_edges', type=int, default=1,
                        help='Minimum number of edges in a graph (default 1)')
    parser.add_argument('--max_edges', '-max', dest='max_edges', type=int, default=None,
                        help='Maximum number of edges in a graph (default takes maximum from reference data)')
    return parser


def main():
    args = get_parser().parse_args()
    reference_data = args.input
    output_data = args.output + '_' + str(args.num) + '/'
    os.makedirs(output_data, exist_ok=True)

    seed(42)

    graphs_to_generate = args.num
    edge_mean = args.edge_mean
    edge_dev = args.edge_dev
    min_edges = args.min_edges
    max_edges = args.max_edges

    # read in reference graphs to get an edge distribution from which we may sample
    graphs_dict = read_multiple_graphs_per_evolution(reference_data, output_data , args.threads, False)

    edges = defaultdict(lambda: 0)
    # for sampling we only use transitive reduced graphs and in POTTR we will add the transitive edges again
    # because we assume that transitive edges exist in any case, so we don't have to simulate them
    graph_count = 0
    edge_count = []
    for evol_id in graphs_dict:
        for graph in graphs_dict[evol_id]:
            red_graph = nx.transitive_reduction(graph)
            graph_count += 1
            edge_count.append(len(red_graph.edges))
            for e in red_graph.edges:
                edges[e] += 1

    if not edge_dev:
        edge_dev = stdev(edge_count)

    if not edge_mean:
        edge_mean = sum(edge_count) / graph_count

    if not min_edges:
        min_edges = min(edge_count)

    if not max_edges:
        max_edges = max(edge_count)

    print('Sampling with deviation %s, mean %s, minimum of %s edges, maximum of %s edges' % (edge_dev, edge_mean, min_edges, max_edges))

    # sample edges with frequency equal to in how many graphs the edge exists
    # only sample from edges not adjacent to two nodes in the graph
    edge_list = [e for e in edges]
    edge_weights = [edges[e] / graph_count for e in edge_list]

    sampled_graphs = []
    for i in range(graphs_to_generate):
        # sample num of edges from normal distribution with user defined or measured mean and stdev
        num_edges = int(normalvariate(edge_mean, edge_dev))
        num_edges = max(min_edges, num_edges)
        num_edges = min(max_edges, num_edges)

        graph = nx.DiGraph()
        while len(graph.edges) < num_edges:
            candidates = [index for index, (u, v) in enumerate(edge_list) if u not in graph.nodes or v not in graph.nodes]
            candidate_weights = [edge_weights[index] for index in candidates]
            if not candidates:
                print('Graph %d Candidate list is empty. Could only sample %d of %d edges.' % (i, len(graph.edges), num_edges))
                break
            # sample edge here
            sample = choices(candidates, weights=candidate_weights)
            e = edge_list[sample[0]]
            u, v = e
            graph.add_edge(u, v)

        graph.add_node('0')
        for n in graph.nodes():
            if n == '0':
                continue
            graph.add_edge('0', n)

        clusters = defaultdict(list)
        for n in graph.nodes():
            preds = frozenset(graph.predecessors(n))
            succs = frozenset(graph.successors(n))
            if preds and succs:
                clusters[(preds, succs)].append(n)

        for cluster in clusters:
            if len(clusters[cluster]) > 1:
                # print('Cluster detected')
                for node in clusters[cluster]:
                    tmp_set = set(clusters[cluster]) - {node}
                    nx.set_node_attributes(graph, {node: {'cluster_nodes': ','.join(sorted(tmp_set))}})

        nx.write_gexf(graph, output_data + str(i) + '_graph.gexf')

    return


if __name__ == '__main__':
    main()