import os
import sys
import numpy as np
import math
import argparse
import networkx as nx
import matplotlib.pyplot as plt
import itertools
import random
from tqdm import tqdm

CODE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if CODE_DIR not in sys.path:
    sys.path.insert(0, CODE_DIR)

from read_input_dags import read_multiple_graphs_per_evolution


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", help="input file with graphs")
    parser.add_argument("-o", help="prefix of output file paths")
    parser.add_argument("-c", type=int, help="number of cores", default=1)
    parser.add_argument("-m", type=int, help="number of resamples to create")
    return parser


def permute_graph_permutation(graph_list, rng):
    debug_permute_graph = 0
    # use first graph of the list to create a new labeling
    # use the same relabelling for all graphs of the patient
    old_g_nodes = set([n for graph in graph_list for n in graph])
    old_g_nodes.remove('0')
    node_names = sorted([node_ for node_ in old_g_nodes])
    shuffled_names = list(node_names)
    rng.shuffle(shuffled_names)
    mapping = dict(zip(node_names, shuffled_names))

    relabeled_graphs = list()
    for graph in graph_list:
        # shuffle the labels of the relabeled graph, map labels back to gene names
        relabeled_g = nx.relabel_nodes(graph,mapping)
        for _, data in relabeled_g.nodes(data=True):
            if data.get('cluster_nodes'):
                data['cluster_nodes'] = {mapping.get(m, m) for m in data['cluster_nodes']}
        relabeled_graphs.append(relabeled_g)

    return relabeled_graphs


def main():
    rng = np.random.default_rng(42)

    args = get_parser().parse_args()
    input_graphs = args.g

    input_graphs_dict = read_multiple_graphs_per_evolution(path=input_graphs, out=args.o,
                                                           parallel_processes=args.c, verbose_flag=False)

    print("creating",args.m,"resamples")
    for i in tqdm(range(int(args.m))):
        new_folder = os.path.join(args.o, 'resample_' + str(i))
        os.makedirs(new_folder, exist_ok=True)
        num_patients_ = len(input_graphs_dict)
        for patient, graph_list in input_graphs_dict.items():
            # for p == 1 we perform the permutation test for POTTRs test
            graphs_perm = permute_graph_permutation(graph_list, rng)
            for index in range(len(graph_list)):
                g = graphs_perm[index]
                for _, data in g.nodes(data=True):
                    if 'cluster_nodes' in data:
                        data['cluster_nodes'] = ','.join(sorted(data['cluster_nodes']))
                nx.write_gexf(g, os.path.join(new_folder, patient + "-" + str(index) + "_resampled.gexf"))


if __name__ == '__main__':
  main()