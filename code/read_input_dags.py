# ---------------------------------------------------------------------------- #
#                                    IMPORTS                                   #
# ---------------------------------------------------------------------------- #
import glob
import os
import pickle
import re

import pandas as pd
import networkx as nx
from create_graphs import get_graphs_parallel, get_graphs_single_thread

# ---------------------------------------------------------------------------- #
#                               GLOBAL VARIABLES                               #
# ---------------------------------------------------------------------------- #
verbose: bool


# ---------------------------------------------------------------------------- #
#                                   FUNCTIONS                                  #
# ---------------------------------------------------------------------------- #
def log(string):
    if verbose:
        print(string)


def get_graphs_from_lines(lines, tree_numbers):
    log('Start parallel creation of graphs from input')
    graphs = {}
    for evol_process in lines:
        patient_graphs = get_graphs_parallel(lines[evol_process], tree_numbers[evol_process])
        graphs[evol_process] = patient_graphs

    return graphs


def file_name_match(file_name):
    # indicator for multiple trees per evolutionary process is '-'
    if '-' in file_name:
        match = re.search(r'((.*)-\d+)\D.*(txt|gexf)', file_name)
        evolution = match.group(2)
        phylo_tree = match.group(1)
    else:
        match = re.search(r'((.*)\d+)\D.*(txt|gexf)', file_name)
        evolution = match.group(1)
        phylo_tree = evolution + '-0'
    return evolution, phylo_tree


def read_multiple_graphs_per_evolution(path: str, out: str, verbose_flag: bool=False):
    global verbose
    verbose = verbose_flag

    if not os.path.isdir(path) and not os.path.isfile(path):
        print('Please provide an existing folder to input trees')
        exit(-1)

    log('Start iterating over multiple input trees for each evolution')
    # iterate over TRACERx input files TODO: generalize to other input data
    evolution_lines = {}
    phylo_tree_numbers = {}

    if path.endswith('.txt'):
        with open(path) as f:
            for i, line in enumerate(f):
                evolution_lines[str(i)] = [line]
                phylo_tree_numbers[str(i)] = [str(i) + '-0']
        graphs = get_graphs_from_lines(evolution_lines, phylo_tree_numbers)
    else:
        filelist = glob.glob(path + '/*.txt')
        for file in sorted(filelist):
            file_name = file.split('/')[-1]
            # requires tree id "-\d+"
            evolution, phylo_tree = file_name_match(file_name)
            with (open(file, 'r') as f):
                for line in f:
                    tmp = line.strip().split(' ')
                    sorted_line = ' '.join(sorted(tmp)) + '\n'
                    if evolution not in evolution_lines:
                        evolution_lines[evolution] = [sorted_line]
                        phylo_tree_numbers[evolution] = [phylo_tree]
                    else:
                        if sorted_line not in evolution_lines[evolution]:
                            evolution_lines[evolution].append(sorted_line)
                            phylo_tree_numbers[evolution].append(phylo_tree)
        graphs = get_graphs_from_lines(evolution_lines, phylo_tree_numbers)

        # allow for graphs in gexf format
        gexflist = glob.glob(path + '/*.gexf')
        for gexf_file in gexflist:
            file_name = gexf_file.split('/')[-1]
            evolution, phylo_tree = file_name_match(file_name)
            # read in graph from gexf file
            graph = nx.read_gexf(gexf_file)
            # add root node and connect it to all other nodes
            graph.add_node('0')
            graph.add_edges_from([('0', node) for node in graph.nodes if node != '0'])
            graph = nx.transitive_closure_dag(graph)
            graph.name = str(phylo_tree)
            if evolution in graphs:
                graphs[evolution].append(graph)
            else:
                graphs[evolution] = [graph]

    df = pd.DataFrame([(p, len(n)) for p, n in phylo_tree_numbers.items()], columns=['evolution', 'distinct trees'])
    df.to_csv(out + '/number_of_distinct_dags_per_sample.csv')

    log('Number of graphs read from input: ' + str(len(graphs)))
    return graphs