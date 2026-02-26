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


def read_multiple_graphs_per_evolution(path: str, out: str, parallel_processes: int, verbose_flag: bool=False):
    global verbose
    verbose = verbose_flag

    if not os.path.isdir(path) and not os.path.isfile(path):
        print('Please provide an existing folder to input trees')
        exit(-1)

    log('Start iterating over multiple input trees for each evolution')
    evol_processes = [] # triplet list with id, tree number and the read in line

    '''
    if provided a single txt file, each line is assumed to be a distinct evolutionary process, 
    i.e. distinct tumor / patient
    '''
    if path.endswith('.txt'):
        with open(path) as f:
            for i, line in enumerate(f):
                evol_processes.append((str(i), str(i) + '-0', line))
        graphs = get_graphs_parallel(evol_processes, parallel_processes)
    else:
        filelist = glob.glob(path + '/*.txt')
        for file in sorted(filelist):
            file_name = file.split('/')[-1]
            # requires tree id "-\d+"
            evolution, phylo_tree = file_name_match(file_name)

            # read in multiple trees per patient, but skip over duplicates
            tmp_duplicates = dict()
            with (open(file, 'r') as f):
                for line in f:
                    tmp = line.strip().split(' ')
                    sorted_line = ' '.join(sorted(tmp)) + '\n'
                    if sorted_line in tmp_duplicates:
                        continue
                    else:
                        evol_processes.append((evolution, phylo_tree, sorted_line))
                        tmp_duplicates[sorted_line] = None

        graphs = get_graphs_parallel(evol_processes, parallel_processes)

        # allow for graphs in gexf format
        gexflist = glob.glob(path + '/*.gexf')
        for gexf_file in gexflist:
            file_name = gexf_file.split('/')[-1]
            evolution, phylo_tree = file_name_match(file_name)
            # read in graph from gexf file
            graph = nx.read_gexf(gexf_file)

            for n, data in graph.nodes(data=True):
                if 'cluster_nodes' in data:
                    val = data.get('cluster_nodes')
                    data['cluster_nodes'] = set(val.split(',')) if val else set()

            # add root node and connect it to all other nodes
            graph.add_node('0')
            graph.add_edges_from([('0', node) for node in graph.nodes if node != '0'])
            graph = nx.transitive_closure_dag(graph)
            graph.name = str(phylo_tree)
            if evolution in graphs:
                graphs[evolution].append(graph)
            else:
                graphs[evolution] = [graph]

    if verbose:
        id_tree_pairs = [(e, len(trees)) for e, trees in graphs.items()]
        df = pd.DataFrame(id_tree_pairs, columns=['evolution', 'distinct trees'])
        df.to_csv(out + '/number_of_distinct_trees_per_patient.csv')
        num = sum([pair[1] for pair in id_tree_pairs])
        log('Number of graphs read from input: ' + str(num))

    return graphs