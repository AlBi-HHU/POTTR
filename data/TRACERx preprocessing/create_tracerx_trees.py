# ---------------------------------------------------------------------------- #
#                                    IMPORTS                                   #
# ---------------------------------------------------------------------------- #
import argparse
import glob
import itertools
import os
import re
import pandas as pd
import networkx as nx
from pyreadr import pyreadr


# ---------------------------------------------------------------------------- #
#                                   FUNCTIONS                                  #
# ---------------------------------------------------------------------------- #
def create_column(df):
    df.loc[(df['func'] == 'splicing') & (df['exonic.func'].isna()), 'exonic.func'] = 'M'
    df['Hugo+Type'] = df['Hugo_Symbol'] + '_' + df['exonic.func'].replace({'nonsynonymous SNV': 'M',
                                                                           'nonsynonymous': 'M',
                                                                           'nonframeshift substitution': 'M',
                                                                           'stopgain SNV': 'N',
                                                                           'immediate-stopgain': 'N',
                                                                           'nonframeshift insertion': 'M',
                                                                           'frameshift substitution': 'N',
                                                                           'frameshift insertion': 'N'})

    return df['Hugo+Type']


def get_clusters(table, patient_id, scheme):
    # remove rows where cluster assignment is nan
    mut_table = pd.DataFrame(table[table['PyCloneCluster_SC'].notna()])
    mut_table['PyCloneCluster_SC'] = mut_table['PyCloneCluster_SC'].astype(int)
    if scheme == 'Hugo+Cluster':
        return mut_table[mut_table['tumour_id'] == patient_id].groupby('PyCloneCluster_SC')['Hugo+Cluster'].apply(list).to_dict()
    elif scheme == 'Hugo+AAChange':
        return mut_table[mut_table['tumour_id'] == patient_id].groupby('PyCloneCluster_SC')['Hugo+AAChange'].apply(list).to_dict()
    elif scheme == 'Hugo':
        return mut_table[mut_table['tumour_id'] == patient_id].groupby('PyCloneCluster_SC')['Hugo_Symbol'].apply(list).to_dict()
    elif scheme == 'Hugo+Type':
        return mut_table[mut_table['tumour_id'] == patient_id].groupby('PyCloneCluster_SC')['Hugo+Type'].apply(list).to_dict()
    else:
        print('Unknown scheme {0}'.format(scheme))
        exit(-1)


def get_edges(delimiter, tuple_list):
    edge_list = []
    for a, b in tuple_list:
        a = str(a)
        b = str(b)
        if a == 'nan' or b == 'nan':
            continue
        elif a == b:
            continue
        else:
            edge_list.append(a + delimiter + b)
    return edge_list


def edges_to_text(G):
    edge_list = []
    for edge in G.edges:
        x, y = edge
        cluster_edges_x = get_edges('-?-', itertools.combinations(G.nodes[x]['mutations'], 2))
        edge_list.extend(cluster_edges_x)
        cluster_edges_y = get_edges('-?-', itertools.combinations(G.nodes[y]['mutations'], 2))
        edge_list.extend(cluster_edges_y)

        if x != y:
            directed_edges_list = list(itertools.product(G.nodes[x]['mutations'], G.nodes[y]['mutations']))
            edges = get_edges('->-', directed_edges_list)
            edge_list.extend(edges)
    return edge_list


def add_node_labels(graph, root, mut_table, patient_id, scheme):
    cluster_dict = get_clusters(mut_table, patient_id, scheme)

    queue = [root]
    all_mutations = set()
    visited = set()
    # some clusters might be empty now due to our filtering, we have to remove these nodes later
    nodes_to_remove = []
    # start breadth-first-search to add labels without duplicates
    while queue:
        # pop first element of queue indexed by 0
        node = queue.pop(0)
        if node in visited:
            continue
        visited.add(node)

        if node in cluster_dict:
            node_cluster = set(cluster_dict[node])
            # get intersectin of all mutations seen so far and mutations of current cluster
            shared = all_mutations & node_cluster
            if shared:
                print('Duplication found in cluster {0} of patient {1}'.format(node, patient_id))
                print(shared)
            cluster_no_dp = node_cluster - shared

            graph.nodes[node]['mutations'] = list(cluster_no_dp)
            all_mutations |= cluster_no_dp
        else:
            print('Node {0} not in cluster dict'.format(node))
            nodes_to_remove.append(node)

        for succ in graph.successors(node):
            if succ not in visited:
                queue.append(succ)

    tv_closed_graph = nx.transitive_closure_dag(graph)
    tv_closed_graph.remove_nodes_from(nodes_to_remove)
    return tv_closed_graph


def build_tree(df, mut_table, patient, scheme):
    tuple_list = [tuple(x) for x in df.itertuples(index=False, name=None)]

    G = nx.DiGraph()
    G.add_edges_from(tuple_list)
    G.remove_edges_from(nx.selfloop_edges(G))
    roots = [n for n in G.nodes if G.in_degree(n) == 0]
    inserted_root = False
    if len(roots) > 1:
        print('Multiple root nodes found')
        exit(-1)
    if len(roots) == 0:
        print('No root found in patient {0}'.format(patient))
        inserted_root = True
        roots = ['root']
        G.add_node('root')
        for node in G.nodes:
            if node == 'root':
                continue
            else:
                G.add_edge('root', node)

    labelled_graph = add_node_labels(G, roots[0], mut_table, patient, scheme)
    if inserted_root:
        labelled_graph.remove_node('root')

    return labelled_graph


def get_parser():
    """Get parser object for combinatorial vaccine design."""
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', '-i', required=True, dest='input', type=str,
                        help='Path to input files')
    parser.add_argument('--output', '-o', required=True, dest='output', type=str,
                        help='Path to store output file')
    parser.add_argument('--mutation-table', '-mut', required=True, dest='mutation_table', type=str,
                        help='Path to mutation table file')
    parser.add_argument('--scheme', '-s', required=True, dest='scheme', type=str,
                        help='Labelling scheme of nodes. Choose from Hugo (only Hugo symbol, first occurrence, no '
                             'duplications), Hugo+Cluster (Hugo symbol and cluster number), or Hugo+AAChange (Hugo '
                             'symbol and mutation level information), Hugo+Type (Hugo symbol and mutation type)')
    parser.add_argument('--only-drivers', '-drivers', dest='drivers', action='store_true',
                        help='If set, mutations will be filtered to only keep driver mutations. '
                             'Default, mutations will not be filtered.')
    parser.add_argument('--gene-list', '-g', required=False, dest='gene_list', type=str,
                        help='Genes of interest, provided as tsv')
    return parser


def main():
    args = get_parser().parse_args()

    # global mut_table
    mut_table_tmp = pyreadr.read_r(args.mutation_table)
    mut_table = pd.DataFrame(mut_table_tmp[None])
    # remove synonymous mutations
    mut_table = mut_table[(mut_table['exonic.func'] != 'synonymous') &
                          (mut_table['exonic.func'] != 'synonymous SNV')]
    if args.drivers:
        mut_table = mut_table[mut_table['DriverMut'] == True]
    mut_table['Hugo+Cluster'] = mut_table['Hugo_Symbol'] + '_' + mut_table['PyCloneCluster_SC'].astype(str)
    mut_table['Hugo+AAChange'] = mut_table['Hugo_Symbol'] + '_' + mut_table['AAChange']
    mut_table['Hugo+Type'] = create_column(mut_table)

    dir = args.output + '/'
    if not os.path.isdir(dir):
        os.mkdir(dir)

    if args.gene_list:
        if os.path.isfile(args.gene_list):
            if args.gene_list.endswith('.tsv'):
                print('Filter for genes of interest provided in file {0}'.format(args.gene_list))
                df = pd.read_csv(args.gene_list, sep='\t')
                df = df[(df['Is Oncogene'] == 'Yes') | (df['Is Tumor Suppressor Gene'] == 'Yes')]
                gene_list = df['Hugo Symbol'].dropna().tolist()
                mut_table = mut_table[mut_table['Hugo_Symbol'].isin(gene_list)]
            else:
                print('Gene list is not in tsv format. Exiting.')
                exit(-1)

    scheme = args.scheme
    for i, pf in enumerate(glob.glob(args.input + '/CRUK*.xlsx')):
        if i % 100 == 0:
            print('Step ', i)
        patient_file = pf.split('/')[-1]
        patient = re.search(r'(CRUK.*)\.xlsx', patient_file).group(1)
        # sheet_name = None creates a dict of dfs, needed to handle multiple sheets
        dfs = pd.read_excel(pf, sheet_name=None)
        tmp_lst = []
        for i, sheet in enumerate(dfs):
            df = dfs[sheet]
            tree = build_tree(df, mut_table, patient, scheme)
            if tree:
                # do not write duplicates if multiple trees are the same
                if set(tree.nodes) in tmp_lst and set(tree.edges) in tmp_lst:
                    continue
                else:
                    tmp_lst.append(set(tree.nodes))
                    tmp_lst.append(set(tree.edges))

                    line_to_write = []
                    # need to use nx.all_neighbors to get all adjacencies in directed graph
                    nodes_wo_edges = [n for n in tree.nodes if not list(nx.all_neighbors(tree, n))]
                    # must add nodes without any neighbors separately to output
                    for node in nodes_wo_edges:
                        if len(tree.nodes[node]['mutations']) == 1:
                            line_to_write.extend(tree.nodes[node]['mutations'])
                        # create "cluster edges" between mutations of the node
                        single_node_cluster = get_edges('-?-', itertools.combinations(tree.nodes[node]['mutations'], 2))
                        line_to_write.extend(single_node_cluster)

                    edges = edges_to_text(tree)
                    line_to_write.extend(edges)
                    if line_to_write:
                        with open(dir + '/' + patient + '-' + str(i) + '_tracerx_tree_' + scheme + '.txt', 'w') as file:
                            file.write(' '.join(line_to_write) + '\n')
                    else:
                        print('No line for patient ', patient)
            else:
                print('No tree for patient ', patient)


if __name__ == '__main__':
    main()
