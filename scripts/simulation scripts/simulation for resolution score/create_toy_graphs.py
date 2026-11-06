import os
from itertools import combinations
from string import ascii_uppercase


def chain(labels):
    return [a + '->-' + b for a, b in zip(labels, labels[1:])]


def cluster(labels):
    # all pairs of the cluster are written, as in POTTR's own output
    return [a + '-?-' + b for a, b in combinations(labels, 2)]


def create_toy_graphs(labels, resolved_edges, reference_graphs, max_graphs):
    """
    The default graphs have all alterations in one cluster, e.g. {A, B, C, D, E}.
    reference_graphs graphs resolve resolved_edges edges: the first len(labels) - resolved_edges alterations stay in a
    cluster, followed by a chain of resolved_edges edges, e.g. {A, B, C} -> D -> E for resolved_edges = 2.
    """
    default_graph = ' '.join(cluster(labels))

    clustered, chained = labels[:len(labels) - resolved_edges], labels[len(labels) - resolved_edges:]
    relations = cluster(clustered)
    relations += [a + '->-' + chained[0] for a in clustered]
    relations += chain(chained)
    reference_graph = ' '.join(relations)

    graphs = []
    for i in range(reference_graphs):
        graphs.append(reference_graph)

    for i in range(max_graphs - reference_graphs):
        graphs.append(default_graph)

    return graphs


n = 10                  # graphs (patients) per file
# alterations per graph; a chain of n_alterations has n_alterations - 1 edges. Keep it small: the significance test
# enumerates all n_alterations! label permutations of the trajectory (compute_num_automorph)
n_alterations = 5
out_path = 'toy_graphs_resolved/'
labels = list(ascii_uppercase[:n_alterations])

# toy_graphs_resolved/<r>_resolved_edges/<i>_reference_graphs.txt:
# i of the n graphs are reference graphs that resolve r edges, the others have all alterations in one cluster
for resolved_edges in range(1, n_alterations):
    sub_path = os.path.join(out_path, str(resolved_edges) + '_resolved_edges')
    os.makedirs(sub_path, exist_ok=True)
    for i in range(1, n + 1):
        with open(os.path.join(sub_path, str(i) + '_reference_graphs.txt'), 'w') as f:
            for graph in create_toy_graphs(labels, resolved_edges, reference_graphs=i, max_graphs=n):
                f.write(graph + '\n')
