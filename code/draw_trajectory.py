import networkx as nx
import pandas as pd


def process_relations(edge_list):
    clusters = {}
    for edge in edge_list:
        if edge.__contains__('-?-'):
            a, b = edge.split('-?-')
            if a not in clusters and b not in clusters:
                cl = {a, b}
                clusters[a] = cl
                clusters[b] = cl
            elif a not in clusters:
                clusters[b].add(a)
                clusters[a] = clusters[b]
            elif b not in clusters:
                clusters[a].add(b)
                clusters[b] = clusters[a]
            else:
                clusters[a].update(clusters[b])
                for node in clusters[b]:
                    clusters[node] = clusters[a]

    nodes = set()
    for cl in clusters:
        label = ', '.join(sorted(clusters[cl]))
        nodes.add(label)

    edges = []
    for edge in edge_list:
        if edge.__contains__('->-'):
            a, b = edge.split('->-')
            label_a = a
            label_b = b
            if a in clusters:
                label_a = ', '.join(sorted(clusters[a]))
            if b in clusters:
                label_b = ', '.join(sorted(clusters[b]))
            nodes.add(label_a)
            nodes.add(label_b)
            edges.append((label_a, label_b))

    nodes = list(nodes)
    return nodes, edges


def draw_trajectory_graph(file_name: str, directory: str):
    with open(file_name, 'r') as file:
        lines = list(enumerate(file))

    file_name_number = 0
    input_graph_list = []
    names = {'file_number': [], 'graph_names': []}
    for i in range(0, len(lines), 2):
        trajectory, support = lines[i][1].split('(')
        trajectory_edges = trajectory.strip().split(' ')

        nodes, edges = process_relations(trajectory_edges)
        G = nx.DiGraph()
        G.add_nodes_from(nodes)
        G.add_edges_from(edges)

        G.add_node('0')
        for node in G.nodes:
            if node == '0':
                continue
            G.add_edge('0', node)

        trajectory_graphs = lines[i + 1][1].strip().split(' ')
        trajectory_graphs_concat = '_'.join(trajectory_graphs)
        input_graph_list.append(trajectory_graphs)

        graph = nx.transitive_reduction(G)
        dot = nx.drawing.nx_pydot.to_pydot(graph)
        dot.write_png(directory + 'trajectory_' + str(file_name_number) + '.png')
        names['file_number'].append(str(file_name_number))
        names['graph_names'].append(str(trajectory_graphs_concat))
        file_name_number += 1

    df = pd.DataFrame.from_dict(names)
    df.to_csv(directory + 'file_name_graph_matching.csv')
    return input_graph_list


def draw_input_graph(file_name: str, out_file: str):
    with open(file_name, 'r') as file:
        lines = list(enumerate(file))

    edges = lines[0][1].strip().split(' ')
    nodes, edges = process_relations(edges)
    G = nx.DiGraph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)

    G.add_node('0')
    for node in G.nodes:
        if node == '0':
            continue
        G.add_edge('0', node)

    graph = nx.transitive_reduction(G)
    dot = nx.drawing.nx_pydot.to_pydot(graph)
    dot.write_png(out_file)


if __name__ == '__main__':
    label_type = 'hugo+type_drivers_only'
    type = 'Hugo+Type'
    path = '/Users/sara/brandt/commonphylogenies/Output/new_selection_tracerx_' + label_type + '/tracerx_k3_pool50000/'
    input_graphs = draw_trajectory_graph(path + 'converted_graphs.txt', path)

    input_graph_path = '/Users/sara/brandt/commonphylogenies/Data/data_tracerx/all_trees_' + label_type + '/'
    for graphs in input_graphs:
        for g in graphs:
            draw_input_graph(file_name=input_graph_path + g + '_tracerx_tree_' + type + '.txt', out_file=path + g + '.png')

    #draw_input_graph('/Users/sara/brandt/commonphylogenies/Data/data_tracerx/all_trees_hugo_drivers_only/CRUK0753-0_tracerx_tree_Hugo.txt', 'test-0.png')
    #draw_input_graph('/Users/sara/brandt/commonphylogenies/Data/data_tracerx/all_trees_hugo_drivers_only/CRUK0753-1_tracerx_tree_Hugo.txt', 'test-1.png')

