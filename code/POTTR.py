# ---------------------------------------------------------------------------- #
#                                    IMPORTS                                   #
# ---------------------------------------------------------------------------- #
import networkx as nx
import gurobipy as gp
from gurobipy import GRB

# ---------------------------------------------------------------------------- #
#                                   FUNCTIONS                                  #
# ---------------------------------------------------------------------------- #


def find_max_k_common_trajectory(union_conflict_graph: nx.MultiGraph, input_graphs: dict, k: int, path: str, cores: int, solution_pool_size: int=5000, verbose: bool=False):
    m = gp.Model('POTTR')
    if cores > 0:
        m.setParam('Threads', cores)
    m.setParam('OutputFlag', verbose)
    m.setParam('Seed', 42)

    if solution_pool_size > 0:
        m.setParam('PoolSearchMode', 2)  # Search for alternative optimal solutions
        m.setParam('PoolGap', 0.00001)
        m.setParam('PoolSolutions', solution_pool_size)

    # define variables
    nodes = gp.tupledict()
    for n in union_conflict_graph.nodes:
        nodes[n] = m.addVar(vtype=GRB.BINARY, name=n)

    edges = gp.tupledict()
    for e in union_conflict_graph.edges:
        a, b, _ = e
        edge_name = a + '->-' + b + ';' + union_conflict_graph.edges[e]['label']
        edges[edge_name] = m.addVar(vtype=GRB.BINARY, name=edge_name)

    graphs = dict()
    for patient in input_graphs:
        graphs[patient] = gp.tupledict()
        for g in input_graphs[patient]:
            name = g.name
            graphs[patient][name] = m.addVar(vtype=GRB.BINARY, name=name)

    # set constraints
    m.setObjective(gp.quicksum(nodes), GRB.MAXIMIZE)

    for edge in edges:
        node_names, graph_names = edge.split(';')
        # if both graphs of a conflict edge are selected, edge must be active
        g1, g2 = graph_names.split(':')
        p1, _ = g1.split('-')
        p2, _ = g2.split('-')
        m.addConstr(edges[edge] >= graphs[p1][g1] + graphs[p2][g2] - 1,
                    'Activate edge ' + edge + ' for graphs ' + g1 + ' ' + g2)
        # independent set constraint
        a, b = node_names.split('->-')
        m.addConstr(nodes[a] + nodes[b] <= 2 - edges[edge], 'Forbidden to include both nodes ' + a + ' ' + b)

    for patient in input_graphs:
        # in case of multiple graphs per patient, we ensure that at most one can be selected per patient
        m.addConstr(gp.quicksum(graphs[patient]) <= 1, 'Select at most one graph per patient')

        for g in input_graphs[patient]:
            # add constraint that selected node must occur in all selected graphs
            for node in union_conflict_graph.nodes:
                constant = int(node in g.nodes)
                m.addConstr(nodes[node] <= constant + (1 - graphs[patient][g.name]),
                            'Selected nodes must occur in all selected graphs')
    
    m.addConstr(sum([gp.quicksum(graphs[patient]) for patient in input_graphs]) >= k, 'Select k ' + str(k) + ' graphs')

    m.optimize()

    max_size = int(m.ObjVal)
    print('Size of a maximum trajectory for this instance is ', max_size)

    node_selection = []
    graph_selection = []
    if m.Status == GRB.OPTIMAL:
        for i in range(m.SolCount):
            m.setParam('SolutionNumber', i)
            if int(m.PoolObjVal) == max_size:
                solution_nodes = sorted([n for n in nodes if nodes[n].Xn > 0.5])
                solution_graphs = []
                for patient in input_graphs:
                    solution_graphs.extend([g for g in input_graphs[patient] if graphs[patient][g.name].Xn > 0.5])
                if solution_nodes not in node_selection:
                    node_selection.append(solution_nodes)
                    graph_selection.append(solution_graphs)

    return node_selection, graph_selection