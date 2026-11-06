# ---------------------------------------------------------------------------- #
#                                    IMPORTS                                   #
# ---------------------------------------------------------------------------- #
import networkx as nx
import gurobipy as gp
from gurobipy import GRB

# ---------------------------------------------------------------------------- #
#                                   FUNCTIONS                                  #
# ---------------------------------------------------------------------------- #


def get_max_support_graphs(node_set: frozenset, union_conflict_graph: nx.MultiGraph, input_graphs: dict):
    """
    For a fixed trajectory node set, select the maximum number of input graphs that contain all trajectory nodes and
    have no conflict among the trajectory nodes. At most one graph per patient, as graphs of the same patient are not
    compared in the conflict graph.
    """
    s = gp.Model('support')
    s.setParam('OutputFlag', 0)

    graph_vars = {g.name: s.addVar(vtype=GRB.BINARY, name=g.name)
                  for patient in input_graphs for g in input_graphs[patient] if node_set.issubset(g.nodes)}

    for patient in input_graphs:
        patient_vars = [graph_vars[g.name] for g in input_graphs[patient] if g.name in graph_vars]
        if len(patient_vars) > 1:
            s.addConstr(gp.quicksum(patient_vars) <= 1, 'Select at most one graph per patient')

    # trajectory nodes are fixed, so every conflict edge among them forbids selecting both of its graphs
    for _, _, label in union_conflict_graph.subgraph(node_set).edges(data='label'):
        g1, g2 = label.split(':')
        if g1 in graph_vars and g2 in graph_vars:
            s.addConstr(graph_vars[g1] + graph_vars[g2] <= 1, 'Conflict between graphs ' + g1 + ' ' + g2)

    s.setObjective(gp.quicksum(graph_vars.values()), GRB.MAXIMIZE)
    s.optimize()

    return [g for patient in input_graphs for g in input_graphs[patient]
            if g.name in graph_vars and graph_vars[g.name].X > 0.5]


def enumerate_max_trajectories(model, where):
    """
    Gurobi callback: record each new trajectory node set and exclude it with a lazy constraint on the node variables,
    which excludes all graph combinations leading to the same node set.
    """
    if where != GRB.Callback.MIPSOL:
        return
    values = model.cbGetSolution(model._node_vars)
    node_set = frozenset(n for n, v in zip(model._node_names, values) if v > 0.5)
    if node_set not in model._found:
        model._found.append(node_set)
    model.cbLazy(gp.quicksum(model._nodes[n] for n in node_set) <= len(node_set) - 1)
    # 0 searches for all trajectories until no feasible solution is left
    if model._max_trajectories and len(model._found) >= model._max_trajectories:
        model.terminate()


def find_max_k_common_trajectory(union_conflict_graph: nx.MultiGraph, input_graphs: dict, k: int, cores: int, solution_pool_size: int=1, verbose: bool=False):
    m = gp.Model('POTTR')
    if cores > 0:
        m.setParam('Threads', cores)
    m.setParam('OutputFlag', verbose)
    m.setParam('Seed', 42)

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
        g1, g2 = graph_names.split(':')
        p1, _ = g1.split('-')
        p2, _ = g2.split('-')
        a, b = node_names.split('->-')
        # if both graphs of a conflict edge are selected, edge must be active
        # m.addConstr(edges[edge] >= graphs[p1][g1] + graphs[p2][g2] - 1, 'Activate edge ' + edge + ' for graphs ' + g1 + ' ' + g2)
        # independent set constraint
        # m.addConstr(nodes[a] + nodes[b] <= 2 - edges[edge], 'Forbidden to include both nodes ' + a + ' ' + b)
        # combined to one constraint:
        m.addConstr(nodes[a] + nodes[b] <= 2 - (graphs[p1][g1] + graphs[p2][g2] - 1),
                    'Forbidden to include both nodes ' + a + ' ' + b + ' if graphs ' + g1 + ' ' + g2 + ' are selected')

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

    node_selection = []
    graph_selection = []
    if m.Status != GRB.OPTIMAL:
        return node_selection, graph_selection

    max_size = round(m.ObjVal)
    print('Size of a maximum trajectory for this instance is ', max_size - 1) # report size without root node
    first_node_set = frozenset(n for n in nodes if nodes[n].X > 0.5)
    found = [first_node_set]

    if solution_pool_size == 0 or solution_pool_size > 1:
        # search further trajectories of maximum size, each node set is reported only once
        m.addConstr(gp.quicksum(nodes) == max_size, 'Trajectory of maximum size')
        m.addConstr(gp.quicksum(nodes[n] for n in first_node_set) <= max_size - 1, 'Exclude first trajectory')
        m.setParam('LazyConstraints', 1)
        m._nodes = nodes
        m._node_names, m._node_vars = list(nodes.keys()), list(nodes.values())
        m._found = found
        m._max_trajectories = solution_pool_size
        # ends infeasible once all trajectories of maximum size are excluded
        m.optimize(enumerate_max_trajectories)
        if m.Status == GRB.INFEASIBLE:
            print('No further trajectories of maximum size exist (the infeasible model above is expected)')

    print('Number of trajectories found:', len(found), '| maximum trajectory size (objective value):', max_size - 1) # report size without root node

    for node_set in found:
        node_selection.append(sorted(node_set))
        graph_selection.append(get_max_support_graphs(node_set, union_conflict_graph, input_graphs))

    return node_selection, graph_selection