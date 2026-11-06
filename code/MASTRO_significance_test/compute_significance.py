import glob
import os
import re

import numpy as np
import math
import networkx as nx
import itertools
import POTTR_compute_significance
from filelock import FileLock
from read_input_dags import read_multiple_graphs_per_evolution


def print_graph(tree):
    nx.draw(tree,with_labels=True)


def write_min_pvalue(file, min_pvalue):
    with FileLock(file + ".lock"):
        fout_ = open(file, "a")
        fout_.write(str(min_pvalue) + "\n")
        fout_.close()


def load_graph(line, name):
    edges_ = line.split(" ")
    G = nx.DiGraph()
    G.name = str(name)
    for edge in edges_:
        if "->-" in edge:
            nodes = edge.split("->-")
            G.add_edge(nodes[0], nodes[1])
        if "-/-" in edge:
            nodes = edge.split("-/-")
            G.add_node(nodes[0])
            G.add_node(nodes[1])
        if "-?-" in edge:
            nodes = edge.split("-?-")
            a = nodes[0]
            b = nodes[1]
            #G.add_edge(nodes[0] , nodes[1], label="?")
            #G.add_edge(nodes[1] , nodes[0], label="?")
            # modifications by Sara
            if a in G and 'cluster_nodes' in G.nodes[a]:
                if b not in G.nodes[a]['cluster_nodes']:
                    G.nodes[a]['cluster_nodes'].append(b)
            else:
                G.add_node(a, cluster_nodes=[b])
            if b in G and 'cluster_nodes' in G.nodes[b]:
                if a not in G.nodes[b]['cluster_nodes']:
                    G.nodes[b]['cluster_nodes'].append(a)
            else:
                G.add_node(b, cluster_nodes=[a])
    # Sara: changed "g" to "0"
    G.add_node("0")
    for node in G.nodes:
        if node != "0":
            G.add_edge("0" , node)
    G = nx.transitive_closure_dag(G)
    return G

def compute_support(traj, graphs_list):
    support = 0
    trans_id = 0
    trans_ids = []
    debug_support = 0
    traj = nx.transitive_closure_dag(traj)
    if debug_support == 1:
        print("traj.edges",traj.edges)
    for graph_ in graphs_list:
        sub_g = graph_.subgraph(traj.nodes).copy()
        if sub_g.edges == traj.edges:#nx.is_isomorphic(sub_g, traj):
            support += 1
            trans_ids.append(trans_id)
            if debug_support == 1:
                print("trans_id",trans_id,"sub_g.edges",sub_g.edges,"traj.edges",traj.edges)
                print("graph_.edges",graph_.edges)
        trans_id += 1
    return support, trans_ids

def findsubsets(set, subset_size):
    return list(itertools.combinations(set, subset_size))

def chernoffbound(p_mean , n , freq_traj):
    if freq_traj <= p_mean:
        return 1.
    q_ = 1 - p_mean
    if freq_traj == 0.:
        freq_traj += p_mean/1000.
    if freq_traj == 1:
        freq_traj -= p_mean/1000.
    pval_chern = ( (p_mean/freq_traj)**(freq_traj) * (q_/(1-freq_traj))**(1-freq_traj) )**n
    return pval_chern

def compute_pvalue_fast(probs , supp_traj , debug_pval=0):
    # consider only non-zero entries
    probs_ = []
    for p_ in probs:
        if p_ > 0:
            probs_.append(p_)
    probs = probs_
    n_traj = len(probs)
    if n_traj < supp_traj:
        print("Problem! n_traj < supp_traj , n_traj",n_traj,"supp_traj",supp_traj)
    indexes = [i for i in range(len(probs))]
    if debug_pval == 1:
        print("indexes",indexes)
        print("probs",probs)
        print("supp_traj",supp_traj)
    pval = 0.
    probs_dyn = np.zeros((n_traj+1,n_traj+1))
    probs_dyn[1,1] = probs_[0]
    probs_dyn[1,0] = (1-probs_[0])
    # fill the entire matrix
    for i in range(2,n_traj+1):
        for j in range(i+1):
            pred_idx = j-1
            if pred_idx >= 0:
                pred_prob = probs_dyn[i-1 , pred_idx]
            else:
                pred_prob = 0.
            probs_dyn[i,j] = probs_[i-1]*pred_prob + (1-probs_[i-1])*probs_dyn[i-1 , j]
    for j in range(supp_traj,n_traj+1):
        pval += probs_dyn[n_traj,j]
    return pval

def compute_pvalue_exact(probs , supp_traj , debug_pval=0):
    # consider only non-zero entries
    probs_ = []
    for p_ in probs:
        if p_ > 0:
            probs_.append(p_)
    probs = probs_
    n_traj = len(probs)
    if n_traj < supp_traj:
        print("ATTENTION! n_traj < supp_traj , n_traj",n_traj,"supp_traj",supp_traj)
        print('This might be due to cluster resolution')
    indexes = [i for i in range(len(probs))]
    if debug_pval == 1:
        print("indexes",indexes)
        print("probs",probs)
        print("supp_traj",supp_traj)
    pval = 0.
    for k in range(supp_traj,n_traj+1):
        subsets = findsubsets(indexes, k)
        for subset in subsets:
            if debug_pval == 1:
                print("considering subset",subset)
            p_this_subset = 1.
            indices_set = set(indexes)
            for idx in subset:
                p_this_subset = p_this_subset * probs[idx]
                indices_set.remove(idx)
            if debug_pval == 1:
                print("other indices",indices_set)
            for idx in indices_set:
                p_this_subset = p_this_subset * (1.-probs[idx])
            if debug_pval == 1:
                print("done, p_this_subset",p_this_subset)
            pval += p_this_subset
            if debug_pval == 1:
                print("pval increased to",pval)
    if debug_pval == 1:
        print("computed pval",pval)
    return pval


def compute_num_automorph(trajectory):
    set_t_nodes_nogerm = set(trajectory.nodes)
    # Sara: change "g" to "0"
    set_t_nodes_nogerm.remove("0")
    # compute number of automorphism of the trajectory
    automorph = 0
    trajectory_original = trajectory.copy()
    trajectory_permuted = trajectory.copy()
    set_t_nodes_nogerm_list = list(set_t_nodes_nogerm)
    #print("set_t_nodes_nogerm_list",set_t_nodes_nogerm_list)
    # Sara: this becomes too big for large trajectories; e.g. not feasible for 13 nodes
    list_permutations = list(itertools.permutations(set_t_nodes_nogerm_list))
    #print("list_permutations",list_permutations)
    num_t_nodes = len(set_t_nodes_nogerm_list)
    for permutation_ in list_permutations:
        # create a mapping
        map_p = dict()
        for i in range(len(permutation_)):
            map_p[set_t_nodes_nogerm_list[i]] = permutation_[i]
        #print("map_p",map_p,"perm",permutation_)
        trajectory_permuted = nx.relabel_nodes(trajectory_original, map_p)
        if trajectory_permuted.edges == trajectory_original.edges:
            automorph += 1
            if automorph > 1:
                #print("found ",automorph," automorphism!",map_p)
                #print_graph(trajectory_original)
                pass
    return automorph

def compute_prob(trajectory , graph_ , automorph_traj):
    t = len(list(graph_.nodes))-1
    k = len(list(trajectory.nodes))-1
    set_g_nodes_nogerm = set(graph_.nodes)
    # Sara: change "g" to "0"
    set_g_nodes_nogerm.remove("0")

    # compute number of subsets of nodes of graph isomorphic to trajectory
    subsets = findsubsets(set_g_nodes_nogerm, k)
    num_isomorph = 0

    debug_prob = 0
    for u, v, label_ in graph_.edges(data="label"):
        if label_ is not None:
            #debug_prob = 1
            debug_prob = 0
    if debug_prob == 1:
        print(graph_.edges(data="label"))
        print_graph(graph_)
        print_graph(trajectory)
    if debug_prob == 1:
        print("subsets",subsets)
    for subset_g in subsets:
        subset_g = list(subset_g)
        # Sara: change "g" to "0"
        subset_g.append("0")
        sub_g = graph_.subgraph(subset_g).copy()
        if nx.is_isomorphic(sub_g, trajectory):
            num_isomorph += 1
            if debug_prob == 1:
                print("found iso!")
        else:
            if debug_prob == 1:
                print("not iso...")
        if debug_prob == 1:
            #print_graph(trajectory)
            print_graph(sub_g)

    prob_perm = num_isomorph*automorph_traj
    for i in range(k):
        if t-i > 0:
            prob_perm = prob_perm/float(t-i)
    prob_indip = num_isomorph*automorph_traj/t**k

    #print(prob_)
    #print_graph(trajectory)
    return prob_indip , prob_perm

def compute_statistics(probs , supp_traj , n):
    var_ = 0.
    avg_ = 0.
    freq_traj = 0.
    t_stat = 0.
    for prob_ in probs:
        var_ += prob_*(1-prob_)
    var_ = var_/n**2
    if var_ == 0.:
        var_ = 1/n**6
    if var_ < 0:
        print("Error! var_",var_,"probs",probs,"supp_traj",supp_traj)
        exit()
    avg_ = np.array(probs).sum()/n
    freq_traj = supp_traj/n
    t_stat = freq_traj - avg_
    t_stat_norm = t_stat/math.sqrt(var_)
    pval_e = 0. #compute_pvalue_exact(probs , supp_traj)
    pval = compute_pvalue_fast(probs , supp_traj)
    if pval_e > 0. and abs(pval_e - pval) > pval_e/100:
        print("error in pvalue computation. pval_e",pval_e,"pval_f",pval)
        exit()
    if pval == 0.:
        compute_pvalue_exact(probs , supp_traj , 1)
        exit()

    non_zero_probs = []
    for p_ in probs:
        if p_ > 0.:
            non_zero_probs.append(p_)
    #pval_chern = chernoffbound(avg_ , n , freq_traj)
    #pval_chern = chernoffbound(np.array(non_zero_probs).mean() , len(non_zero_probs) , supp_traj/len(non_zero_probs))
    #print("pval",pval,"pval_chern",pval_chern,"non_zero_probs",non_zero_probs)
    return freq_traj , t_stat , t_stat_norm , avg_ , var_ , pval


def run_stat_significancce_test(support_file: str, graphs_dict: dict, output_file: str, minp: str=''):
    verbose = 0
    reading_ = True
    seps = ["->-", "-/-", "-?-"]

    # Updated by Sara
    # out_dir will actually not be used
    # out_dir = output_file.replace('/significance_output.txt', '')
    # graphs_list = [graph for id in graphs_dict for graph in graphs_dict[id]]

    # Sara added range and readlines here
    with open(support_file,"r") as file:
        fin_results = file.readlines()
        trajectories_list = []
        for i in range(0, len(fin_results), 2):
            line = fin_results[i]
            line = line.replace("\n","")
            if "(" in line:
                index_p = line.find("(")
                line_edges_only = line[:index_p]
                traj_rep_supp = line[index_p:]
                traj_rep_supp = traj_rep_supp.replace("(","")
                traj_rep_supp = int(traj_rep_supp.replace(")",""))
                # Sara changed this to reported support from own support computation
                next_line = fin_results[i+1].strip()
                trans_ids = next_line.split(' ')
                trajectory = load_graph(line_edges_only, '_'.join(trans_ids))
                trajectories_list.append((trajectory, traj_rep_supp, trans_ids))

        if verbose == 1:
            print("loaded",len(trajectories_list),"trajectories")

    min_pval = 1.
    fout_ = open(output_file,"w")
    fout_.write("edges_traj;traj_occ_list;traj_alt_occ_list;traj_supp;traj_freq;traj_exp_ind;traj_var_ind;t_stat_ind;t_stat_norm_ind;pval_ind;traj_exp_perm;traj_var_perm;t_stat_perm;t_stat_norm_perm;pval_perm;traj_POTTR;traj_var_POTTR;t_stat_POTTR;t_stat_norm_POTTR;pval_POTTR;resolution_evidence;support_score;\n")
    for (trajectory , supp_traj , trans_ids) in trajectories_list:
        id_index_map = {}
        for patient in graphs_dict:
            for i, graph_ in enumerate(graphs_dict[patient]):
                if set(trajectory.nodes).issubset(set(graph_.nodes)):
                    # Changes made by Sara
                    if patient in id_index_map:
                        id_index_map[patient].append(i)
                    else:
                        id_index_map[patient] = [i]

        if verbose == 1:
            print(supp_traj,len(id_index_map.keys()),trajectory.edges)

        probs_ind = []
        probs_perm = []
        automorph_traj = compute_num_automorph(trajectory)
        # Modifications by Sara to account for new orders and multiple trees per patient
        for patient_id, indices in id_index_map.items():
            patient_probs_ind = []
            patient_probs_perm = []
            for index in indices:
                graph = graphs_dict[patient_id][index]
                if str(graph.name) in trans_ids and (not set(trajectory.edges).issubset(set(graph.edges))):
                    ordered_graph = graph.copy()
                    for edge in trajectory.edges:
                        if edge not in ordered_graph.edges:
                            ordered_graph.add_edge(edge[0], edge[1])
                            ordered_graph = nx.transitive_closure_dag(ordered_graph)
                    prob_graph_indip, prob_graph_perm = compute_prob(trajectory, ordered_graph, automorph_traj)
                else:
                    prob_graph_indip , prob_graph_perm = compute_prob(trajectory, graph, automorph_traj)
                patient_probs_ind.append(prob_graph_indip)
                patient_probs_perm.append(prob_graph_perm)
            # Sara: Average over probablities computed per patient to not treat patient trees independently
            patient_avg_probs_ind = sum(patient_probs_ind) / len(patient_probs_ind)
            patient_avg_probs_perm = sum(patient_probs_perm) / len(patient_probs_perm)
            probs_ind.append(patient_avg_probs_ind)
            probs_perm.append(patient_avg_probs_perm)

        n = len(graphs_dict.keys())
        # compute all pvalues
        freq_traj , t_stat_ind , t_stat_ind_norm , avg_ind , var_ind , pval_indip = compute_statistics(probs_ind , supp_traj , n)
        freq_traj , t_stat_perm , t_stat_perm_norm , avg_perm , var_perm , pval_perm = compute_statistics(probs_perm , supp_traj , n)

        # Changes made by Sara
        # add POTTR statistical test, based on MASTRO's test, returns probability for current trajectory
        per_patient_probs, n_POTTR, candidate_graphs = POTTR_compute_significance.compute_stat_significance_POTTR(trajectory, graphs_dict)
        freq_traj_POTTR, t_stat_POTTR, t_stat_norm_POTTR, avg_POTTR, var_POTTR, pval_POTTR = compute_statistics(per_patient_probs, supp_traj, n_POTTR)
        #print('Freq MASTRO', freq_traj, 'Freq POTTR', freq_traj_POTTR, 'Supp mastro', supp_traj, 'supp pottr', supp_traj_POTTR)
        #print('n MASTRO', n, 'n POTTR', n)
        # compute the probabilities per patient that a random order of the patient's clusters fits the trajectory
        per_patient_support_probs, true_support, edge_support = POTTR_compute_significance.compute_support_probs(trajectory, graphs_dict)
        # evidence score is 1 if all patients fit the trajectory and smaller the more orders rely on resolved clusters
        evidence_support = sum(per_patient_support_probs) / true_support

        if verbose == 1:
            print("obs freq",freq_traj, "exp" , avg_ind ,"var", var_ind ,"t-stat", t_stat_ind ,"t-stat-norm", t_stat_ind_norm,"pval", pval_indip )

        min_pval_cand = pval_POTTR
        if min_pval_cand < min_pval:
            min_pval = min(min_pval,min_pval_cand)

        str_to_output = "["
        edge_id = 0
        for edge in trajectory.edges:
            if edge_id > 0:
                str_to_output = str_to_output+","
            str_to_output = str_to_output+str(edge[0])+"->-"+str(edge[1])
            edge_id += 1
        str_to_output = str_to_output+"];"
        str_to_output = str_to_output+str(trans_ids)+";"
        str_to_output = str_to_output+str(candidate_graphs)+";"
        str_to_output = str_to_output+str(supp_traj)+";"
        str_to_output = str_to_output+str(freq_traj)+";"
        # independent test
        str_to_output = str_to_output+str(avg_ind)+";"
        str_to_output = str_to_output+str(var_ind)+";"
        str_to_output = str_to_output+str(t_stat_ind)+";"
        str_to_output = str_to_output+str(t_stat_ind_norm)+";"
        str_to_output = str_to_output+str(pval_indip)+";"
        # permutation test
        str_to_output = str_to_output+str(avg_perm)+";"
        str_to_output = str_to_output+str(var_perm)+";"
        str_to_output = str_to_output+str(t_stat_perm)+";"
        str_to_output = str_to_output+str(t_stat_perm_norm)+";"
        str_to_output = str_to_output+str(pval_perm)+";"
        # POTTR test
        str_to_output = str_to_output + str(avg_POTTR) + ";"
        str_to_output = str_to_output + str(var_POTTR) + ";"
        str_to_output = str_to_output + str(t_stat_POTTR) + ";"
        str_to_output = str_to_output + str(t_stat_norm_POTTR) + ";"
        str_to_output = str_to_output + str(pval_POTTR) + ";"
        str_to_output = str_to_output + str(evidence_support) + ";"
        str_to_output = str_to_output + str(edge_support) + ";"
        str_to_output = str_to_output+"\n"
        fout_.write(str_to_output)

    fout_.close()

    if verbose:
        print("min pvalue",min_pval)

    if minp:
        write_min_pvalue(minp, min_pval)

