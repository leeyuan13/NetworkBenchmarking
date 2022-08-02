import numpy as np
import matplotlib.pyplot as plt
from itertools import chain, combinations
from ortools.linear_solver import pywraplp
import scipy.optimize
import time
import seaborn as sns
import os

from sim7_util import get_ecal, gen_all_coalitions
from sim7_util_chain import gen_params_chain, gen_coalitions_chain
from sim7_util_grid import gen_params_grid, gen_coalitions_grid, \
                             gen_coalitions_grid_sq
from sim7_util_dumbbell import gen_params_dumbbell, gen_coalitions_dumbbell

# Sets:
#   nodes = list of nodes
#   ecal = list of possible entangled links, superset of elementary links

# Parameters:
#   p_edges [dict] = link: physical-layer rate
#   q_nodes [dict] = node: entanglement swap efficiency
#   vol_scaling [int] = base for coalition volume contribution
#                       (i.e. volume = vol_scaling**size)

# Variables:
#   y[s] = rate vector = amount of state s (in ecal) to cash out
#   w[s1][s2] = amount of state s1 (in ecal) used to produce state s2 (in ecal)

## Data folders:
# sim7_data - true data
# sim7_data_old - missing the e2c(d, n) factor

get_edge = lambda i, j: (i, j) if i < j else (j, i)

# Implicit parameter: error of generated entanglement.
eps_edge = 0

def run_LP_volume(nodes, ecal, p_edges, q_nodes, vol_scaling, \
                  gen_coalitions = gen_all_coalitions):
    get_y_index = {ecal[i]: i * len(ecal) + i for i in range(len(ecal))}
    get_w_index = {(ecal[i], ecal[j]): i * len(ecal) + j \
                   for i in range(len(ecal)) for j in range(len(ecal)) if j != i}
    get_rr_index_list = dict() # edge -> list of coalition indices
    coalitions = gen_coalitions(nodes)
    coal_index = len(ecal)**2
    for coal in coalitions:
        for i1 in range(len(coal)):
            for i2 in range(i1+1, len(coal)):
                get_rr_index_list.setdefault(get_edge(coal[i1], coal[i2]), []).\
                                                        append(coal_index)
        coal_index += 1

    num_y = len(ecal)
    num_w = len(ecal)**2 - len(ecal)
    num_coalitions = coal_index - len(ecal)**2
    num_vars = len(ecal)**2 + num_coalitions
    print('Number of variables:', num_vars)

    # In the case where the network only produces entanglement of fixed
    # fidelity, i.e. where any entangled link can only be produced with some
    # given error rate, the effective size of a coalition is a pre-determined
    # constant. Note that the error rate can be specific to the link;
    # the end user must just have no choice as to what the error rate can be.
    # Then the optimization problem can be written as a linear program.
    eps_coalition = eps_edge # from parent scope --> error rate for any coalition
    coalitions = gen_coalitions(nodes)  # recall this is a generator
                                        # --> one use only
    # Depth of computation associated with a coalition.
    # Note that we should always have depth <= number of memories in a
    # coalition, so the depth is also the exponent in the volume calculation.
    # When a coalition is too big, the maximum depth is limited by fidelity.
    m_sol = []
    # Number of nodes in coalition.
    coal_len = []
    # We want to decide the optimal depth for each coalition.
    # The relevant factor in the objective function is:
    obj_m = lambda d: d * np.log(vol_scaling) - np.log(d)
    for coal in coalitions:
        if eps_coalition == 0 or len(coal) < 1/np.sqrt(eps_coalition):
            # Pick the optimal depth of computation.
            if obj_m(1) <= obj_m(len(coal)):
                m_sol.append(len(coal))
            else:
                m_sol.append(1)
        else:
            if obj_m(1) <= obj_m(int(1/eps_coalition/len(coal))):
                m_sol.append(int(1/eps_coalition/len(coal)))
            else:
                m_sol.append(1)
        coal_len.append(len(coal))

    solver = pywraplp.Solver.CreateSolver('GLOP')
    infinity = solver.infinity()
    variables = [solver.NumVar(0, infinity, 'x['+str(i)+']') for i in range(num_vars)]

    # Note: we allow for exactly one computation memory per node,
    # but we allow for infinitely many communication memories per node.
    # Hence, rate region is given by the p-q model;
    # the number of memoriesin a cluster of nodes is equal to the number of nodes.

    # Reaction ratio constraints.
    for s in ecal:
        for c in nodes:
            if c == s[0] or c == s[1]:
                continue
            else:
                w1 = (get_edge(s[0], c), s)
                w2 = (get_edge(s[1], c), s)
                solver.Add(variables[get_w_index[w1]] == variables[get_w_index[w2]])

    # Balance.
    for s in ecal:
        terms_in = []
        terms_out = []
        a, b = s
        for c in nodes:
            if c == a or c == b:
                continue
            else:
                # Incoming.
                wACab = variables[get_w_index[(get_edge(a, c), get_edge(a, b))]]
                wBCab = variables[get_w_index[(get_edge(b, c), get_edge(a, b))]]
                wABac = variables[get_w_index[(get_edge(a, b), get_edge(a, c))]]
                wABbc = variables[get_w_index[(get_edge(a, b), get_edge(b, c))]]
                terms_in.append((wACab + wBCab)/2 * q_nodes[c])
                terms_out.append(wABac + wABbc)
        yab = variables[get_y_index[s]]
        pab = p_edges.get(s, 0) # default to zero if s is not elementary
        solver.Add(yab <= sum(terms_in) + pab - sum(terms_out))

    # Coalition rate requirements.
    for s in ecal:
        yab = variables[get_y_index[s]]
        rri = []
        involved_coalitions = get_rr_index_list[s]
        for coal_index in involved_coalitions:
            rri.append(variables[coal_index])
        solver.Add(sum(rri) <= yab)

    # Objective.
    # Convert coalition entanglement rate to coalition computation rate, including
    # the depth d and the size of the coalition n.
    e2c = lambda d, n: (n-1)/d if n % 2 == 0 else n/d
    scale_factor = 0.01 # rough fix to avoid overflow issues
    microvolumes = [scale_factor * (vol_scaling**m_sol[i]) * variables[len(ecal)**2+i] * e2c(m_sol[i], coal_len[i]) for i in range(num_coalitions)]
    solver.Maximize(sum(microvolumes))
    # Solve.
    print('Solving...')
    status = solver.Solve()
    if status != pywraplp.Solver.OPTIMAL:
        print(status, solver.Objective().Value())
        # Output codes: https://google.github.io/or-tools/python/ortools/linear_solver/pywraplp.html#Solver.OPTIMAL
        ## OPTIMAL = 0
        ## optimal.
        ## FEASIBLE = 1
        ## feasible, or stopped by limit.
        ## INFEASIBLE = 2
        ## proven infeasible.
        ## UNBOUNDED = 3
        ## proven unbounded.
        ## ABNORMAL = 4
        ## abnormal, i.e., error of some kind.
        ## NOT_SOLVED = 6
        ## not been solved yet.
    assert status == pywraplp.Solver.OPTIMAL
    # Extract values.
    obj_value = solver.Objective().Value() / scale_factor
    y_sol = {s: variables[get_y_index[s]].solution_value() for s in ecal}
    w_sol = {w: variables[get_w_index[w]].solution_value() for w in get_w_index}
    rr_sol = [variables[len(ecal)**2+i].solution_value() for i in range(num_coalitions)]
    print('Objective value =', obj_value)

    return obj_value, y_sol, w_sol, rr_sol, m_sol

def simulate(nodes, ecal, p_edges, q_nodes, vol_scaling, gen_coalitions, \
             coalition_rule = ''):
    results = run_LP_volume(nodes, ecal, p_edges, q_nodes, vol_scaling, \
                            gen_coalitions)
    return (nodes, ecal, p_edges, q_nodes, vol_scaling), results, coalition_rule

def print_results(params, results, coalition_rule):
    nodes, _, p_edges, q_nodes, vol_scaling = params
    obj_value, y_sol, w_sol, rr_sol, m_sol = results
    text = '# PARAMS\n## EPS_EDGE\n'
    text += str(eps_edge) + '\n'
    # Nodes: each edge occupies a line "node q_node",
    # with entries separated by \t.
    text += '## NODES\n'
    for node in nodes:
        text += str(node) + '\t' + str(q_nodes[node]) + '\n'
    # Edges: each edge occupies a line, with "node1 node2 p_edge".
    # As before, entries are separated by \t.
    text += '## EDGES\n'
    for edge in p_edges:
        text += str(edge[0]) + '\t' + str(edge[1]) + '\t' + \
                str(p_edges[edge]) + '\n'
    # Volume scaling: float encoded as string.
    text += '## VOL_SCALING\n' + str(vol_scaling) + '\n'
    # Optimal value: float encoded as string.
    text += '# RESULTS\n## OBJ_VALUE\n' + str(obj_value) + '\n'
    # Solutions for raw entanglement output.
    # y_sol is a dict: edge -> output rate.
    # Each edge occupies a line "node1 node2 rate", with entries separated
    # by \t.
    text += '## Y_SOL\n'
    for edge in y_sol:
        text += str(edge[0]) + '\t' + str(edge[1]) + '\t' + \
                str(y_sol[edge]) + '\n'
    # Solutions for entanglement transfer.
    # w_sol is a dict: (edge1, edge2) -> rate of edge1 used to generate edge2.
    # Each pair of edges occupies a line "node1_1 node1_2 node2_1 node2_2 rate",
    # with entries separated by \t.
    text += '## W_SOL\n'
    for edge1, edge2 in w_sol:
        text += str(edge1[0]) + '\t' + str(edge1[1]) + '\t' + \
                str(edge2[0]) + '\t' + str(edge2[1]) + '\t' + \
                str(w_sol[(edge1, edge2)]) + '\n'
    # Solutions for coalition entanglement output.
    # rr_sol is a list of coalition output rates, where coalitions are
    # generated in lexicographic order (by the order in which nodes are
    # ordered in the nodes list). Use the generator produced by
    # gen_coalitions(nodes) to recover this order.
    # All coalition rates are listed in order in a single line.
    text += '## RR_SOL\n'
    for coal_rate in rr_sol:
        text += str(coal_rate) + '\t'
    text = text[:-1] + '\n'
    # Solutions for coalition effective sizes.
    # m_sol is a list of coalition effective sizes.
    # All effective sizes are listed in order in a single line.
    text += '## M_SOL\n'
    for coal_size in m_sol:
        text += str(coal_size) + '\t'
    text = text[:-1] + '\n'
    # Coalition rule: string.
    text += '# COALITION_RULE\n' + coalition_rule + '\n'
    text = text[:-1]
    return text
    
def save_results(filename, params, results, coalition_rule):
    text = print_results(params, results, coalition_rule)
    with open(filename, 'w') as fn:
        fn.write(text)
    return

def gen_filename(network_type, coalition_rule, suffix = '', folder = ''):
    return 'sim7_data/'+folder+'data_'+network_type+'('+coalition_rule+')'+suffix+'.txt'

if __name__ == '__main__' and False:
    # Limited by state fidelity.
    max_coalsize_fidelity = int(1/np.sqrt(eps_edge))
    print('Fidelity limit:', max_coalsize_fidelity)
    # Value of each link if the link is in a coalition of size n.
    # Distillation factor = approximate fractional loss in raw number of
    # entanglements per link when going from a coalition of size 2 to a
    # coalition of size n.
    # Choose the minmax loss: to decide the loss of a coalition, pick the maximum
    # loss over all links in a coalition, but to decide the loss at size n,
    # pick the minimum loss over all coalitions of size n.
    # This is only approximate since different links in a coalition involve
    # different fractional losses.
    # Note that you lose this fraction for each link in the coalition.
    num_nodes = max_coalsize_fidelity
    vol_scaling = 2
    distillation_factor = lambda n: 0.99**(np.log2(n-1))
    link_value = [round(vol_scaling**n / scipy.special.comb(n, 2) \
                        * (distillation_factor(n) ** scipy.special.comb(n, 2)), 2) \
                  for n in range(2, num_nodes+1)]
    # Value considerations.
    print('Link value:', np.argmax(link_value) + 2)
    print('\t2 ->', link_value, '->', num_nodes)

    # Note that this calculation assumes that the rate vector y_sol is
    # already fixed -- the question is whether we would want to build
    # larger coalitions for a given rate vector.
    # Of course, there is an additional disadvantage of generating bigger
    # coalitions: we need to generate entanglement (y_sol) between faraway
    # nodes, which costs a lot of short-distance entanglement.

if __name__ == '__main__' and False:
    Ns = [2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 40]
    Qs = [0.9, 0.99]
    VSs = [1.5, 2, 3, 4, 5]

    for chain_length in Ns:
        for q_node in Qs:
            for vol_scaling in VSs:
                coalition_rule = 'contiguous'
                suffix = '_N_'+str(chain_length) + \
                         '_Q_'+str(q_node) + \
                         '_VS_'+str(vol_scaling)
                print('('+coalition_rule+')'+suffix)

                folder = 'eps_edge_'+str(eps_edge)+'/'
                if not os.path.exists('sim7_data/'+folder):
                    os.makedirs('sim7_data/'+folder)

                params = gen_params_chain(chain_length, q_node, vol_scaling)
                if coalition_rule == 'all':
                    gen_coalitions = gen_all_coalitions
                elif coalition_rule == 'contiguous':
                    gen_coalitions = gen_coalitions_chain

                try:
                    params, results, coalition_rule = simulate(*params, \
                                                               gen_coalitions, \
                                                               coalition_rule)
                    
                    save_results(gen_filename('chain', coalition_rule, suffix, folder), \
                                 params, results, coalition_rule)
                except AssertionError as e:
                    print(e)
                    continue

if __name__ == '__main__' and True:
    Ns = [2, 3, 4, 5,] # 6, 7]
    Bs = [4, 6, 8, 10] #[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
    Qs = [0.9, 0.99]
    VSs = [1.5, 2, 3, 4, 5]

    for side_num in Ns:
        for bar_factor in Bs:
            for q_node in Qs:
                for vol_scaling in VSs:
                    coalition_rule = 'contiguous'
                    suffix = '_N_'+str(side_num) + \
                             '_B_'+str(bar_factor) + \
                             '_Q_'+str(q_node) + \
                             '_VS_'+str(vol_scaling)
                    print('('+coalition_rule+')'+suffix)

                    folder = 'eps_edge_'+str(eps_edge)+'/'
                    if not os.path.exists('sim7_data/'+folder):
                        os.makedirs('sim7_data/'+folder)

                    params = gen_params_dumbbell(side_num, side_num, bar_factor, q_node, vol_scaling)
                    if coalition_rule == 'all':
                        gen_coalitions = gen_all_coalitions
                    elif coalition_rule == 'contiguous':
                        gen_coalitions = gen_coalitions_dumbbell

                    try:
                        params, results, coalition_rule = simulate(*params, \
                                                                   gen_coalitions, \
                                                                   coalition_rule)
                        
                        save_results(gen_filename('dumbbell', coalition_rule, suffix, folder), \
                                     params, results, coalition_rule)
                    except AssertionError as e:
                        print(e)
                        continue
                
##if __name__ == '__main__' and False:
##    half_lengths = [3, 4]
##    p_edge = 0.3
##    Qs = [0.9, 0.99]
##    VSs = [5,] # [1.5, 2, 3, 4, 5,]
##
##    for half_length in half_lengths:
##        for q_node in Qs:
##            for vol_scaling in VSs:
##                coalition_rule = 'square'
##                suffix = '_N_'+str(half_length) + \
##                         '_P_'+str(p_edge) + \
##                         '_Q_'+str(q_node) + \
##                         '_VS_'+str(vol_scaling)
##                print('('+coalition_rule+')'+suffix)
##
##                params = gen_params_grid(half_length, p_edge, q_node, vol_scaling)
##                if coalition_rule == 'all':
##                    gen_coalitions = gen_all_coalitions
##                elif coalition_rule == 'rectangle':
##                    gen_coalitions = gen_coalitions_grid
##                elif coalition_rule == 'square':
##                    gen_coalitions = gen_coalitions_grid_sq
##
##                params, results, coalition_rule = simulate(*params, \
##                                                           gen_coalitions, \
##                                                           coalition_rule)
##
##                save_results(gen_filename('grid', coalition_rule, suffix), \
##                             params, results, coalition_rule)
    
    

    
