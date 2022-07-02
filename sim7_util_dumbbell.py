from itertools import chain, combinations
import time

from sim7_util import get_edge, get_ecal

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

get_edge = lambda i, j: (i, j) if i < j else (j, i)

def get_params_dumbbell(num_per_side, bar_factor, q_node, vol_scaling):
    nodes = list(range(2 + 2*num_per_side))
    ecal = get_ecal(nodes)
    p_edges = {(0, 2*i): 0.3 for i in range(1, num_per_side+1)} | \
              {(1, 2*i+1): 0.3 for i in range(1, num_per_side+1)} | \
              {(0, 1): 0.3 * bar_factor}
    q_nodes = {node: q_node for node in nodes}
    return nodes, ecal, p_edges, q_nodes, vol_scaling

def gen_coalitions_dumbbell(nodes):
    pass
