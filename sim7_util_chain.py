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

def gen_params_chain(chain_length, q_node, vol_scaling):
    nodes = list(range(chain_length))
    edges = [(i, i+1) for i in range(chain_length-1)]
    p_edges = {s: 0.3 for s in edges}
    ecal = get_ecal(nodes)
    q_nodes = {node: q_node for node in nodes}
    return nodes, ecal, p_edges, q_nodes, vol_scaling

# Returns a generator of coalitions, using a heuristic.
def gen_coalitions_chain(nodes):
    # Heuristic: only contiguous nodes can form a coalition.
    # Note that coalitions are produced in a different order from
    # gen_all_coalitions. Here, we order coalitions like a dictionary;
    # in gen_all_coalitions we order coalitions by size first.
    for i in range(len(nodes)):
        for j in range(i+1, len(nodes)):
            yield tuple(range(i, j+1))
