from itertools import chain, combinations
import time

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

# Construct indices for LP variables.
# Get possible entangled links.
def get_ecal(nodes):
    ecal = [get_edge(nodes[i1], nodes[i2]) for i1 in range(len(nodes)) \
                                    for i2 in range(i1+1, len(nodes))]
    return ecal

# Get default generator of coalitions.
def gen_all_coalitions(nodes):
    # Number of potential coalitions with >= 2 nodes is
    # 2**num_nodes-1-num_nodes.
    # If necessary, we could further restrict the size of coalitions.
    max_coalition_size = len(nodes) + 1
    coalitions = chain.from_iterable(combinations(nodes, r) \
                                     for r in range(2, max_coalition_size))
    # Note that itertools.combinations produces combinations in lexicographic
    # order (according to their positions in nodes).
    # Hence, different versions of this generator produce coalitions
    # in the same order.
    # Each combination is also emitted in lexicographic order.
    return coalitions    
