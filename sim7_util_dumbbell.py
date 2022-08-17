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

def gen_params_dumbbell(num_left, num_right, bar_factor, q_node, vol_scaling):
    nodes = [0, 1] + [2*i for i in range(1, num_left+1)] + [2*i+1 for i in range(1, num_right+1)]
    nodes.sort()
    ecal = get_ecal(nodes)
    # For Python 3.9 and above.
    # p_edges = {(0, 2*i): 0.3 for i in range(1, num_left+1)} | \
    #           {(1, 2*i+1): 0.3 for i in range(1, num_right+1)} | \
    #           {(0, 1): 0.3 * bar_factor}
    # For backward compatibility.
    p_edges = {(0, 2*i): 0.3 for i in range(1, num_left+1)}
    p_edges.update({(1, 2*i+1): 0.3 for i in range(1, num_right+1)})
    p_edges.update({(0, 1): 0.3 * bar_factor})
    q_nodes = {node: q_node for node in nodes}
    return nodes, ecal, p_edges, q_nodes, vol_scaling

# Generate nonempty subsets of [1, 2, ..., r].
def gen_nonempty_subsets(r):
    if r < 1:
        return iter(())
    else:
        return chain.from_iterable(combinations(range(1, r+1), q) \
                                   for q in range(1, r+1))

def gen_coalitions_dumbbell(nodes):
    # Heuristic: only contiguous nodes can form a coalition.
    # To generate contiguous nodes in a dumbbell, note that at least one of
    # the cores (0 and/or 1) must be in the coalition.
    # Why? Recall that coalitions must have at least two nodes.
    # We generate coalitions in this way:
    # 1.    The first coalition is (0, 1).
    # 2.    Next consider coalitions involving the left core (0) and its
    #       neighbors. This excludes the right periphery nodes. We also
    #       exclude (0, 1), since this is already accounted for. We generate
    #       the power set of left periphery nodes, excluding the empty set.
    #       Let X be a nonempty subset of left periphery nodes. The set of
    #       relevant coalitions includes (0, X) and (0, 1, X).
    # 3.    Next consider coalitions involving the right core (1) and its
    #       neighbors. This excludes the left periphery nodes. We also exclude
    #       (0, 1). The procedure is the same as before.
    # 4.    Finally, consider coalitions involving the left and right cores,
    #       at least one of the left periphery nodes and at least one of the
    #       right periphery nodes. We generate the power sets of left
    #       and right periphery nodes, excluding empty sets. Let X and Y be
    #       nonempty subsets of left and right periphery nodes respectively.
    #       Then (0, 1, X, Y) is a coalition.

    # Get num_left, num_right.
    nodes_descending = sorted(nodes, reverse = True)
    num_left, num_right = None, None
    for node in nodes_descending:
        if num_left is None and (node % 2 == 0):
            num_left = node // 2
        if num_right is None and (node % 2 == 1):
            num_right = (node - 1) // 2

    # Generate coalitions.
    yield (0, 1)

    pset_left = gen_nonempty_subsets(num_left)
    for subset_left in pset_left:
        yield (0,) + tuple(2*x for x in subset_left)
        yield (0, 1,) + tuple(2*x for x in subset_left)

    pset_right = gen_nonempty_subsets(num_right)
    for subset_right in pset_right:
        yield (1,) + tuple(2*x+1 for x in subset_right)
        yield (0, 1,) + tuple(2*x+1 for x in subset_right)

    pset_left = gen_nonempty_subsets(num_left)
    for subset_left in pset_left:
        pset_right = gen_nonempty_subsets(num_right)
        for subset_right in pset_right:
            yield (0, 1,) + tuple(sorted([2*x for x in subset_left] + \
                                         [2*x+1 for x in subset_right]))

    
