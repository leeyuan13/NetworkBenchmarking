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

def gen_params_grid(half_length, p_edge, q_node, vol_scaling):
    # Only square grids, with an odd number of grid points per side.
    nodes = [(i, j) for i in range(-half_length, half_length+1) \
             for j in range(-half_length, half_length+1)]
    ecal = get_ecal(nodes)
    p_edges = {((i, j), (i+1, j)): p_edge for i in range(-half_length, half_length) \
               for j in range(-half_length, half_length+1)} | \
               {((i, j), (i, j+1)): p_edge for i in range(-half_length, half_length+1) \
                for j in range(-half_length, half_length)}
    q_nodes = {node: q_node for node in nodes}
    return nodes, ecal, p_edges, q_nodes, vol_scaling

# Generates coalitions, using a heuristic.
def gen_coalitions_grid(nodes):
    # Assumes nodes form a grid with cell length 1.
    # Heuristic: only nodes in a rectangle can form a coalition.
    xs = [node[0] for node in nodes]
    ys = [node[1] for node in nodes]
    xmin = min(xs)
    xmax = max(xs)
    ymin = min(ys)
    ymax = max(ys)

    for xi in range(xmin, xmax+1):
        for xj in range(xi, xmax+1):
            for yi in range(ymin, ymax+1):
                for yj in range(yi, ymax+1):
                    # Coalition must have at least two nodes.
                    if xi == xj and yi == yj: continue
                    coal_nodes = tuple((xx, yy) for xx in range(xi, xj+1) \
                                  for yy in range(yi, yj+1))
                    yield coal_nodes

# Does not reduce the number of variables by much relative to
# gen_coalitions_grid.
def gen_coalitions_grid_sq(nodes):
    # Assumes nodes form a grid with cell length 1.
    # Heuristic: only nodes in a square can form a coalition.
    xs = [node[0] for node in nodes]
    ys = [node[1] for node in nodes]
    xmin = min(xs)
    xmax = max(xs)
    ymin = min(ys)
    ymax = max(ys)

    for xi in range(xmin, xmax+1):
        for xj in range(xi+1, xmax+1): # coalition must have at least two nodes
            side_length = xj - xi + 1
            for yi in range(ymin, ymax+1-(side_length-1)):
                yj = yi + (side_length-1)
                coal_nodes = tuple((xx, yy) for xx in range(xi, xj+1) \
                              for yy in range(yi, yj+1))
                yield coal_nodes
