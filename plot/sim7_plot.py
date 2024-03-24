import numpy as np
import matplotlib.pyplot as plt
from itertools import chain, combinations
import scipy.optimize
import time
import seaborn as sns
import networkx as nx

from sim7 import get_edge, get_ecal, gen_all_coalitions, gen_filename

plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.it'] = 'Arial:italic'
plt.rcParams['mathtext.rm'] = 'Arial'

plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'

def read_node(text):
    try:
        return int(text)
    except ValueError:
        return tuple(map(int, text[1:-1].split(', ')))

def read_results(text):
    lines = text.strip().split('\n')

    key_words = ['## EPS_EDGE', '## NODES', '## EDGES', '## VOL_SCALING', \
                 '## OBJ_VALUE', '## Y_SOL', '## W_SOL', '## RR_SOL', \
                 '## M_SOL', '# COALITION_RULE']

    indices = [lines.index(word) for word in key_words]

    toprint = False

    eps_edge = float(lines[indices[0]+1])
    if toprint: print(eps_edge)

    lines_nodes = lines[indices[1]+1:indices[2]]
    nodes = []
    q_nodes = dict()
    for line in lines_nodes:
        entries = line.strip().split('\t')
        node = read_node(entries[0])
        q = float(entries[1])
        nodes.append(node)
        q_nodes[node] = q
    if toprint: print(nodes, q_nodes)

    lines_edges = lines[indices[2]+1:indices[3]]
    p_edges = dict()
    for line in lines_edges:
        entries = line.strip().split('\t')
        node1 = read_node(entries[0])
        node2 = read_node(entries[1])
        p = float(entries[2])
        p_edges[get_edge(node1, node2)] = p
    if toprint: print(p_edges)

    vol_scaling = float(lines[indices[3]+1])
    if toprint: print(vol_scaling)

    obj_value = float(lines[indices[4]+1])
    if toprint: print(obj_value)

    lines_ysol = lines[indices[5]+1:indices[6]]
    y_sol = dict()
    for line in lines_ysol:
        entries = line.strip().split('\t')
        node1 = read_node(entries[0])
        node2 = read_node(entries[1])
        y = float(entries[2])
        y_sol[get_edge(node1, node2)] = y
    if toprint: print(y_sol)

    lines_wsol = lines[indices[6]+1:indices[7]]
    w_sol = dict()
    for line in lines_wsol:
        entries = line.strip().split('\t')
        node11 = read_node(entries[0])
        node12 = read_node(entries[1])
        node21 = read_node(entries[2])
        node22 = read_node(entries[3])
        w = float(entries[4])
        w_sol[(get_edge(node11, node12), get_edge(node21, node22))] = w
    if toprint: print(w_sol)

    line = lines[indices[7]+1]
    rr_sol = []
    entries = line.strip().split('\t')
    for entry in entries:
        rr_sol.append(float(entry))
    if toprint: print(rr_sol)

    line = lines[indices[8]+1]
    m_sol = []
    entries = line.strip().split('\t')
    for entry in entries:
        m_sol.append(float(entry))
    if toprint: print(m_sol)

    line = lines[indices[9]+1]
    coalition_rule = line.strip()

    # Return (eps_edge, params, results, coalition_rule).
    return eps_edge, (nodes, p_edges, q_nodes, vol_scaling), \
           (obj_value, y_sol, w_sol, rr_sol, m_sol), coalition_rule
    
def extract_results(network_type, coalition_rule, suffix, folder = ''):
    with open(gen_filename(network_type, coalition_rule, suffix, folder), 'r') as fn:
        text = fn.read()
    return read_results(text)
