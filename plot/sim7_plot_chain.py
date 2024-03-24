from sim7_plot import *
from sim7_util_chain import gen_coalitions_chain

### FOR PAPER

cmap = plt.get_cmap('tab10')

if __name__ == '__main__' and True:
    Ns = [2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 40]
    q_node = 0.9
    VSs = [1.5, 2, 3, 4]
    VS_to_plot = 2

    eps_edges = [0, 0.01]

    OBJs = dict()
    OBJscaled = dict()
    graphs = dict()
    graph_pos = dict()
    max_coal_size = dict()

    for eps_edge in eps_edges:
        folder = 'eps_edge_'+str(eps_edge)+'/'
        for vol_scaling in VSs:
            values = []
            values_scaled = []
            maxcoals = []
            for chain_length in Ns:
                coalition_rule = 'contiguous'
                if coalition_rule == 'all':
                    gen_coalitions = gen_all_coalitions
                elif coalition_rule == 'contiguous':
                    gen_coalitions = gen_coalitions_chain
                suffix = '_N_'+str(chain_length) + \
                         '_Q_'+str(q_node) + \
                         '_VS_'+str(vol_scaling)
                try: # maybe optimization was not completed
                    _, _, results, coalition_rule = extract_results('chain', \
                                                                    coalition_rule, \
                                                                    suffix, folder)
                    obj_value, y_sol, _, rr_sol, _ = results
                    values.append(obj_value)
                    values_scaled.append(obj_value / 0.3 / (chain_length-1) / max(vol_scaling**2 * 1/2, vol_scaling * 1))

                    edges = [(e[0], e[1], y_sol[e]) for e in y_sol]
                    G = nx.Graph()
                    G.add_weighted_edges_from(edges, weight = 'y')
                    graphs[(eps_edge, vol_scaling, chain_length)] = G
                    graph_pos[(eps_edge, vol_scaling, chain_length)] = nx.circular_layout(G)

                    # Get maximum coalition size.
                    nodes = list(range(chain_length))
                    coalitions = gen_coalitions(nodes)
                    ci = 0
                    maxcoal = None
                    for coal in coalitions:
                        if rr_sol[ci] > 0:
                            if maxcoal is None or len(coal) > maxcoal:
                                maxcoal = len(coal)
                        ci += 1
                    maxcoals.append(maxcoal)

                    # Show coalitions.
                    if vol_scaling == VS_to_plot:
                        nodes = list(range(chain_length))
                        coalitions = gen_coalitions(nodes)
                        ci = 0
                        print('-----\n', chain_length)
                        for coal in coalitions:
                            if rr_sol[ci] > 0: print(rr_sol[ci], coal)
                            ci += 1
                            
                except FileNotFoundError as e:
                    print(e)
                    values.append(None)
                    values_scaled.append(None)
                    maxcoals.append(None)
                
            OBJs[(eps_edge, vol_scaling)] = values
            OBJscaled[(eps_edge, vol_scaling)] = values_scaled

            max_coal_size[(eps_edge, vol_scaling)] = maxcoals

    # Plot volume.
    plt.figure(figsize = (6, 4.5))
    ax = plt.axes()
    ax.set_yscale('log')
    for i in range(len(VSs)):
        vol_scaling = VSs[i]
        plt.plot(Ns, OBJscaled[(0, vol_scaling)], '-', alpha = 0.5, \
                      color = cmap(i), \
                      label = r'$\beta =$'+str(vol_scaling))
        plt.plot(Ns, OBJscaled[(0.01, vol_scaling)], '.:', \
                      color = cmap(i), \
                      label = None)

    plt.xlabel('Chain length')
    plt.ylabel('Utility, scaled to no-swap value')
    leg = plt.legend()
    leg.set_draggable(True)
    plt.grid()
    plt.ylim(0.9, 5e3)

    plt.savefig('sim7_data/plot/chain_volume.pdf')

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex = 'col', figsize = (6, 4.5))
    for i in range(len(VSs)):
        vol_scaling = VSs[i]
        ax1.plot(Ns, max_coal_size[(0, vol_scaling)], '-', alpha = 0.5, \
                 color = cmap(i), \
                 label = r'$\beta =$'+str(vol_scaling))
        ax2.plot(Ns, max_coal_size[(0.01, vol_scaling)], '.:', \
                 color = cmap(i), \
                 label = r'$\beta =$'+str(vol_scaling))
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.xlabel('Chain length')
    plt.ylabel('Maximum coalition size')
    ax1.legend()
    ax2.legend()
    ax1.grid()
    ax2.grid()
    ax1.set_ylim(1, 12)
    ax2.set_ylim(1, 12)
    ax1.set_xlim(1, 25)

    plt.savefig('sim7_data/plot/chain_maxcoalsize.pdf')

    vol_scaling = 2
    eps_edge = 0
    plt.figure(figsize = (8, 4))
    for i in range(8):
        chain_length = Ns[i+3]
        plt.subplot(2, 4, i+1)
        G = graphs.get((eps_edge, vol_scaling, chain_length))
        if G is None:
            continue
        else:
            Gpos = graph_pos[(eps_edge, vol_scaling, chain_length)]
            nx.draw_networkx_nodes(G, Gpos, node_size = 25, node_color = '#FB9E9D')
            nx.draw_networkx_labels(G, Gpos, \
                                    labels = {node: str(node) for node in G.nodes}, \
                                    font_color = 'k', font_size = 5)
            edgelist = []
            edgecolors = []
            for edge in G.edges:
                edgelist.append(edge)
                edgecolors.append(G.edges[edge]['y'])
            max_y = max(edgecolors)
            edgecolors = [(0, 0, 0, y/max_y) for y in edgecolors]
            nx.draw_networkx_edges(G, Gpos, edgelist = edgelist, edge_color = edgecolors)
            plt.title(r'$M =$' + str(chain_length))

    plt.savefig('sim7_data/plot/chain_eps0_ratevector.pdf')

    plt.figure(figsize = (8, 4))
    eps_edge = 0.01
    for i in range(8):
        chain_length = Ns[i+3]
        plt.subplot(2, 4, i+1)
        G = graphs.get((eps_edge, vol_scaling, chain_length))
        if G is None:
            continue
        else:
            Gpos = graph_pos[(eps_edge, vol_scaling, chain_length)]
            nx.draw_networkx_nodes(G, Gpos, node_size = 25, node_color = '#FB9E9D')
            nx.draw_networkx_labels(G, Gpos, \
                                    labels = {node: str(node) for node in G.nodes}, \
                                    font_color = 'k', font_size = 5)
            edgelist = []
            edgecolors = []
            for edge in G.edges:
                edgelist.append(edge)
                edgecolors.append(G.edges[edge]['y'])
            max_y = max(edgecolors)
            edgecolors = [(0, 0, 0, y/max_y) for y in edgecolors]
            nx.draw_networkx_edges(G, Gpos, edgelist = edgelist, edge_color = edgecolors)
            plt.title(r'$M =$' + str(chain_length))

    plt.savefig('sim7_data/plot/chain_eps0.01_ratevector.pdf')

    plt.show()


### VISUALIZATION

if __name__ == '__main__' and False:
    eps_edge = 0.01
    folder = 'eps_edge_'+str(eps_edge)+'/'
    
    Ns = [2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 40]
    q_node = 0.9
    VSs = [1.5, 2, 3, 4, 5]
    VS_to_plot = 2

    OBJs = dict()
    OBJscaled = dict()
    graphs = dict()
    graph_pos = dict()
    max_coal_size = dict()

    for vol_scaling in VSs:
        values = []
        values_scaled = []
        maxcoals = []
        for chain_length in Ns:
            coalition_rule = 'contiguous'
            if coalition_rule == 'all':
                gen_coalitions = gen_all_coalitions
            elif coalition_rule == 'contiguous':
                gen_coalitions = gen_coalitions_chain
            suffix = '_N_'+str(chain_length) + \
                     '_Q_'+str(q_node) + \
                     '_VS_'+str(vol_scaling)
            try: # maybe optimization was not completed
                _, _, results, coalition_rule = extract_results('chain', \
                                                                coalition_rule, \
                                                                suffix, folder)
                obj_value, y_sol, _, rr_sol, _ = results
                values.append(obj_value)
                values_scaled.append(obj_value / (chain_length-1))

                edges = [(e[0], e[1], y_sol[e]) for e in y_sol]
                G = nx.Graph()
                G.add_weighted_edges_from(edges, weight = 'y')
                graphs[(vol_scaling, chain_length)] = G
                graph_pos[(vol_scaling, chain_length)] = nx.circular_layout(G)

                # Get maximum coalition size.
                nodes = list(range(chain_length))
                coalitions = gen_coalitions(nodes)
                ci = 0
                maxcoal = None
                for coal in coalitions:
                    if rr_sol[ci] > 0:
                        if maxcoal is None or len(coal) > maxcoal:
                            maxcoal = len(coal)
                    ci += 1
                maxcoals.append(maxcoal)

                # Show coalitions.
                if vol_scaling == VS_to_plot:
                    nodes = list(range(chain_length))
                    coalitions = gen_coalitions(nodes)
                    ci = 0
                    print('-----\n', chain_length)
                    for coal in coalitions:
                        if rr_sol[ci] > 0: print(rr_sol[ci], coal)
                        ci += 1
                        
            except FileNotFoundError as e:
                print(e)
                values.append(None)
                values_scaled.append(None)
                maxcoals.append(None)
            
        OBJs[vol_scaling] = values
        OBJscaled[vol_scaling] = values_scaled

        max_coal_size[vol_scaling] = maxcoals

    print('plotting')

##    plt.figure()
##    for vol_scaling in VSs:
##        plt.plot(Ns, OBJs[vol_scaling], '.--', label = str(vol_scaling))
##
##    plt.legend()
##    plt.grid()

    plt.figure()
    ax = plt.axes()
    ax.set_yscale('log')
    for vol_scaling in VSs:
        plt.plot(Ns, OBJscaled[vol_scaling], '.:', label = str(vol_scaling))

    # Threshold to go from shallow to deep computations.
    xl = plt.xlim()
    plotNs = np.linspace(xl[0], xl[1], 100)
##    # Compute objective for size-3 coalitions with depth-1 computations.
##    plt.plot(plotNs, np.exp(np.log(plotNs)/(plotNs-1)) * 0.3 * 3/1, 'k-', \
##             alpha = 0.2, label = 'depth threshold')
##    # Compute objective for size-3 coalitions with depth-3 computations.
##    plt.plot(plotNs, np.exp(np.log(plotNs*(plotNs+1)*(plotNs+2)/60)/(plotNs-3)) * 0.3 * 3/3, 'r-', \
##             alpha = 0.2, label = 'size threshold')
    plt.xlim(xl)

    plt.xlabel('Chain length')
    plt.ylabel('Scaled objective')
    plt.legend()
    plt.grid()

    plt.figure()
    for vol_scaling in VSs:
        plt.plot(Ns, max_coal_size[vol_scaling], '.:', label = str(vol_scaling))

    plt.xlabel('Chain length')
    plt.ylabel('Max coalition size')
    plt.legend()
    plt.grid()

    vol_scaling = VS_to_plot
    plt.figure(figsize = (8, 4))
    for i in range(8):
        chain_length = Ns[i+3]
        plt.subplot(2, 4, i+1)
        G = graphs.get((vol_scaling, chain_length))
        if G is None:
            continue
        else:
            Gpos = graph_pos[(vol_scaling, chain_length)]
            nx.draw_networkx_nodes(G, Gpos, node_size = 25)
            nx.draw_networkx_labels(G, Gpos, \
                                    labels = {node: str(node) for node in G.nodes}, \
                                    font_color = 'w', font_size = 5)
            edgelist = []
            edgecolors = []
            for edge in G.edges:
                edgelist.append(edge)
                edgecolors.append(G.edges[edge]['y'])
            max_y = max(edgecolors)
            edgecolors = [(0, 0, 0, y/max_y) for y in edgecolors]
            nx.draw_networkx_edges(G, Gpos, edgelist = edgelist, edge_color = edgecolors)
            plt.title('V = ' + str(vol_scaling) + ', N = ' + str(chain_length))

    for i in range(len(Ns)):
        # Show objective value.
        chain_length = Ns[i]
        print('-----\n', chain_length)
        print(coalition_rule, 'obj =', OBJscaled[VS_to_plot][i])

    plt.show()
