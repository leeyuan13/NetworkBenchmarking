from sim7_plot import *
from sim7_util_dumbbell import gen_params_dumbbell, gen_coalitions_dumbbell

### FOR PAPER

cmap = plt.get_cmap('tab10')

if __name__ == '__main__' and True:
    Ns = [2, 3, 4, 5, 6, 7]
    N_to_plot = 4
    Bs = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 4.0, 6.0, 8.0, 10.0]
    q_node = 0.9
    VSs = [2,]
    VS_to_plot = 2

    eps_edges = [0, 0.01]

    OBJs = dict()
    OBJscaledNS = dict()
    graphs = dict()
    graph_pos = dict()
    max_coal_size = dict()

    for eps_edge in eps_edges:
        folder = 'eps_edge_'+str(eps_edge)+'/'
        for vol_scaling in VSs:
            for bar_factor in Bs:
                for num_side in Ns:
                    coalition_rule = 'contiguous'
                    if coalition_rule == 'all':
                        gen_coalitions = gen_all_coalitions
                    elif coalition_rule == 'contiguous':
                        gen_coalitions = gen_coalitions_dumbbell
                    suffix = '_N_'+str(num_side) + \
                             '_B_'+str(bar_factor) + \
                             '_Q_'+str(q_node) + \
                             '_VS_'+str(vol_scaling)
                    try: # maybe optimization was not completed
                        try:
                            _, _, results, coalition_rule = extract_results('dumbbell', \
                                                                            coalition_rule, \
                                                                            suffix, folder)
                        except FileNotFoundError as e: # sometimes bar_factor is saved as int
                            if int(bar_factor) == bar_factor:
                                suffix = '_N_'+str(num_side) + \
                                         '_B_'+str(int(bar_factor)) + \
                                         '_Q_'+str(q_node) + \
                                         '_VS_'+str(vol_scaling)
                                _, _, results, coalition_rule = extract_results('dumbbell', \
                                                                            coalition_rule, \
                                                                            suffix, folder)
                        obj_value, y_sol, _, rr_sol, _ = results
                        OBJs[(eps_edge, vol_scaling, bar_factor, num_side)] = obj_value
                        OBJscaledNS[(eps_edge, vol_scaling, bar_factor, num_side)] = obj_value / (2*num_side*0.3/2*vol_scaling**2 + 0.3*bar_factor/2*vol_scaling**2)

                        nodes, _, _, _, _ = gen_params_dumbbell(num_side, num_side, bar_factor, q_node, vol_scaling)

                        edges = [(e[0], e[1], y_sol[e]) for e in y_sol]
                        G = nx.Graph()
                        G.add_weighted_edges_from(edges, weight = 'y')
                        nx.set_node_attributes(G, {0: 1, 1: 2} | {2*i: 0 for i in range(1, num_side+1)} | {2*i+1: 3 for i in range(1, num_side+1)}, name = 'subset')
                        graphs[(eps_edge, vol_scaling, bar_factor, num_side)] = G
                        graph_pos[(eps_edge, vol_scaling, bar_factor, num_side)] = nx.multipartite_layout(G, subset_key = 'subset')

                        # Get maximum coalition size.
                        coalitions = gen_coalitions(nodes)
                        ci = 0
                        maxcoal = None
                        for coal in coalitions:
                            if rr_sol[ci] > 0:
                                if maxcoal is None or len(coal) > maxcoal:
                                    maxcoal = len(coal)
                            ci += 1
                        max_coal_size[(eps_edge, vol_scaling, bar_factor, num_side)] = maxcoal

                        # Show coalitions.
                        if vol_scaling == VS_to_plot:
                            coalitions = gen_coalitions(nodes)
                            ci = 0
                            print('-----\n', num_side)
                            for coal in coalitions:
                                if rr_sol[ci] > 0: print(rr_sol[ci], coal)
                                ci += 1
                                
                    except FileNotFoundError as e:
                        print(e)

    print('plotting')

    vol_scaling = VS_to_plot

    plt.figure(figsize = (6, 4.5))
    ax = plt.axes()
    ax.set_yscale('linear')
    for i in range(len(Ns)):
        num_side = Ns[i]
        obj_vec1 = [OBJscaledNS.get((0.01, vol_scaling, bar_factor, num_side)) for bar_factor in Bs]
        plt.plot(Bs, obj_vec1, '.:', color = cmap(i), \
                 label = None)
        obj_vec0 = [OBJscaledNS.get((0, vol_scaling, bar_factor, num_side)) for bar_factor in Bs]
        plt.plot(Bs, obj_vec0, '-', color = cmap(i), alpha = 0.5, \
                 label = r'$M_{side} =$'+str(num_side))

    plt.xlabel('Bar rate / spoke rate')
    plt.ylabel('Utility, scaled to no-swap value')
    leg = plt.legend(loc='upper right')
    leg.set_draggable(True)
    plt.grid()
    plt.ylim(0.9, 5)

    plt.savefig('sim7_data/plot/dumbbell_volume.pdf')

    num_side = 3
    eps_edge = 0.01
    Bs_to_plot = [0, 0.5, 1.0, 1.6, 2.0, 4.0, 6.0, 8.0]
    plt.figure(figsize = (8, 4))
    for i in range(min(8, len(Bs_to_plot))):
        bar_factor = Bs_to_plot[i]
        plt.subplot(2, 4, i+1)
        G = graphs.get((eps_edge, vol_scaling, bar_factor, num_side))
        if G is None:
            continue
        else:
            Gpos = graph_pos[(eps_edge, vol_scaling, bar_factor, num_side)]
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
            plt.title('Bar / spoke = ' + str(bar_factor))

    plt.savefig('sim7_data/plot/dumbbell_eps0.01_ratevector.pdf')

    eps_edge = 0
    plt.figure(figsize = (8, 4))
    for i in range(min(8, len(Bs_to_plot))):
        bar_factor = Bs_to_plot[i]
        plt.subplot(2, 4, i+1)
        G = graphs.get((eps_edge, vol_scaling, bar_factor, num_side))
        if G is None:
            continue
        else:
            Gpos = graph_pos[(eps_edge, vol_scaling, bar_factor, num_side)]
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
            plt.title('Bar / spoke = ' + str(bar_factor))

    plt.savefig('sim7_data/plot/dumbbell_eps0_ratevector.pdf')

    plt.show()

### VISUALIZATION

if __name__ == '__main__' and False:
    eps_edge = 0
    folder = 'eps_edge_'+str(eps_edge)+'/'
    
    Ns = [2, 3, 4, 5, 6, 7,]
    Bs = [0.5, 1, 2,]
    B_to_plot = 0.5
    q_node = 0.9
    VSs = [1.5, 2, 3, 4, 5]
    VS_to_plot = 2

    OBJs = dict()
    OBJscaled = dict()
    graphs = dict()
    graph_pos = dict()
    max_coal_size = dict()

    for vol_scaling in VSs:
        for bar_factor in Bs:
            for num_side in Ns:
                coalition_rule = 'contiguous'
                if coalition_rule == 'all':
                    gen_coalitions = gen_all_coalitions
                elif coalition_rule == 'contiguous':
                    gen_coalitions = gen_coalitions_dumbbell
                suffix = '_N_'+str(num_side) + \
                         '_B_'+str(bar_factor) + \
                         '_Q_'+str(q_node) + \
                         '_VS_'+str(vol_scaling)
                try: # maybe optimization was not completed
                    _, _, results, coalition_rule = extract_results('dumbbell', \
                                                                    coalition_rule, \
                                                                    suffix, folder)
                    obj_value, y_sol, _, rr_sol, _ = results
                    OBJs[(vol_scaling, bar_factor, num_side)] = obj_value
                    OBJscaled[(vol_scaling, bar_factor, num_side)] = obj_value / (2*num_side+1)

                    nodes, _, _, _, _ = gen_params_dumbbell(num_side, num_side, bar_factor, q_node, vol_scaling)

                    edges = [(e[0], e[1], y_sol[e]) for e in y_sol]
                    G = nx.Graph()
                    G.add_weighted_edges_from(edges, weight = 'y')
                    nx.set_node_attributes(G, {0: 1, 1: 2} | {2*i: 0 for i in range(1, num_side+1)} | {2*i+1: 3 for i in range(1, num_side+1)}, name = 'subset')
                    graphs[(vol_scaling, bar_factor, num_side)] = G
                    graph_pos[(vol_scaling, bar_factor, num_side)] = nx.multipartite_layout(G, subset_key = 'subset')

                    # Get maximum coalition size.
                    coalitions = gen_coalitions(nodes)
                    ci = 0
                    maxcoal = None
                    for coal in coalitions:
                        if rr_sol[ci] > 0:
                            if maxcoal is None or len(coal) > maxcoal:
                                maxcoal = len(coal)
                        ci += 1
                    max_coal_size[(vol_scaling, bar_factor, num_side)] = maxcoal

                    # Show coalitions.
                    if vol_scaling == VS_to_plot:
                        coalitions = gen_coalitions(nodes)
                        ci = 0
                        print('-----\n', num_side)
                        for coal in coalitions:
                            if rr_sol[ci] > 0: print(rr_sol[ci], coal)
                            ci += 1
                            
                except FileNotFoundError as e:
                    print(e)

    print('plotting')

    bar_factor = B_to_plot

    plt.figure()
    ax = plt.axes()
    ax.set_yscale('log')
    for vol_scaling in VSs:
        obj_vec = [OBJscaled.get((vol_scaling, bar_factor, num_side)) for num_side in Ns]
        plt.plot(Ns, obj_vec, '.:', label = str(vol_scaling))

    plt.xlabel('Number of periphery nodes per side')
    plt.ylabel('Scaled objective')
    plt.legend()
    plt.grid()

    plt.figure()
    for vol_scaling in VSs:
        maxcoal_vec = [max_coal_size.get((vol_scaling, bar_factor, num_side)) for num_side in Ns]
        plt.plot(Ns, maxcoal_vec, '.:', label = str(vol_scaling))

    plt.xlabel('Number of periphery nodes per side')
    plt.ylabel('Maximum coalition size')
    plt.legend()
    plt.grid()

    vol_scaling = VS_to_plot
    plt.figure(figsize = (4, 8))
    for i in range(min(8, len(Ns))):
        num_side = Ns[i]
        plt.subplot(4, 2, i+1)
        G = graphs.get((vol_scaling, bar_factor, num_side))
        if G is None:
            continue
        else:
            Gpos = graph_pos[(vol_scaling, bar_factor, num_side)]
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
            plt.title('V = ' + str(vol_scaling) + ', B = ' + str(bar_factor) + ', N = ' + str(num_side))

    for i in range(len(Ns)):
        # Show objective value.
        num_side = Ns[i]
        print('-----\n', num_side)
        print(coalition_rule, 'obj =', OBJscaled[(VS_to_plot, B_to_plot, num_side)])

    plt.show()

if __name__ == '__main__' and False:
    eps_edge = 0
    folder = 'eps_edge_'+str(eps_edge)+'/'
    
    Ns = [2, 3, 4, 5, 6, 7,]
    Bs = [0.5, 1, 2,]
    B_to_plot = 0.5
    q_node = 0.9
    VSs = [1.5, 2, 3, 4, 5]
    VS_to_plot = 2

    OBJs = dict()
    OBJscaled = dict()
    graphs = dict()
    graph_pos = dict()
    max_coal_size = dict()

    for vol_scaling in VSs:
        for bar_factor in Bs:
            for num_side in Ns:
                coalition_rule = 'contiguous'
                if coalition_rule == 'all':
                    gen_coalitions = gen_all_coalitions
                elif coalition_rule == 'contiguous':
                    gen_coalitions = gen_coalitions_dumbbell
                suffix = '_N_'+str(num_side) + \
                         '_B_'+str(bar_factor) + \
                         '_Q_'+str(q_node) + \
                         '_VS_'+str(vol_scaling)
                try: # maybe optimization was not completed
                    _, _, results, coalition_rule = extract_results('dumbbell', \
                                                                    coalition_rule, \
                                                                    suffix, folder)
                    obj_value, y_sol, _, rr_sol, _ = results
                    OBJs[(vol_scaling, bar_factor, num_side)] = obj_value
                    OBJscaled[(vol_scaling, bar_factor, num_side)] = obj_value / (2*num_side+1)

                    nodes, _, _, _, _ = gen_params_dumbbell(num_side, num_side, bar_factor, q_node, vol_scaling)

                    edges = [(e[0], e[1], y_sol[e]) for e in y_sol]
                    G = nx.Graph()
                    G.add_weighted_edges_from(edges, weight = 'y')
                    nx.set_node_attributes(G, {0: 1, 1: 2} | {2*i: 0 for i in range(1, num_side+1)} | {2*i+1: 3 for i in range(1, num_side+1)}, name = 'subset')
                    graphs[(vol_scaling, bar_factor, num_side)] = G
                    graph_pos[(vol_scaling, bar_factor, num_side)] = nx.multipartite_layout(G, subset_key = 'subset')

                    # Get maximum coalition size.
                    coalitions = gen_coalitions(nodes)
                    ci = 0
                    maxcoal = None
                    for coal in coalitions:
                        if rr_sol[ci] > 0:
                            if maxcoal is None or len(coal) > maxcoal:
                                maxcoal = len(coal)
                        ci += 1
                    max_coal_size[(vol_scaling, bar_factor, num_side)] = maxcoal

                    # Show coalitions.
                    if vol_scaling == VS_to_plot:
                        coalitions = gen_coalitions(nodes)
                        ci = 0
                        print('-----\n', num_side)
                        for coal in coalitions:
                            if rr_sol[ci] > 0: print(rr_sol[ci], coal)
                            ci += 1
                            
                except FileNotFoundError as e:
                    print(e)

    print('plotting')

    vol_scaling = VS_to_plot

    plt.figure()
    ax = plt.axes()
    ax.set_yscale('log')
    for bar_factor in Bs:
        obj_vec = [OBJscaled.get((vol_scaling, bar_factor, num_side)) for num_side in Ns]
        plt.plot(Ns, obj_vec, '.:', label = str(bar_factor))

    plt.xlabel('Number of periphery nodes per side')
    plt.ylabel('Scaled objective')
    plt.legend()
    plt.grid()

    plt.figure()
    for bar_factor in Bs:
        maxcoal_vec = [max_coal_size.get((vol_scaling, bar_factor, num_side)) for num_side in Ns]
        plt.plot(Ns, maxcoal_vec, '.:', label = str(bar_factor))

    plt.xlabel('Number of periphery nodes per side')
    plt.ylabel('Maximum coalition size')
    plt.legend()
    plt.grid()

    bar_factor = B_to_plot
    plt.figure(figsize = (4, 8))
    for i in range(min(8, len(Ns))):
        num_side = Ns[i]
        plt.subplot(4, 2, i+1)
        G = graphs.get((vol_scaling, bar_factor, num_side))
        if G is None:
            continue
        else:
            Gpos = graph_pos[(vol_scaling, bar_factor, num_side)]
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
            plt.title('V = ' + str(vol_scaling) + ', B = ' + str(bar_factor) + ', N = ' + str(num_side))

    for i in range(len(Ns)):
        # Show objective value.
        num_side = Ns[i]
        print('-----\n', num_side)
        print(coalition_rule, 'obj =', OBJscaled[(VS_to_plot, B_to_plot, num_side)])

    plt.show()

if __name__ == '__main__' and False:
    eps_edge = 0.01
    folder = 'eps_edge_'+str(eps_edge)+'/'
    
    Ns = [2, 3, 4, 5, 7]
    N_to_plot = 4
    Bs = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 4.0, 6.0, 8.0, 10.0]
    q_node = 0.9
    VSs = [1.5, 2, 3, 4, 5]
    VS_to_plot = 2

    OBJs = dict()
    OBJscaledNS = dict()
    graphs = dict()
    graph_pos = dict()
    max_coal_size = dict()

    for vol_scaling in VSs:
        for bar_factor in Bs:
            for num_side in Ns:
                coalition_rule = 'contiguous'
                if coalition_rule == 'all':
                    gen_coalitions = gen_all_coalitions
                elif coalition_rule == 'contiguous':
                    gen_coalitions = gen_coalitions_dumbbell
                suffix = '_N_'+str(num_side) + \
                         '_B_'+str(bar_factor) + \
                         '_Q_'+str(q_node) + \
                         '_VS_'+str(vol_scaling)
                try: # maybe optimization was not completed
                    _, _, results, coalition_rule = extract_results('dumbbell', \
                                                                    coalition_rule, \
                                                                    suffix, folder)
                    obj_value, y_sol, _, rr_sol, _ = results
                    OBJs[(vol_scaling, bar_factor, num_side)] = obj_value
                    OBJscaledNS[(vol_scaling, bar_factor, num_side)] = obj_value / (2*num_side*0.3/2*vol_scaling**2 + 0.3*bar_factor/2*vol_scaling**2)

                    nodes, _, _, _, _ = gen_params_dumbbell(num_side, num_side, bar_factor, q_node, vol_scaling)

                    edges = [(e[0], e[1], y_sol[e]) for e in y_sol]
                    G = nx.Graph()
                    G.add_weighted_edges_from(edges, weight = 'y')
                    nx.set_node_attributes(G, {0: 1, 1: 2} | {2*i: 0 for i in range(1, num_side+1)} | {2*i+1: 3 for i in range(1, num_side+1)}, name = 'subset')
                    graphs[(vol_scaling, bar_factor, num_side)] = G
                    graph_pos[(vol_scaling, bar_factor, num_side)] = nx.multipartite_layout(G, subset_key = 'subset')

                    # Get maximum coalition size.
                    coalitions = gen_coalitions(nodes)
                    ci = 0
                    maxcoal = None
                    for coal in coalitions:
                        if rr_sol[ci] > 0:
                            if maxcoal is None or len(coal) > maxcoal:
                                maxcoal = len(coal)
                        ci += 1
                    max_coal_size[(vol_scaling, bar_factor, num_side)] = maxcoal

                    # Show coalitions.
                    if vol_scaling == VS_to_plot:
                        coalitions = gen_coalitions(nodes)
                        ci = 0
                        print('-----\n', num_side)
                        for coal in coalitions:
                            if rr_sol[ci] > 0: print(rr_sol[ci], coal)
                            ci += 1
                            
                except FileNotFoundError as e:
                    print(e)

    print('plotting')

    vol_scaling = VS_to_plot

    plt.figure()
    ax = plt.axes()
    ax.set_yscale('linear')
    for num_side in Ns:
        obj_vec = [OBJscaledNS.get((vol_scaling, bar_factor, num_side)) for bar_factor in Bs]
        plt.plot(Bs, obj_vec, '.:', label = 'N = '+str(num_side))

    plt.xlabel('Bar rate / spoke rate')
    plt.ylabel('Objective value / no-swap value')
    plt.legend()
    plt.grid()

    plt.figure()
    for num_side in Ns:
        maxcoal_vec = [max_coal_size.get((vol_scaling, bar_factor, num_side)) for bar_factor in Bs]
        plt.plot(Bs, maxcoal_vec, '.:', label = 'N = '+str(num_side))

    plt.xlabel('Bar rate / spoke rate')
    plt.ylabel('Maximum coalition size')
    plt.legend()
    plt.grid()

    num_side = N_to_plot
    plt.figure(figsize = (8, 4))
    for i in range(min(8, len(Bs))):
        bar_factor = Bs[2*i]
        plt.subplot(2, 4, i+1)
        G = graphs.get((vol_scaling, bar_factor, num_side))
        if G is None:
            continue
        else:
            Gpos = graph_pos[(vol_scaling, bar_factor, num_side)]
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
            plt.title('V = ' + str(vol_scaling) + ', B = ' + str(bar_factor) + ', N = ' + str(num_side))

    for i in range(len(Bs)):
        # Show objective value.
        bar_factor = Bs[i]
        print('-----\n', bar_factor)
        print(coalition_rule, 'obj =', OBJscaledNS[(VS_to_plot, bar_factor, N_to_plot)])

    plt.show()

        
