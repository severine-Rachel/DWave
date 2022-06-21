from dwave.system import DWaveSampler, EmbeddingComposite
from collections import defaultdict
import networkx as nx
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
from collections import defaultdict
import greedy
import dwave.inspector
# Size of the grid
size = (3,4)

lagrangian_objective = 1
# lagrangien structure
lagrangian_struct = 5
lagrangian_struct_start = lagrangian_struct
lagrangian_struct_end = lagrangian_struct
lagrangian_struct_path = lagrangian_struct

def coord_to_index(coord):
    return coord[0] + coord[1]*size[0]

def index_to_coord(index):
    return (index % size[0], index // size[0])

def unique_edge(left_node_index, right_node_index):
    return (min(left_node_index, right_node_index), max(left_node_index, right_node_index))

def unique_edge_edge(edge1, edge2):
    if (edge1[0] < edge2[0]):
        return (edge1, edge2)
    elif (edge1[0] == edge2[0]):
        if (edge1[1] < edge2[1]):
            return (edge1, edge2)
        else:
            return (edge2, edge1)
    else:
        return (edge2, edge1)

assert (2,3) == index_to_coord(coord_to_index((2,3)))

# Start and end node
start_point = (1,0)
end_point = (2,3)

start_point_index = coord_to_index(start_point)
end_point_index = coord_to_index(end_point)

for j in range(0, size[1]):
    for i in range(0, size[0]):
        if (i, j) == start_point:
            print("S", end="")
        elif (i, j) == end_point:
            print("E", end="")
        else:
            print(".", end="")
    print("")

def plot_graph(graph, filename="graph_plot.png", path=[]):
    pos = nx.spring_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_color='grey')
    #nx.draw_networkx_nodes(graph, pos, nodelist=[node for node in graph.nodes if node == coord_to_index(start_point[0], start_point[1])], node_color='red')
    #nx.draw_networkx_nodes(graph, pos, nodelist=[node for node in graph.nodes if node == coord_to_index(end_point[0], end_point[1])], node_color='black')
    nx.draw_networkx_nodes(graph, pos, nodelist=[node for node in G2.nodes if node[0] == start_point_index or node[1] == start_point_index], node_color='r')
    nx.draw_networkx_nodes(graph, pos, nodelist=[node for node in G2.nodes if node[0] == end_point_index or node[1] == end_point_index], node_color='c')
    if len(path) > 0:
        nx.draw_networkx_nodes(graph, pos, nodelist=path, node_color='g')
    nx.draw_networkx_edges(graph, pos)
    
    #nx.draw_networkx_labels(graph, pos)
    plt.show()
    #plt.savefig(filename)

# Initial graph
G1 = nx.Graph()
G1.add_nodes_from([(coord_to_index((i, j)), {"x":i,"y":j}) for i in range (0, size[0]) for j in range(0, size[1])])
for i in range(0, size[0]-1):
    G1.add_edges_from([unique_edge(coord_to_index((i, j)), coord_to_index((i+1, j))) for j in range(0, size[1])])
for j in range(0, size[1]-1):
    G1.add_edges_from([unique_edge(coord_to_index((i, j)), coord_to_index((i, j+1))) for i in range(0, size[0])])

# Edge graph : edge becomes a node and path from edge to edge become vertex
G2 = nx.Graph()
G2.add_nodes_from([unique_edge(u, v) for (u, v) in G1.edges])
for (u, v) in G1.edges:
    for (u_i, v_i) in G1.edges:
        if ((u in [u_i, v_i] or v in [u_i, v_i])):
            G2.add_edge(unique_edge(u, v), unique_edge(u_i, v_i))

# plot_graph(G2)

Q = defaultdict(int)
nodes_from_start = [node for node in G2.nodes if node[0] == start_point_index or node[1] == start_point_index]
nodes_from_end = [node for node in G2.nodes if node[0] == end_point_index or node[1] == end_point_index]

# ajout des noeuds dans le graphe
for node in G2.nodes:
    Q[(node, node)] = 1 * lagrangian_objective

# contrainte sur le dÃ©part
for node in nodes_from_start:
    Q[(node, node)] += -1 * lagrangian_struct_start
    for node2 in nodes_from_start:
        if (node != node2):
            Q[(node, node2)] += 1 * lagrangian_struct_start

# contrainte sur l'arrivÃ©e
for node in nodes_from_end:
    Q[(node, node)] += -1 * lagrangian_struct_end
    for node2 in nodes_from_end:
        if (node != node2):
            Q[(node, node2)] += 1 * lagrangian_struct_end


for node in G2.nodes:
    # arcs dits entrants
    entring_node = [enode for enode in G2.adj[node] if enode[0] == node[0] or enode[1] == node[0]]
    exiting_node = [enode for enode in G2.adj[node] if enode[0] == node[1] or enode[1] == node[1]]
    #print ("Les points entrants de ", node," sont ", entring_node)
    #print ("Les points sorants de ", node," sont ", exiting_node)
    ## contrainte sur les arcs de passage, reprÃ©sentÃ©s par des noeuds dans G2
    # On vÃ©rifie que le nombre total de voisins est bien Ã©gal Ã  la somme des sortants et entrants
    assert G2.degree[node] == len(entring_node) + len(exiting_node)
    # Maintenant on dÃ©finit les contraintes : 
    # Sum_{entrants}(xi) + Sum_{sortants}(xi) + 2*Sum_{entrants,j>i}(xi*xj) + 2*Sum_{sortants,j>i}(xi*xj) - 2*Sum_{entrants}(xi)*Sum_{sortants}(xi)
    # Somme des xi entrants
#for node in G2.nodes:
    for enode in entring_node:
        Q[(enode, enode)] += 1 * lagrangian_struct_path
    #for enode in entring_node:
    #    print("valeur points entrant :", (enode, enode), Q[(enode, enode)], "pour la valeur des noeurs entrants : ", entring_node)

    # Somme des xi sortants
    for enode in exiting_node:
        Q[(enode, enode)] += 1 * lagrangian_struct_path
#    for enode in exiting_node:
#        print("valeur points sortants :", (enode, enode), Q[(enode, enode)])
    # 2*Sum_{entrants,j>i}(xi*xj)

    for enode in entring_node:
        for enode2 in entring_node: # Avec la double boucle on compte deux fois, donc on divise le coeff par 2
             if (enode != enode2):
                 Q[(enode, enode2)] += 1 * lagrangian_struct_path
    for enode in entring_node:
        for enode2 in entring_node: 
            if (enode != enode2):
                print("valeur arcs entrant :", (enode, enode2), Q[(enode, enode2)])
    # 2*Sum_{sortants,j>i}(xi*xj)


    for enode in exiting_node:
        for enode2 in exiting_node: # Avec la double boucle on compte deux fois, donc on divise le coeff par 2
             if (enode != enode2):
                 Q[(enode, enode2)] += 1 * lagrangian_struct_path
                 print(Q[(enode, enode2)])
    # -2*Sum_{entrants}(xi)*Sum_{sortants}(xi)

    for enode in entring_node:
        for enode2 in exiting_node:
            Q[(enode, enode2)] += -2 * lagrangian_struct_path

    ## Contrainte sur le nombre d'arcs sortants
    for i in range(0, len(exiting_node)-1):
        for j in range(i+1, len(exiting_node)):
            Q[(exiting_node[i], exiting_node[j])] += 1 * lagrangian_struct_path

print(Q)
plot_graph(G2)


solver = greedy.SteepestDescentSolver()

sampler = EmbeddingComposite(DWaveSampler())
sampleset = sampler.sample_qubo(Q,
                               chain_strength=1,
                               num_reads=100,
                               label='spath')

print(sampleset)

# On affiche les trois premiers chemins
print(sampleset.variables)
for i in range(0, 3):
    sample = sampleset.record.sample[i]
    print(sample)
    path = [sampleset.variables[i] for i in range(0, len(sample)) if sample[i] == 1]
    print("path: " + str(path))
    #plot_graph(G2, path=path)


    graph = [[0 for x in range(size[1])] for y in range(size[0]) ]
    for node in path:
        coord0 = index_to_coord(node[0])
        coord1 = index_to_coord(node[1])
        graph[coord0[0]][coord0[1]] = 1
        graph[coord1[0]][coord1[1]] = 1

    for j in range(0, size[1]):
        for i in range(0, size[0]):
            if (i, j) == start_point:
                print("S", end="")
            elif (i, j) == end_point:
                print("E", end="")
            elif graph[i][j] == 1:
                print("X", end="")
            else:
                print(".", end="")
        print("")
filename = "spath.png"
plt.savefig(filename, bbox_inches='tight')
print("\nYour plot is saved to {}".format(filename))
#dwave.inspector.show(sampleset)