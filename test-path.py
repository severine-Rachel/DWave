print ("etape 1 projet de recherche de chemin le plus court")
import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict
import greedy

# Size of the grid
size = (3,4)
start_point = (0,0)
end_point = (2,3)

lagrangian_objective = 60
lagrangian_extrem = 100
lagrangian_struct_path = 100
lagrangien_sortie = 100
"""
def coord_to_index(coord):
    return coord[0] + coord[1]*size[0]

def index_to_coord(index):
    return (index % size[0], index // size[0])

def unique_edge(left_node_index, right_node_index):
    return (min(left_node_index, right_node_index), max(left_node_index, right_node_index))

start_point_index = coord_to_index(start_point)
end_point_index = coord_to_index(end_point)
"""
def plot_graph(graph, filename="graph_plot2.png", path=[]):
    pos = nx.spring_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_color='grey')
    nx.draw_networkx_nodes(graph, pos, nodelist=[node for node in graph.nodes if (node == start_point)], node_color='red')
    nx.draw_networkx_nodes(graph, pos, nodelist=[node for node in graph.nodes if (node == end_point)], node_color='green')
    nx.draw_networkx_edges(graph, pos)
    if len(path) > 0:
        nx.draw_networkx_edges(graph, pos, edgelist=path, edge_color='r')
    plt.show()

G1 = nx.Graph()
G1.add_nodes_from([(i, j) for i in range (0, size[0]) for j in range(0, size[1])])
for i in range(0, size[0]-1):
    G1.add_edges_from([((i, j), (i+1, j)) for j in range(0, size[1])])
    #G1.add_edges_from([((i, j), (i+2, j)) for j in range(0, size[1])])
for j in range(0, size[1]-1):
    G1.add_edges_from([((i, j), (i, j+1)) for i in range(0, size[0])]) 
    #G1.add_edges_from([((i, j), (i, j+2)) for i in range(0, size[0])]) 

#plot_graph(G1)
Q = defaultdict(int)
for j in range(0, size[1]):
    for i in range(0, size[0]):
        while((i,j) != (size[0], size[1])):
            S = [(i,j)]
print (S)
entring_edge_t =[(node, end_point) for node in G1.adj[end_point]]
entring_edge_s =[(node, start_point) for node in G1.adj[start_point]]
exiting_edge_t =[(end_point, node) for node in G1.adj[end_point]]
exiting_edge_s=[(start_point, node) for node in G1.adj[start_point]]

#CONTRAINTE (1) min(SUM(Xij))
for edge in G1.edges:
    Q[(edge,edge)] += 1 * lagrangian_objective

#CONTRAINTE (2) SUM(Xsj)-SUM(Xjs)=1    
for edge in exiting_edge_s:
    Q[(edge,edge)] -= 1 *lagrangian_extrem
    for edge2 in exiting_edge_s:
        if (edge != edge2):
            Q[(edge, edge2)] += 2 *lagrangian_extrem
    for edge2 in entring_edge_s:
        Q[edge, edge2] -= 2 *lagrangian_extrem
for edge in entring_edge_s:
    Q[(edge,edge)] += 3 *lagrangian_extrem
    for edge2 in entring_edge_s:
        if (edge != edge2):
            Q[(edge, edge2)] += 2 *lagrangian_extrem  
    


#CONTRAINTE (2) SUM(Xtj)-SUM(Xjt)=-1     
for edge in exiting_edge_t:
    Q[(edge,edge)] += 3 *lagrangian_extrem
    for edge2 in exiting_edge_t:
        if (edge != edge2):
            Q[(edge, edge2)] += 2 *lagrangian_extrem
    for edge2 in entring_edge_t:
        Q[(edge, edge2)] -= 2 *lagrangian_extrem
for edge in entring_edge_t:
    Q[(edge,edge)] -= 1 *lagrangian_extrem
    for edge2 in entring_edge_t:
        if (edge != edge2):
            Q[(edge, edge2)] += 2 *lagrangian_extrem  

for node in G1.nodes:
    
    if ((node != end_point) & (node != start_point)):
        #print("nodes :", node)
        entring_edge = [(enode, node) for enode in G1.adj[node]] 
        exiting_edge = [(node, enode) for enode in G1.adj[node]] 
        

        #CONTRAINTE (3) SUM(Xij) <=1
        for edge in exiting_edge:
            for edge2 in exiting_edge:
                if (edge != edge2):
                    Q[(edge,edge2)] +=1 * lagrangien_sortie
        #CONTRAINTE (2) SUM(Xij)-SUM(Xji)=0
        for edge in exiting_edge:
            Q[(edge,edge)] += 1 * lagrangian_struct_path
            for edge2 in exiting_edge:
                if (edge != edge2):
                    Q[(edge,edge2)] += 2 * lagrangian_struct_path
            for edge2 in entring_edge:
                    Q[(edge, edge2)] -= 2 * lagrangian_struct_path
        
        for edge in entring_edge:
            Q[(edge,edge)] += 1 * lagrangian_struct_path
            for edge2 in entring_edge:
                if (edge != edge2):
                    Q[(edge, edge2)] += 2 * lagrangian_struct_path


"""
for edge in G1.edges:
    for edge1 in G1.edges:
        print((edge,edge1), " = ", Q[(edge,edge1)])

print(Q)
"""
type(G1.edges)
solver = greedy.SteepestDescentSolver()
sampleset = sampler.sample_qubo(Q,
                               chain_strength=chainstrength,
                               num_reads=numruns)
print(sampleset)

# On affiche les trois premiers chemins
#print(sampleset.variables)
for i in range(0, 5):
    sample = sampleset.record.sample[i]
    print(sample)
    path = [sampleset.variables[i] for i in range(0, len(sample)) if sample[i] == 1]
    print("path: " + str(path))
    plot_graph(G1, path=path)
 