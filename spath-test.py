import networkx as nx
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

#convertion tuple coordonné en int index
def coord_to_index(coord):
    return coord[0] + coord[1]*size[0]

#convertion int index en tuple coordonné
def index_to_coord(index):
    return (index % size[0], index // size[0])

#affiche coordonné avec la position minimal sur x et maximal sur y, ça donne juste un graphe en longueur
def unique_edge(left_node_index, right_node_index): #prend un tuple en entré
    return (min(left_node_index, right_node_index), max(left_node_index, right_node_index)) #retourne un tuple

#n'est pas utilisé, remet le bord inférieur toujours a gauche de la liaison de G2 (composé de deux tuples)
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

#convertion coord en index vrai que si x : min et y : max
assert (unique_edge(3,2)) == index_to_coord(coord_to_index(unique_edge(3,2)))

# Start and end node
start_point = (1,1)
end_point = (2,3)

start_point_index = coord_to_index(start_point) #4
end_point_index = coord_to_index(end_point)#11

for j in range(0, size[1]):      #for(int j = 0, j < 4, i++)
    for i in range(0, size[0]):  #for(int i = 0, i < 3, i++)
        if (i, j) == start_point:
            print("S", end="")
        elif (i, j) == end_point:
            print("E", end="")
        else:
            print(".", end="")
    print("")

#fonction appliquée a G2
def plot_graph(graph, filename="graph_plot.png", path=[]):
    pos = nx.spring_layout(graph)   
    nx.draw_networkx_nodes(graph, pos, node_color='grey')   #color tous les points en gris
    #ne sert a rien car seul les points de G1 sont caractérisés par un index
    #nx.draw_networkx_nodes(graph, pos, nodelist=[node for node in graph.nodes if node == coord_to_index(start_point[0], start_point[1])], node_color='red')
    #nx.draw_networkx_nodes(graph, pos, nodelist=[node for node in graph.nodes if node == coord_to_index(end_point[0], end_point[1])], node_color='black')

    #dessine en rouge la liste de point tel que la coordonné x du point ou y soit egale au coordonné du point start
    # start_point_index est le point de G1 qui commence, relié par 4 liaison devenu des points dans G2, il y a donc 4 points ayant la coordonnée de départ
    nx.draw_networkx_nodes(graph, pos, nodelist=[node for node in G2.nodes if node[0] == start_point_index or node[1] == start_point_index], node_color='r') 
    #dessine en cyan la liste de point tel que la coordonné x du point ou y soit egale au coordonné du point end   
    #de la même manière, il y a deux points ayant les coordonnées d'arrivé. 
    nx.draw_networkx_nodes(graph, pos, nodelist=[node for node in G2.nodes if node[0] == end_point_index or node[1] == end_point_index], node_color='c')

    if len(path) > 0:
        nx.draw_networkx_nodes(graph, pos, nodelist=path, node_color='g') #tous les points du chemin sont colorié en vert
    nx.draw_networkx_edges(graph, pos) #dessine les liaisons
    
    #nx.draw_networkx_labels(graph, pos)

    plt.show()  #affiche le graphe
    #plt.savefig(filename)

# Initial graph
G1 = nx.Graph()
#les points de G1 ne servent a rien dans le programme
G1.add_nodes_from([(coord_to_index((i, j)), {"x":i,"y":j}) for i in range (0, size[0]) for j in range(0, size[1])])

#tous les points sont reliés sur l'axe x
for i in range(0, size[0]-1):
    G1.add_edges_from([unique_edge(coord_to_index((i, j)), coord_to_index((i+1, j))) for j in range(0, size[1])]) 

#tous les points sont reliés sur l'axe y
for j in range(0, size[1]-1):
    G1.add_edges_from([unique_edge(coord_to_index((i, j)), coord_to_index((i, j+1))) for i in range(0, size[0])]) 


G2 = nx.Graph()
#les chemins de G1 deviennent les points de G2
G2.add_nodes_from([unique_edge(u, v) for (u, v) in G1.edges])
#si deux coté de G1 sont relié par un point, alors on créé une liaison de G2 reliant les deux liaisons de G1 (devenu des points de G2)
for (u, v) in G1.edges:
    for (u_i, v_i) in G1.edges:
        if ((u in [u_i, v_i] or v in [u_i, v_i])):  #concordance d'une extremité entre deux liaisons
            G2.add_edge(unique_edge(u, v), unique_edge(u_i, v_i))

plot_graph(G2)

Q = defaultdict(int)
nodes_from_start = [node for node in G2.nodes if node[0] == start_point_index or node[1] == start_point_index]
nodes_from_end = [node for node in G2.nodes if node[0] == end_point_index or node[1] == end_point_index]

# QUBO trouvée : SUM(xii) - SUM(x)

# ajout des noeuds dans le graphe
print('-' * 60 + "TOUTES LES VALEURS BINAIRES DU GRAPHE" + '-' * 60)
#Constraint (1)
for node in G2.nodes:
    Q[(node, node)] = 1 * lagrangian_objective
    print("les noeuds ", (node, node), " ont comme valeur ", Q[(node, node)])

# contrainte sur le départ
print('-' * 60 + "TOUTES LES VALEURS DE DEPART DU GRAPHE" + '-' * 60)
# CONSTRAINT (2)
for node in nodes_from_start:
    Q[(node, node)] += -1 * lagrangian_struct_start
    print("les noeuds ", (node, node), " ont comme valeur ", Q[(node, node)])
    for node2 in nodes_from_start:
        if (node != node2):
            Q[(node, node2)] += 1 * lagrangian_struct_start
            print("les noeuds ", (node, node2), " ont comme valeur ", Q[(node, node2)])

# contrainte sur l'arrivée
print('-' * 60 + "TOUTES LES VALEURS D'ARRIVE DU GRAPHE" + '-' * 60)
# CONSTRAINT (2)
for node in nodes_from_end:
    Q[(node, node)] += -1 * lagrangian_struct_end
    print("les noeuds ", (node, node), " ont comme valeur ", Q[(node, node)])
    for node2 in nodes_from_end:
        if (node != node2):
            Q[(node, node2)] += 1 * lagrangian_struct_end
            print("les noeuds ", (node, node2), " ont comme valeur ", Q[(node, node2)])

print('-' * 60 + "TOUS LES POINTS ENTRANTS" + '-' * 60)
for node in G2.nodes:
    # arcs dits entrants
    #point voisin dont leur coordonnée partage une coordonnée entrante node[0] du point G2 qui correspond a l'extrémité de la liaison de G1 par lequel on entre
    entring_node = [enode for enode in G2.adj[node] if enode[0] == node[0] or enode[1] == node[0]]
    print ("Les points entrants de ", node," sont ", entring_node)

print('-' * 60 + "TOUS LES POINTS SORTANTS" + '-' * 60)
for node in G2.nodes:
    #point voisin dont leur coordonnée partage une coordonnée sortante node[0] du point G2 qui correspond au point de la liaison de G1 sortant
    exiting_node = [enode for enode in G2.adj[node] if enode[0] == node[1] or enode[1] == node[1]]
    print ("Les points sorants de ", node," sont ", exiting_node)

    ## contrainte sur les arcs de passage, représentés par des noeuds dans G2
    # On vérifie que le nombre total de voisins est bien égal à la somme des sortants et entrants
    #assert G2.degree[node] == len(entring_node) + len(exiting_node)

print('-' * 60 + "VALEURS DES POINTS ENTRANTS" + '-' * 60)    
for node in G2.nodes:
    # Maintenant on définit les contraintes : 
    # Sum_{entrants}(xi) + Sum_{sortants}(xi) + 2*Sum_{entrants,j>i}(xi*xj) + 2*Sum_{sortants,j>i}(xi*xj) - 2*Sum_{entrants}(xi)*Sum_{sortants}(xi)
    # Somme des xi entrants 
    # CONSTRAINT (2)
    for enode in entring_node:
        Q[(enode, enode)] += 1 * lagrangian_struct_path
for enode in entring_node:
    print("valeur points entrant :", (enode, enode), Q[(enode, enode)], "pour la valeur des noeurs entrants : ", entring_node)

print('-' * 60 + "VALEURS DES POINTS SORTANTS" + '-' * 60)    
for node in G2.nodes:
    # Somme des xi sortants 
    # CONSTRAINT (2)
    for enode in exiting_node:
        Q[(enode, enode)] += 1 * lagrangian_struct_path
for enode in exiting_node:
    print("valeur points sortants :", (enode, enode), Q[(enode, enode)])

print('-' * 60 + "VALEURS DES ARCS ENTRANTS" + '-' * 60)         
for node in G2.nodes:
    # 2*Sum_{entrants,j>i}(xi*xj) 
    # CONSTRAINT (2)
    for enode in entring_node:
        for enode2 in entring_node: # Avec la double boucle on compte deux fois, donc on divise le coeff par 2
            if (enode != enode2):
                Q[(enode, enode2)] += 1 * lagrangian_struct_path
for enode in entring_node:
    for enode2 in entring_node: 
        if (enode != enode2):
            print("valeur arcs entrant :", (enode, enode2), Q[(enode, enode2)])
print('-' * 60 + "VALEURS DES ARCS SORTANTS" + '-' * 60)  
for node in G2.nodes:
    # 2*Sum_{sortants,j>i}(xi*xj) 
    # CONSTRAINT (2)
    for enode in exiting_node:
        for enode2 in exiting_node: # Avec la double boucle on compte deux fois, donc on divise le coeff par 2
            if (enode != enode2):
                Q[(enode, enode2)] += 1 * lagrangian_struct_path
for enode in exiting_node:
        for enode2 in exiting_node: # Avec la double boucle on compte deux fois, donc on divise le coeff par 2
            if (enode != enode2):
                print("valeur arcs sortant :", (enode, enode2), Q[(enode, enode2)])
print('-' * 60 + "VALEURS DES ARCS ENTRANTS ET SORTANTS D'UN MÊME POINT" + '-' * 60) 
for node in G2.nodes:
    # -2*Sum_{entrants}(xi)*Sum_{sortants}(xi) 
    # CONSTRAINT (2)
    for enode in entring_node:
        for enode2 in exiting_node:
            Q[(enode, enode2)] += -2 * lagrangian_struct_path
for enode in entring_node:
        for enode2 in exiting_node:
            print("valeur des points entrants et sortants d'une liaison :", (enode, enode2), Q[(enode, enode2)])

for node in G2.nodes:
    ## Contrainte sur le nombre d'arcs sortants
    for i in range(0, len(exiting_node)-1):
        for j in range(i+1, len(exiting_node)):
            Q[(exiting_node[i], exiting_node[j])] += 1 * lagrangian_struct_path

#MANQUE CONSTRAINT (3)

for i in range(0, len(exiting_node)-1):
        for j in range(i+1, len(exiting_node)):
            print (Q[(exiting_node[i], exiting_node[j])])

print(Q)
#plot_graph(G2)

solver = greedy.SteepestDescentSolver()
sampleset = solver.sample_qubo(Q, num_reads=100)

print(sampleset)

# On affiche les trois premiers chemins
print(sampleset.variables)
for i in range(0, 3):
    sample = sampleset.record.sample[i]
    print("sample:")
    print(sample)
    path = [sampleset.variables[i] for i in range(0, len(sample)) if sample[i] == 1]
    print("path: " + str(path))
    plot_graph(G2, path=path)

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