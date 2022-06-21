import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict
import greedy
import dimod
import neal
from dimod import CQM

# Size of the grid
size = (3,3)

# DiGraph of the problem
G = nx.DiGraph()
for i in range(0, size[0]):
    for j in range(0, size[1]):
        G.add_node((i, j))

for i in range(0, size[0]):
    for j in range(0, size[1]):
        if i < size[0]-1:
            G.add_edge((i, j), (i+1, j))
            G.add_edge((i+1, j), (i, j))
        if j < size[1]-1:
            G.add_edge((i, j), (i, j+1))
            G.add_edge((i, j+1), (i, j))

##pos = nx.spring_layout(G)
##nx.draw_networkx_nodes(G, pos, node_color='grey')
##nx.draw_networkx_edges(G, pos)
##plt.show()

## Creation de la CQM
# On créer une variable binaire par arc du graphe
e = [dimod.Binary(f'{i}') for i in G.edges]


cqm = dimod.ConstrainedQuadraticModel()
cqm.set_objective(sum(e))

start_point = (0,0)
end_point = (2,2)

# Contrainte de structure
for node in G.nodes:
    # Traitement du point de départ
    if node == start_point:
        start_out_edges = [dimod.Binary(f'{i}') for i in G.out_edges(start_point)]
        start_in_edges = [dimod.Binary(f'{i}') for i in G.in_edges(start_point)]
        #print(start_out_edges)
        cqm.add_constraint(sum(start_out_edges) == 1)
        cqm.add_constraint(sum(start_in_edges) == 0)
    # Traitement du point d'arrivé
    elif node == end_point:
        end_in_edges = [dimod.Binary(f'{i}') for i in G.in_edges(end_point)]
        end_out_edges = [dimod.Binary(f'{i}') for i in G.in_edges(end_point)]
        #print(end_in_edges)
        cqm.add_constraint(sum(end_in_edges) == 1)
        cqm.add_constraint(sum(end_out_edges) == 0)
    # Traitement des autres points
    else:
        in_edges = [dimod.Binary(f'{i}') for i in G.in_edges(node)]
        print(G.in_edges(node))
        out_edges = [dimod.Binary(f'{i}') for i in G.out_edges(node)]
        cqm.add_constraint(sum(in_edges) - sum(out_edges) == 0)

    # Nombre maximum d'arc sortant <= 1
    in_edges = [dimod.Binary(f'{i}') for i in G.in_edges(node)]
    cqm.add_constraint(sum(in_edges) <= 1)

bqm, invert = dimod.cqm_to_bqm(cqm)
qubo, c = bqm.to_qubo()
print(qubo)
solver = greedy.SteepestDescentSolver()
#solver = neal.SimulatedAnnealingSampler()
sampleset = solver.sample_qubo(qubo, num_reads=5)
print(sampleset)

# Traitement des résultats
# TODO traiter les résultats invalides

pos = nx.spring_layout(G)
nx.draw_networkx_nodes(G, pos, node_color='grey')
nx.draw_networkx_nodes(G, pos, [start_point], node_color='green')
nx.draw_networkx_nodes(G, pos, [end_point], node_color='red')
#nx.draw_networkx_edges(G, pos)
sample = sampleset.record.sample[0]
binary_path = [sampleset.variables[i] for i in range(0, len(sample)) if sample[i] == 1]
edges_to_draw = [edge for edge in G.edges if f'{edge}' in binary_path]
nx.draw_networkx_edges(G, pos, edges_to_draw, edge_color='green')
plt.show()

# Coût d'un chemins
#print(len(G.edges))
#for edge in G.edges:
#    print(edge)