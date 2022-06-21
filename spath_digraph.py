import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict
import greedy
import dimod
import neal
from dimod import CQM

# Size of the grid
size = (3,4)

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

pos = nx.spring_layout(G)
nx.draw_networkx_nodes(G, pos, node_color='grey')
nx.draw_networkx_edges(G, pos)
#plt.show()


## Creation de la CQM
# On créer une variable binaire par arc du graphe
e = [dimod.Binary(f'{i}') for i in G.edges]
cqm = dimod.ConstrainedQuadraticModel()
cqm.set_objective(sum(e))

bqm, invert = dimod.cqm_to_bqm(cqm)
qubo, c = bqm.to_qubo()
print("creation \n", qubo)

# On impose le point de départ
start_out_edges = [dimod.Binary(f'{i}') for i in G.out_edges((0,0))]
print("start_out_edges \n", start_out_edges)
cqm.add_constraint(sum(start_out_edges) == 1)
print("start \n", qubo)

# On impose le point d'arrivée
end_in_edges = [dimod.Binary(f'{i}') for i in G.in_edges((2,2))]
print("end_in_edges \n", end_in_edges)
cqm.add_constraint(sum(end_in_edges) == 1)

print("end \n", qubo)

bqm, invert = dimod.cqm_to_bqm(cqm)
qubo, c = bqm.to_qubo()
print("end \n", qubo)
print("now \n",qubo)
#solver = greedy.SteepestDescentSolver()
solver = neal.SimulatedAnnealingSampler()
sampleset = solver.sample_qubo(qubo, num_reads=5)
print(sampleset)

# Coût d'un chemins
#print(len(G.edges))
#for edge in G.edges:
#    print(edge)