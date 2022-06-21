import networkx as nx
import dwave_networkx as dnx
import matplotlib.pyplot as plt
from collections import defaultdict
from dwave.system import DWaveSampler, EmbeddingComposite
import greedy
import dimod
import dwave.inspector

start=1
end=5

#----- set up the graph -----
def plot_graph(graph, filename="graph_plot.png", path=[]):
    pos = nx.spring_layout(graph)   
    nx.draw_networkx_nodes(graph, pos, node_color='grey')
    if len(path) > 0:
        nx.draw_networkx_nodes(graph, pos, nodelist=path, node_color='g')
    nx.draw_networkx_edges(graph, pos)
    nx.draw_networkx_nodes(graph, pos, nodelist=[node for node in graph.nodes if (node == start)], node_color='r') 
    nx.draw_networkx_nodes(graph, pos, nodelist=[node for node in graph.nodes if (node == end)], node_color='c')
    plt.show()


G = nx.Graph()
#G.add_nodes_from([1, 12, 13, 2, 3, 23, 24, 34, 4])
#G.add_edges_from([(1,12),(1,13),(12,2),(13,3),(2,23),(23,3),(2,24),(3,34),(24,4),(34,4)])

G.add_nodes_from([1, 12, 13, 2, 3, 23, 24, 34, 4, 45, 5, 6, 46, 56, 36])
G.add_edges_from([(1,12),(1,13),(12,2),(13,3),(2,23),(23,3),(2,24),(3,34),(24,4),(34,4), (4,45), (45,5),(4,46),(46,6),(5,56),(56,6),(3,36),(36,6)])

#graph_size = 10
#G = nx.gnp_random_graph(graph_size, 0.30)
def get_token():
    '''Returns personal access token. Only required if submitting to autograder.'''
    # TODO: Enter your token here
    return 'qGuo-e7acca8272621f4c81e8d1b07bc6c79972d6b92d'

sampler = EmbeddingComposite(DWaveSampler(solver={'topology__type__eq': 'pegasus'}))
#sampler = DWaveSampler(solver=dict(qpu=True))
S = dnx.min_vertex_cover(G, sampler, lagrange=2.0)
print (S)
D=2


print (" D = ", D)

chain_strength = 2
num_reads = 100
"""
Q1=defaultdict(int)
for node in G.nodes:
    Q1[node,node] -= 1
    for anode in G.adj[node]:
        Q1[anode,anode] -= 0.5
        Q1[node,node] -= 0.5
for edge in G.edges:
    Q1[edge] += 2

sampler = EmbeddingComposite(DWaveSampler(solver={'topology__type__eq': 'pegasus'}))
sampleset = sampler.sample_qubo(Q1,
                               chain_strength=chain_strength,
                               num_reads=num_reads,
                               label='Training - Minpath',
                               return_embedding=True)

E = [sampleset.variables[i] for i in range (0, len(sampleset.record.sample[0])) if sampleset.record.sample[0][i] ==  1]
"""
#sampler = dimod.ExactSolver()

#----- set up QUBO dict -----

Q = defaultdict(int)

for i in S:
    Q[(i,i)] += 1
    if i == start:
        for adjs in G.adj[i]:
            Q[(adjs, adjs)] += 1*D
            Q[(i, adjs)] -= 2*D
            for adjs1 in G.adj[i]:
                if adjs != adjs1:
                    Q[(adjs, adjs1)] += 1*D
    elif i == end:
        for adjt in G.adj[i]:
            Q[(adjt, adjt)] += 1*D
            Q[(i, adjt)] -= 2*D
            for adjt1 in G.adj[i]:
                if adjt != adjt1:
                    Q[(adjt, adjt1)] += 1*D
    else:
        Q[(i,i)] += 4*D
        for adji in G.adj[i]:
            Q[(i,adji)] -= 4*D
            Q[(adji,adji)] += 1*D
            for adji1 in G.adj[i]:
                if adji != adji1:
                    Q[(adji,adji1)] +=1*D
"""
for i in S:
    if i == start:
        for adjs in G.adj[i]:
            print("Q[(",adjs,",", adjs,")] =\t", Q[(adjs, adjs)])
            print("Q[(",i,",", adjs,")] =\t", Q[(i, adjs)])
            for adjs1 in G.adj[i]:
                if adjs != adjs1:
                    print("Q[(",adjs,",", adjs1,")] =\t", Q[(adjs, adjs1)])
    elif i == end:
        for adjt in G.adj[i]:
            print("Q[(",adjt,",", adjt,")] =\t", Q[(adjt, adjt)])
            print("Q[(",i,",", adjt,")] =\t", Q[(i, adjt)])
            for adjt1 in G.adj[i]:
                if adjt != adjt1:
                    print("Q[(",adjt,",", adjt1,")] =\t", Q[(adjt, adjt1)])
    else:
        print("Q[(",i,",", i,")] =\t", Q[(i, i)])
        for adji in G.adj[i]:
            print("Q[(",i,",", adji,")] =\t", Q[(i, adji)])
            print("Q[(",adji,",", adji,")] =\t", Q[(adji, adji)])
            for adji1 in G.adj[i]:
                if adji != adji1:
                    print("Q[(",adji,",", adji1,")] =\t", Q[(adji, adji1)])
"""
sampler = EmbeddingComposite(DWaveSampler(solver={'topology__type__eq': 'pegasus'}))
sampleset = sampler.sample_qubo(Q,
                               chain_strength=chain_strength,
                               num_reads=num_reads,
                               label='Training - Minpath',
                               return_embedding=True)
print(sampleset.info["timing"])  
#print("\nEmbedding found:\n", sampleset.info['embedding_context']['embedding'])

#print("\nSampleset:")
#print(sampleset)
print(sampleset.first)
sampler = dimod.ExactSolver()
path = [index for index in G if sampleset.first.sample[index] == 1]
print(path)
"""
for i in range(0, 3):
    sample = sampleset.record.sample[i]
    print(sample)
    path = [sampleset.variables[i] for i in range(0, len(sample)) if sample[i] == 1]
    print("path: " + str(path))

# ----- Print results to user -----
print("\nSolutions:")
print('-' * 60)
print('{:>30s}{:>30s}{:^15s}{:^15s}'.format('Set 0','Set 1','Energy','Cut Size'))
print('-' * 60)
for sample, E in sampleset.data(fields=['sample','energy']):
    S0 = [k for k,v in sample.items() if v == 0]
    S1 = [k for k,v in sample.items() if v == 1]
    print('{:>30s}{:>30s}{:^15s}{:^15s}'.format(str(S0),str(S1),str(E),str(int(-1*E))))
"""

dwave.inspector.show(sampleset)

plot_graph(G, path=path)

filename = "node_based_approach.png"
plt.savefig(filename, bbox_inches='tight')
print("\nYour plot is saved to {}".format(filename))