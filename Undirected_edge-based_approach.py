import networkx as nx
import dwave_networkx as dnx
import matplotlib.pyplot as plt
from collections import defaultdict
from dwave.system import DWaveSampler, EmbeddingComposite
import greedy
import dimod
import dwave.inspector

start=0
end=5

chain_strength = 1
num_reads = 10
#----- set up the graph -----
G = nx.Graph()
graph_size = 6
G = nx.gnp_random_graph(graph_size, 0.50)

#G.add_nodes_from([1, 2, 3, 4, 5, 6])
#G.add_edges_from([(1,2), (1,3), (2,3),(2,4), (3,4), (4,5),(4,6), (5,6), (3,6)])
for g in G.edges:
    print(g)

def get_token():
    '''Returns personal access token. Only required if submitting to autograder.'''
    # TODO: Enter your token here
    return 'qGuo-e7acca8272621f4c81e8d1b07bc6c79972d6b92d'

def plot_graph(graph, filename="graph_plot.png", path=[]):
    pos = nx.spring_layout(graph)   
    nx.draw_networkx_edges(graph, pos, edge_color='grey')
    nx.draw_networkx_nodes(graph, pos, node_color='grey')
    if len(path) > 0:
        nx.draw_networkx_edges(graph, pos, edgelist=path, edge_color='g')
    nx.draw_networkx_nodes(graph, pos, nodelist=[node for node in graph.nodes if (node == start)], node_color='r') 
    nx.draw_networkx_nodes(graph, pos, nodelist=[node for node in graph.nodes if (node == end)], node_color='c')
    plt.show()
"""

Q1=defaultdict(int)
for node in G.nodes:
    Q1[node,node] += 1
    for anode in G.adj[node]:
        Q1[anode,anode] -= 1
        Q1[node,node] -= 1
for edge in G.edges:
    Q1[edge] += 2

sampler = EmbeddingComposite(DWaveSampler(solver={'topology__type__eq': 'pegasus'}))
sampleset = sampler.sample_qubo(Q1,
                               chain_strength=chain_strength,
                               num_reads=num_reads,
                               label='Training - Minpath',
                               return_embedding=True)
print(sampleset)
E = [sampleset.variables[i] for i in range (0, len(sampleset.record.sample[0])) if sampleset.record.sample[0][i] ==  1]
"""
Q = defaultdict(int)

for node in G.nodes:
    Q[(node,node)] += 0.5
    if node == start:
        for edgeS in G.edges(node):
            Q[(node, edgeS)] -= 2
            Q[(edgeS, edgeS)] += 1
            for edgeSadj in G.edges(node):
                if (edgeSadj != edgeS):
                    Q[(edgeS, edgeSadj)] += 1
    elif node == end:
        for edgeT in G.edges(node):
            Q[(node, (min(edgeT),max(edgeT)))] -= 2
            Q[((min(edgeT),max(edgeT)), (min(edgeT),max(edgeT)))] += 1
            for edgeTadj in G.edges(node):
                if ( (min(edgeTadj),max(edgeTadj)) != (min(edgeT),max(edgeT))):
                    Q[((min(edgeT),max(edgeT)),  (min(edgeTadj),max(edgeTadj)))] += 1
    else:
        Q[(node,node)] += 4
        for edgeI in G.edges(node):
            Q[(node, (min(edgeI),max(edgeI)))] -= 4
            Q[((min(edgeI),max(edgeI)), (min(edgeI),max(edgeI)))] += 1
            for edgeIadj in G.edges(node):
                if (min(edgeI),max(edgeI)) != (min(edgeIadj),max(edgeIadj)):
                    Q[((min(edgeI),max(edgeI)), (min(edgeIadj),max(edgeIadj)))] += 1
"""
for node in G.nodes:
    if node == 1:
        for edgeS in G.edges(node):
            print("Q[(",node,",", edgeS,")]",Q[(node, edgeS)])
            print("Q[(",edgeS,",", edgeS,")]",Q[(edgeS, edgeS)])
            for edgeSadj in G.edges(node):
                if (edgeSadj != edgeS):
                    print("Q[(",edgeS,",", edgeSadj,")]",Q[(edgeS, edgeSadj)])
    elif node == 4:
        for edgeT in G.edges(node):
            print("Q[(",node,",", (min(edgeT),max(edgeT)),")]",Q[(node, (min(edgeT),max(edgeT)))])
            print("Q[(",(min(edgeT),max(edgeT)),",", (min(edgeT),max(edgeT)),")]",Q[((min(edgeT),max(edgeT)),(min(edgeT),max(edgeT)))])
            for edgeTadj in G.edges(node):
                if ((min(edgeTadj),max(edgeTadj)) != (min(edgeT),max(edgeT))):
                    print("Q[(",(min(edgeT),max(edgeT)),",", (min(edgeTadj),max(edgeTadj)),")]",Q[((min(edgeT),max(edgeT)), (min(edgeTadj),max(edgeTadj)))])
    else:
        print("Q[(",node,",", node,")]",Q[(node, node)])
        for edgeI in G.edges(node):
            print("Q[(",node,",", (min(edgeI),max(edgeI)),")]",Q[(node, (min(edgeI),max(edgeI)))])
            print("Q[(",(min(edgeI),max(edgeI)),",", (min(edgeI),max(edgeI)),")]",Q[((min(edgeI),max(edgeI)),(min(edgeI),max(edgeI)))])
            for edgeIadj in G.edges(node):
                if (min(edgeI),max(edgeI)) != (min(edgeIadj),max(edgeIadj)):
                    print("Q[(",(min(edgeI),max(edgeI)),",", (min(edgeIadj),max(edgeIadj)),")]",Q[((min(edgeI),max(edgeI)), (min(edgeIadj),max(edgeIadj)))])
    """


sampler = EmbeddingComposite(DWaveSampler(solver={'topology__type__eq': 'pegasus'}))
sampleset = sampler.sample_qubo(Q,
                               chain_strength=chain_strength,
                               num_reads=num_reads,
                               label='Training - Minpath',
                               return_embedding=True)

print("\nEmbedding found:\n", sampleset.info['embedding_context']['embedding'])

print("\nSampleset:")
print(sampleset)

sampler = dimod.ExactSolver()
path = [index for index in G.edges if sampleset.first.sample[index] == 1]
print("path :\t",path)
# ----- Print results to user -----
print("\nSolutions:")
print('-' * 60)
print('{:>30s}{:>30s}{:^15s}{:^15s}'.format('Set 0','Set 1','Energy','Cut Size'))
print('-' * 60)
for sample, E in sampleset.data(fields=['sample','energy']):
    S0 = [k for k,v in sample.items() if v == 0]
    S1 = [k for k,v in sample.items() if v == 1]
    print('{:>15s}{:>15s}{:^15s}{:^15s}'.format(str(S0),str(S1),str(E),str(int(-1*E))))

dwave.inspector.show(sampleset)

plot_graph(G, path=path)

filename = "edge_based_approach.png"
plt.savefig(filename, bbox_inches='tight')
print("\nYour plot is saved to {}".format(filename))