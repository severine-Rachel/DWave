import networkx as nx
import dwave_networkx as dnx
import matplotlib.pyplot as plt
from collections import defaultdict
from dwave.system import DWaveSampler, EmbeddingComposite
import dimod
import random

import greedy
import dwave.inspector

#----- set up the graph -----
G = nx.Graph()
G.add_edges_from([(1,12),(1,13),(12,2),(13,3),(2,23),(23,3),(2,24),(3,34),(24,4),(34,4)])

def get_token():
    '''Returns personal access token. Only required if submitting to autograder.'''
    # TODO: Enter your token here
    return 'qGuo-e7acca8272621f4c81e8d1b07bc6c79972d6b92d'

Q = defaultdict(int)
for i,j in G.edges:
    S = random.choice([1,-1])
    Q[(i,i)] += 1 * S
    Q[(i,j)] -= 2 * S
    Q[(j,j)] += 1 * S
sampler = EmbeddingComposite(DWaveSampler(solver={'topology__type__eq': 'chimera'}))
#sampler = DWaveSampler(solver=dict(qpu=True))

sampleset = sampler.sample_qubo(Q,
                               chain_strength=1,
                               num_reads=1,
                               label='Training - test',
                               return_embedding=True)

print("\nEmbedding found:\n", sampleset.info['embedding_context']['embedding'])

print("\nSampleset:")
print(sampleset)

sampler = dimod.ExactSolver()

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