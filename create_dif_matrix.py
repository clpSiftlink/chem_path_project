import pandas as pd
import numpy as np
import networkx as nx
from scipy import sparse
import pickle
import math

ppi_df = pd.read_csv('/home/clp/nclc/ppi_uni.csv')
ppi_df = ppi_df[['uni_protein1', 'uni_protein2', 'combined_score']]
ppi_df = ppi_df.dropna()
ppi_df.combined_score = ppi_df.combined_score/1000


ppi_net = nx.from_pandas_edgelist(ppi_df, source='uni_protein1', target='uni_protein2', edge_attr='combined_score')


def ap_diff_matrix_approximation(network, power, lamda):
    
    lapl =  nx.linalg.laplacianmatrix.normalized_laplacian_matrix(ppi_net, weight='probs')
    nonzero_mask = np.array(abs(lapl[lapl.nonzero()]) < 0.006)[0]
    rows = lapl.nonzero()[0][nonzero_mask]
    cols = lapl.nonzero()[1][nonzero_mask]
    lapl[rows, cols] = 0
    
    
    final_lapl = lapl.copy()
    tmp_lapl = lapl.copy()
    
    for i in range(1,power):
        tmp_lapl = sparse.csr_matrix(tmp_lapl).dot(sparse.csr_matrix(lapl.T))
        final_lapl += (math.pow(lamda,i+1)/math.factorial(i+1)) * tmp_lapl
        nonzero_mask = np.array(abs(tmp_lapl[tmp_lapl.nonzero()]) < 0.006)[0]
        rows = tmp_lapl.nonzero()[0][nonzero_mask]
        cols = tmp_lapl.nonzero()[1][nonzero_mask]
        tmp_lapl[rows, cols] = 0
        print(i)
        
    return final_lapl


def ap_diff_matrix(network, power, lamda):
    
    lapl =  nx.linalg.laplacianmatrix.normalized_laplacian_matrix(ppi_net, weight='prob')
 
    
    final_lapl = lapl.copy()
    tmp_lapl = lapl.copy()
    
    for i in range(1,power):
        tmp_lapl = sparse.csr_matrix(tmp_lapl).dot(sparse.csr_matrix(lapl.T))
        final_lapl += (math.pow(lamda,i+1)/math.factorial(i+1)) * tmp_lapl
        
    return final_lapl


power = 10
lamda = 0.2
dif_matrix = ap_diff_matrix(ppi_net, power, lamda)

file_to_store = open("ap_diff_matrix_total.pickle", "wb")
pickle.dump(dif_matrix, file_to_store)
file_to_store.close()