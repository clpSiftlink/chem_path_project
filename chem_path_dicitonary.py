import pandas as pd
import numpy as np
import networkx as nx
from scipy import sparse
import pickle


dif_matrix = sparse.load_npz('dif_matrix_all_chems_fixed.npz')

net_01_df = pd.read_csv('/home/shared/projects/siftlink/nclc/net01_uni.csv', header = 0)
prot_path_df = pd.read_csv('/home/shared/projects/siftlink/nclc/protein_pathways_HSA.tsv', sep = '\t' ,header = None)
ppi_df = pd.read_csv('/home/clp/csvfiles/stack_df_fixed.csv')



net_01_df = net_01_df.dropna()
prot_path_df = prot_path_df.dropna()

ppi_net = nx.from_pandas_edgelist(ppi_df, source='source', target='target', edge_attr='probs')

nodes = list(ppi_net.nodes)
snodes = set(nodes)

def protein_protein_score(dif_matrix, nodes, p1, p2):
    return dif_matrix[nodes.index(p1), nodes.index(p2)]

def chem_path_proteins_scores(dif_matrix, snodes, protein_to_num_map, chem_proteins, pathway_proteins):
    """
    dif_matrix: diffusion matrix
    snodes: set of nodes in ppi_network
    protein_to_num_map: a dictionary that maps the proteins to the indices of the diffusion matrix
    chem_proteins: all the proteins connected to the specific chemical we are interested in
    pathway_proteins: all the proteins connected to the specific pathway we are interested in
    """
    
    chem_proteins_num = []
    pathway_proteins_num = []
    
    for cprot in chem_proteins:
        if(cprot in snodes):
            chem_proteins_num.append(protein_to_num_map[cprot])
        
    for pprot in pathway_proteins:
        if(pprot in snodes):
            pathway_proteins_num.append(protein_to_num_map[pprot])
        
    return dif_matrix.tocsr()[chem_proteins_num, :].tocsc()[:,pathway_proteins_num].sum()

chems = net_01_df.chem.unique()
pathways = prot_path_df[1].unique()

nodes_num = [i for i in range(len(nodes))]
zip_iterator = zip(nodes, nodes_num)
dic_nodes = dict(zip_iterator)

chem_proteins = {}
for chem in chems:
    chem_proteins[chem] = list(net_01_df.loc[net_01_df.chem==chem].uniprot.values)
pathway_proteins = {}
for pathway in pathways:
    pathway_proteins[pathway] = list(prot_path_df.loc[prot_path_df[1]==pathway][0].values)


def path_chem_score_dic(dif_matrix, nodes, dic_nodes, chems, pathways, chem_proteins, pathway_proteins):
    dic1 = {}
    snodes = set(nodes)
    counter = 0
    for chem in chems:
        for pathway in pathways:
            counter +=1
            if(counter%100==0):
                print(counter)
            tupple = (chem,pathway)
            dic1[tupple] = chem_path_proteins_scores(dif_matrix, snodes, dic_nodes, chem_proteins[chem], pathway_proteins[pathway])
    return dic1


chem_path_dictionary = path_chem_score_dic(dif_matrix, nodes, dic_nodes, chems, pathways, chem_proteins, pathway_proteins)

file_to_store = open("chem_path_dictionary_fixed.pickle", "wb")
pickle.dump(chem_path_dictionary, file_to_store)
file_to_store.close()