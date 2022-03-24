import pandas as pd
import numpy as np
import networkx as nx
from scipy import sparse
import pickle

df = pd.read_csv('/shared/dev01_clp_env/csvfiles/chem_path_score_df.csv')
df = df.drop('Unnamed: 0', axis=1)

df.scores = abs(df.scores)
chems = list(df.chems.unique())
print(len(chems))

def get_positive_negative_scores(df, pathway):
    positive_scores = []
    negative_scores = []

    for chem in chems:
        tmp_df = df[df.chems == chem]
        positive_score = tmp_df.loc[tmp_df.pathways == pathway].scores.values[0]
        negative_score = tmp_df.loc[tmp_df.pathways != pathway].scores.sum()
        positive_scores.append(positive_score)
        negative_scores.append(negative_score)
    return np.array(positive_scores), np.array(negative_scores)

def get_chem_score(pathway):
    positive_scores, negative_scores = get_positive_negative_scores(df, pathway)
    scores = positive_scores - negative_scores
    
    pzip_iterator = zip(chems, positive_scores)
    pchem_score = dict(pzip_iterator)
    
    nzip_iterator = zip(chems, negative_scores)
    nchem_score = dict(nzip_iterator)
    
    zip_iterator = zip(chems, scores)
    chem_score = dict(zip_iterator)
    
    
    return chem_score, pchem_score, nchem_score


file_to_read = open('/shared/dev01_clp_env/test/knap_sack_chem_sets_dic.pickle', 'rb')
knap_sack_chem_sets_dic = pickle.load(file_to_read)
file_to_read.close()

file_to_read = open('/shared/dev01_clp_env/test/set_of_chem_sets_dic.pickle', 'rb')
set_of_chem_sets_dic = pickle.load(file_to_read)
file_to_read.close()

file_to_read = open('/shared/dev01_clp_env/test/max_chem_sets_dic.pickle', 'rb')
max_chem_sets_dic = pickle.load(file_to_read)
file_to_read.close()




pathways = list(knap_sack_chem_sets_dic.keys())

pscores_k = []
pscores_c = []
pscores_m = []
nscores_k = []
nscores_c = []
nscores_m = []

counter = 0
for pathway in pathways:
    counter +=1
    print(counter)
    chem_score, pchem_score, nchem_score = get_chem_score(pathway)
    pkscore = 0
    pcscore = 0
    pmscore = 0
    nkscore = 0
    ncscore = 0
    nmscore = 0
    for chem in knap_sack_chem_sets_dic[pathway]:
        pkscore += pchem_score[chem]
        nkscore += nchem_score[chem]
    for chem in set_of_chem_sets_dic[pathway]:
        pcscore += pchem_score[chem]
        ncscore += nchem_score[chem]
    for chem in max_chem_sets_dic[pathway]:
        pmscore += pchem_score[chem]
        nmscore += nchem_score[chem]

    pscores_k.append(pkscore)
    nscores_k.append(nkscore)
    pscores_c.append(pcscore)
    nscores_c.append(ncscore)
    pscores_m.append(pmscore)
    nscores_m.append(nmscore)

final_df = pd.DataFrame({'knapsack_p':pscores_k, 'knapsack_n':nscores_k, 'cutoff_p':pscores_c, 'cutoff_n':nscores_c, 'max_p':pscores_m, 'max_n':nscores_m})
final_df.to_csv('/shared/dev01_clp_env/csvfiles/methods_score_df.csv', index=False)

