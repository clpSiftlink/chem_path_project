from importlib.resources import path
import pandas as pd
import numpy as np


df = pd.read_csv("chem_path_scores_df_fixed.csv")
df.drop('Unnamed: 0', axis=1, inplace=True)
df.scores = abs(df.scores.round(3))
df = df.sort_values(['chems', 'scores'], ascending=False)


pathways = df.pathways.unique()
chems = df.chems.unique()

def find_chems_with_pathway_as_max(df, chems, pathway):
    chems_of_interest = []
    for chem in chems:
        tmp_pway = df.loc[df.chems == chem].loc[df.loc[df.chems == chem].scores == df.loc[df.chems == chem].scores.max()].pathways.iloc[0]
        if tmp_pway == pathway:
            chems_of_interest.append(chem)

    score_plus = 0
    score_minus = 0
    if len(proper_chem_set) > 0:
        for chem in proper_chem_set:
            tmp_df = df[df.chems == chem]
            score_plus += tmp_df.loc[tmp_df.pathways == pathway].scores.values[0]
            score_minus += tmp_df.loc[tmp_df.pathways != pathway].scores.max()
    score = score_plus - score_minus
    return chems_of_interest, score




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


def knapSack_d(negative_score_limit, negative_score, positive_score, chems, number_of_chems):
    negative_score = negative_score * 1000
    positive_score = positive_score * 1000
    K = [[0 for x in range(negative_score_limit + 1)] for x in range(number_of_chems + 1)]
    chem_set = []
    # Build table K[][] in bottom up manner
    for i in range(number_of_chems + 1):
        for w in range(negative_score_limit + 1):
            if i == 0 or w == 0:
                K[i][w] = 0
            elif negative_score[i-1] <= w:
                K[i][w] = max(positive_score[i-1] + K[i-1][w-negative_score[i-1]],  K[i-1][w])                         
            else:
                K[i][w] = K[i-1][w]
                
    res = K[n][W]
    w = negative_score_limit
    
    for i in range(number_of_chems, 0, -1):
        if res <= 0:
            break
        # either the result comes from the
        # top (K[i-1][w]) or from (val[i-1]
        # + K[i-1] [w-wt[i-1]]) as in Knapsack
        # table. If it comes from the latter
        # one/ it means the item is included.
        print('i =',i, 'w =',w)
        if res == K[i - 1][w]:
            continue
        else: 
            # This item is included.
            chem_set.append(chems[i-1])
            res = res - positive_score[i - 1]
            w = w - negative_score[i - 1] 
  
    return chem_set, K[number_of_chems][negative_score_limit]
 
# end of function knapSack
 
pathway = pathways[4]

def get_chem_score(pathway):
    positive_scores, negative_scores = get_positive_negative_scores(df, pathway)
    scores = positive_scores - negative_scores
    zip_iterator = zip(chems, scores)
    chem_score = dict(zip_iterator)
    return chem_score

chem_score = get_chem_score(pathway)

def find_proper_chem_sets_cutoff_method(lamdas):
    chem_sets = []
    for lamda in lamdas:
        chem_set = []
        for chem in chems:
            if chem_score[chem]+lamda > 0:
                chem_set.append(chem)
        chem_sets.append(chem_set)
    
    return chem_sets


def get_chem_set_score(chem_set):
    total_score = 0
    for chem in chem_set:
        total_score += chem_score[chem]
    return total_score


lamdas = [7]

chem_sets = find_proper_chem_sets_cutoff_method(lamdas)
print(len(chem_sets))
print(chem_sets[0])
score = 0
newdf = df.loc[df.pathways==pathway]
for chem in chem_sets[0]:
    score += newdf.loc[newdf.chems==chem].scores.values[0]
print(score)

