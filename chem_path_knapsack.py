import pandas as pd
import numpy as np

df = pd.read_csv("chem_path_scores_df_fixed.csv")
df.drop('Unnamed: 0', axis=1, inplace=True)
df.scores = abs(df.scores.round(3))
df = df.sort_values(['chems', 'scores'], ascending=False)

pathways = df.pathways.unique()
chems = df.chems.unique()

def get_positive_negative_scores(df, pathway):
    positive_scores = []
    negative_scores = []

    for chem in chems:
        tmp_df = df[df.chems == chem]
        positive_score = tmp_df.loc[tmp_df.pathways == pathway].scores.values[0]
        negative_score = tmp_df.loc[tmp_df.pathways != pathway].scores.sum()
        positive_scores.append(int(positive_score*1000))
        negative_scores.append(int(negative_score*1000))
    return positive_scores, negative_scores



def knapSack_d(negative_score_limit, negative_score, positive_score, chems, number_of_chems):
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
positive_scores, negative_scores = get_positive_negative_scores(df, pathway)

W = 500000
n = len(positive_scores)

chem_set, score = knapSack_d(W, negative_scores, positive_scores, chems, n)
print(len(chem_set))
print(score)
score = 0
newdf = df.loc[df.pathways==pathway]
for chem in chem_set:
    score += newdf.loc[newdf.chems==chem].scores.values[0]
print(score)