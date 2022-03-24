import pandas as pd
from d3graph import d3graph
import networkx as nx
import hnet


df = hnet.import_example('cancer')

print(df.head())

from hnet import hnet
hn = hnet(black_list=['tsneX','tsneY','PC1','PC2'])
results = hn.association_learning(df)

print(hn.results['rules'])

G = hn.d3graph()