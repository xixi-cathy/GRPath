print("************************************* Generating figs *************************************")
import pandas as pd
import os
import matplotlib.pyplot as plt
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
import warnings
warnings.filterwarnings("ignore")
init_path = os.getcwd()+ '/..'
os.chdir(init_path)

chf_all = pd.read_csv("regulation_path.csv")

## Plot alternate-allele-pathogenic results
print("Generating figure for alternate-allele-pathogenic regions...")
chf_sub = chf_all.loc[chf_all["type"]=="alternate-allele-pathogenic",:]
rels = []
for i in range(chf_sub.shape[0]):
    rels.append((chf_sub.iloc[i, 3], chf_sub.iloc[i, 4]))
    rels.append((chf_sub.iloc[i, 2], chf_sub.iloc[i, 3]))
    rels.append((chf_sub.iloc[i, 1], chf_sub.iloc[i, 2]))
    if not chf_sub.iloc[i, 1]==chf_sub.iloc[i, 0]:
        rels.append((chf_sub.iloc[i, 0], chf_sub.iloc[i, 1]))

        
rels = list(set(rels))
G=nx.DiGraph()
# G.add_nodes_from(nodes)
G.add_edges_from(rels)
pos=graphviz_layout(G, prog='dot')
f = plt.figure()
#nx.draw(G, pos=pos, ax=f.add_subplot(111), with_labels=True, node_color='#A0E7E5', node_size=400, font_size=4, font_weight='bold', width=5)
nx.draw(G, pos=pos, ax=f.add_subplot(111), with_labels=True, node_color='#A0E7E5', font_size=5, font_weight='bold')
f.set_size_inches(50, 10)
f.savefig("regulation_paths-alternate_pathogenic.pdf")
print("Done.")

## Plot reference-allele-pathogenic results
print("Generating figure for reference-allele-pathogenic regions...")
chf_sub = chf_all.loc[chf_all["type"]=="reference-allele-pathogenic",:]
rels = []
for i in range(chf_sub.shape[0]):
    rels.append((chf_sub.iloc[i, 3], chf_sub.iloc[i, 4]))
    rels.append((chf_sub.iloc[i, 2], chf_sub.iloc[i, 3]))
    rels.append((chf_sub.iloc[i, 1], chf_sub.iloc[i, 2]))
    if not chf_sub.iloc[i, 1]==chf_sub.iloc[i, 0]:
        rels.append((chf_sub.iloc[i, 0], chf_sub.iloc[i, 1]))

        
rels = list(set(rels))
G=nx.DiGraph()
# G.add_nodes_from(nodes)
G.add_edges_from(rels)
pos=graphviz_layout(G, prog='dot')
f = plt.figure()
#nx.draw(G, pos=pos, ax=f.add_subplot(111), with_labels=True, node_color='#A0E7E5', node_size=400, font_size=4, font_weight='bold', width=5)
nx.draw(G, pos=pos, ax=f.add_subplot(111), with_labels=True, node_color='#FFAEBC', font_size=5, font_weight='bold')
f.set_size_inches(20, 10)
f.savefig("regulation_paths-reference_pathogenic.pdf")
print("Done.")

