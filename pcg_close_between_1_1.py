"""
Protein Contact Graphs (PCG) represent proteins as networks and enable topological analysis of complex 3D structures, regardless of 
secondary structure or fold type. They provide insight into structure-function relationships and can be used in conjunction with other 
methods. PCGs are unweighted graphs where nodes represent amino acids and edges are drawn between nodes if the Cα-Cα distance is within 
a certain threshold, allowing for comparison and visualization of structural features and identification of long-range interactions in proteins.


PCG helps to:

-identify the putative allosteric paths and regions in protein structures,
-Understand functional effects of mutations,
-generate information about ligand binding,
-Identify functional domains (communities or clusters) in proteins.

Centrality measures, such as betweenness and closeness, are useful in identifying crucial residues for protein folding and function. 
Betweenness centrality is based on communication flow, and nodes with high betweenness centrality control information flow and are 
often located in regions of experimentally validated hot spots in protein-protein interactions. However, choosing a centrality measure 
is not straightforward as there is lack of consensus on which measure performs best in identifying relevant binding residues. Code can 
be used to apply these measures to protein contact graphs to identify important amino acid residues in cavities, which can be useful 
for protein folding.(del Sol and O’Meara, 2005)
"""

import numpy as np
import subprocess
import networkx as nx
import utils
import pandas as pd
import pymol_utils

from Bio.PDB import *
import Bio
from Bio.PDB import PDBList

## Download pdb files from PDB
pdbl = PDBList()
PDBlist2 = ['6vie', '6n9o']
for i in PDBlist2:
    pdbl.retrieve_pdb_file(i,pdir='.', file_format ='pdb')

import glob, os
folder = os.getcwd()
for filename in glob.iglob(os.path.join(folder, '*.ent')):
    os.rename(filename, filename[:-4] + '.pdb')

## define pdb in local folder
protein = "pdb6vie"
protein_path = "{}.pdb".format(protein)
atoms = utils.readPDBFile(protein_path) #read 
coordinates = utils.getResidueCoordinates(atoms)
residue_names = np.array(list (dict_residue_name.items()))

## Compute adjacnet matrix with the coordinates between 4-8A

output_path = "" ## Local directory
A = utils.adjacent_matrix(output_path, coordinates, protein,4, 8)


## Generation of different centrality measures

G = nx.from_numpy_array(A)
print(G)

res_name = np.array(residue_names[:, 1], dtype = str)  
centrality_measures = utils.betweenness(G, res_names=res_name,n=30) 
#centrality_measures = utils.closeness(G, res_names=res_name,n=30)
#centrality_measures = utils.pagerank_ct(G, res_names=res_name,n=30,alpha=0.8)
#centrality_measures = utils.eigenvector_ct(G, res_names=res_name,n=30)
df = pd.DataFrame.from_dict(centrality_measures,orient='index')
df.sort_values(by=df.columns[0],ascending=False)
df.rename(columns={ df.columns[0]: "betweenness" }, inplace = True)
df.to_csv("betweenness.csv")

## Generation of Pymol session file with centrality scores

pymol_utils.pymol_centralities(output_path, centrality_measures, protein_path, "closeness")


"""
SUMMARY
PCG Network analysis can help in :

-Identifying the major residues in the binding of 2 or more proteins is a crucial task to understand their function, plan mutagenesis 
experiments, or target drug design.
-Network analysis can reveal major binding residues of PPI and Ligand Binding.
-To identify residues highly correlate with functional relevance and constitute a good set of targets for mutagenesis experiments.
-To detect allosteric sites and functional regions activating upon binding as described in (Di Paola et al J Prot Res 2020) that demonstrated 
that network clusters (group of residues) correspond to functional regions in the protein. This method is based on network spectral clustering.
-Tools like laplacian Eigenmaps, Graph Embedding Algorithms can be studied in protein behavior .
"""