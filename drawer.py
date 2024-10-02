import networkx as nx
import matplotlib.pyplot as plt
from rdkit.Chem import Draw

import consts
from rdkit.Chem.rdchem import Mol
from builder import graph2molecule

def draw_graph(g: nx.Graph):
    nodes_colors = [consts.COLOR_MAP[g.nodes[atom].symbol] for atom in g.nodes]
    labels = {atom: g.nodes[atom].symbol for atom in g.nodes}
    edge_labels = {(u, v): g.edges[u, v]['bond_type'] for u, v in g.edges}
    pos = nx.spring_layout(g)
    nx.draw(g, pos, with_labels=True, labels=labels, node_color=nodes_colors)
    nx.draw_networkx_edge_labels(g, pos, edge_labels=edge_labels)
    plt.show()

def draw_graph_as_mol(g: nx.Graph):
    mol = graph2molecule(g)
    draw_mol(mol)

def draw_mol(mol: Mol):
    img = Draw.MolToImage(mol)
    plt.imshow(img)
    plt.axis('off')
    plt.show()