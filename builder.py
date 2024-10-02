import networkx as nx
from rdkit import Chem
from objects import Atom

def molecule2graph(mol) -> nx.Graph:
    g = nx.Graph()
    idx2atom = {}
    for chem_atom in mol.GetAtoms():
        atom = Atom.from_chem(chem_atom)
        idx2atom[atom.idx] = atom
        g.add_node(atom)

    for bond in mol.GetBonds():
        start_idx = bond.GetBeginAtomIdx()
        start = idx2atom[start_idx]
        end_idx = bond.GetEndAtomIdx()
        end = idx2atom[end_idx]
        bond_type = str(bond.GetBondType())
        g.add_edge(start, end, bond_type=bond_type)

    return g

def graph2molecule(g: nx.Graph):
    mol = Chem.RWMol()
    node_to_idx = {}
    for atom in g.nodes():
        chem_atom = atom.to_chem()
        idx = mol.AddAtom(chem_atom)
        node_to_idx[atom] = idx
    for start, end, data in g.edges(data=True):
        start_idx, end_idx = node_to_idx[start], node_to_idx[end]
        bond_type = getattr(Chem.rdchem.BondType, data['bond_type'])
        mol.AddBond(int(start_idx), int(end_idx), order=bond_type)

    mol = mol.GetMol()
    Chem.SanitizeMol(mol)
    return mol