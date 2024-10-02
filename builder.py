import networkx as nx
from rdkit import Chem

def molecule2graph(mol) -> nx.Graph:
    g = nx.Graph()
    for atom in mol.GetAtoms():
        g.add_node(atom.GetIdx(),
                   atomic_num=atom.GetAtomicNum(),
                   symbol=atom.GetSymbol(),
                   charge=atom.GetFormalCharge())

    for bond in mol.GetBonds():
        start = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        bond_type = bond.GetBondType()
        g.add_edge(start, end, bond_type=str(bond_type))

    return g

def graph2molecule(g: nx.Graph):
    mol = Chem.RWMol()
    node_to_idx = {}
    for node, data in g.nodes(data=True):
        atom = Chem.Atom(data['atomic_num'])
        atom.SetFormalCharge(data['charge'])
        idx = mol.AddAtom(atom)
        node_to_idx[node] = idx
    for start, end, data in g.edges(data=True):
        start, end = node_to_idx[start], node_to_idx[end]
        bond_type = getattr(Chem.rdchem.BondType, data['bond_type'])
        mol.AddBond(start, end, bond_type=bond_type)

    mol = mol.GetMol()
    Chem.SanitizeMol(mol)
    return mol