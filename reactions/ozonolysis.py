import networkx as nx
from typing import List, Tuple
from objects import Atom
import copy
import networkx_utils as nx_utils

def find_double_bonds(molecule: nx.Graph) -> List[Tuple[Atom, Atom]]:
    return [(u, v) for (u, v, data) in molecule.edges(data=True)
            if data['bond_type'] == "DOUBLE" and u.symbol == "C" and v.symbol == "C"]

def ozonolysis(molecule: nx.Graph, oxidative: bool) -> List[nx.Graph]:
    molecule = copy.deepcopy(molecule)
    double_bonds = find_double_bonds(molecule)
    if not double_bonds:
        return [molecule]
    res = []
    for bond in double_bonds:
        cleaved_molecules = cleave_bond(molecule, bond, oxidative)
        res.extend(cleaved_molecules)
    return res

def cleave_bond(molecule: nx.Graph, bond: Tuple[Atom, Atom], oxidative: bool) -> List[nx.Graph]:
    molecule = copy.deepcopy(molecule)
    atom1, atom2 = bond
    molecule.remove_edge(atom1, atom2)
    for atom in (atom1, atom2):
        add_oxygen(molecule, atom, "DOUBLE")
        if oxidative:
            add_oxygen(molecule, atom, "SINGLE")
    return nx_utils.get_connected_components(molecule)

def add_oxygen(molecule: nx.Graph, carbon: Atom, bond_type:str):
    oxygen = Atom(
        idx=max(atom.idx for atom in molecule.nodes) + 1,
        symbol="O",
        atomic_num=8,
        charge=0
    )
    molecule.add_node(oxygen)
    molecule.add_edge(carbon, oxygen, bond_type=bond_type)
