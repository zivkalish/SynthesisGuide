import networkx as nx
from typing import List
from objects import Atom
import copy

def count_available_hydrogens(molecule: nx.Graph, atom: Atom) -> int:
    if atom.symbol != "C":
        raise ValueError("intended to work only on carbon atoms")
    total = 4
    for _, neighbor, data in molecule.edges(atom, data=True):
        if neighbor.symbol == 'H':
            continue
        if data['bond_type'] == 'SINGLE':
            total -= 1
        elif data['bond_type'] == 'DOUBLE':
            total -= 2
        elif data['bond_type'] == 'TRIPLE':
            total -= 3
    if total < 0:
        raise ValueError("can't have negative number of hydrogens")
    return total

def count_carbon_bonds(molecule: nx.Graph, atom: Atom) -> int:
    if atom.symbol != 'C':
        raise ValueError("intended to work only on carbon atoms")
    total = 0
    for _, neighbor, data in molecule.edges(atom, data=True):
        if neighbor.symbol != 'C':
            continue
        if data['bond_type'] == 'SINGLE':
            total += 1
        elif data['bond_type'] == 'DOUBLE':
            total += 2
        elif data['bond_type'] == 'TRIPLE':
            total += 3
    if total > 4:
        raise ValueError("carbon can't have more than 4 effective bonds")
    return total

def get_candidate_carbons(g: nx.Graph) -> List:
    """
    Find all carbon atoms with hydrogen
    """
    candidates_carbons = []
    for atom in g.nodes():
        if not atom.symbol == "C":
            continue
        available_hydrogens = count_available_hydrogens(g, atom)
        if available_hydrogens == 0:
            continue
        carbon_bonds = count_carbon_bonds(g, atom)
        if carbon_bonds >= 4:
            continue
        candidates_carbons.append({
            'atom': atom,
            'carbon_bonds': carbon_bonds,
        })

    return candidates_carbons

def halogenation(molecule: nx.Graph, halogen: str) -> nx.Graph:
    halogen_idx = max([atom.idx for atom in molecule.nodes()]) + 1
    if halogen.lower() == "cl":
        halogen = Atom(
            idx=halogen_idx,
            symbol="Cl",
            atomic_num=17,
            charge=0,
        )
    elif halogen.lower() == "br":
        halogen = Atom(
            idx=halogen_idx,
            symbol="Br",
            atomic_num=35,
            charge=0,
        )
    else:
        raise ValueError("halogen must be either Cl or Br")
    molecule = copy.deepcopy(molecule)
    candidate_carbons = get_candidate_carbons(molecule)
    if not candidate_carbons:
        return molecule
    if halogen.symbol == "Cl":
        candidate_carbons.sort(key=lambda carbon: carbon['carbon_bonds'], reverse=False)
    elif halogen.symbol == "Br":
        candidate_carbons.sort(key=lambda carbon: carbon['carbon_bonds'], reverse=True)
    # TODO: add logic what to do in a case where more than one atom have the same available hydrogens
    target_carbon = candidate_carbons[0]['atom']
    molecule.add_node(halogen)
    molecule.add_edge(target_carbon, halogen, bond_type="SINGLE")
    return molecule