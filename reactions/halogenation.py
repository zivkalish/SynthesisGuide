import networkx as nx
from typing import List
from objects import Atom
import copy
from reactions.utils import count_available_hydrogens, count_carbon_bonds


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