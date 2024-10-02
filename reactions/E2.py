import networkx as nx
from typing import Tuple, Optional
from objects import Atom
import consts
from reactions.utils import count_available_hydrogens, count_carbon_bonds
import copy

def find_leaving_group(molecule: nx.Graph) -> Tuple[Atom, Atom] | Tuple[None, None]:
    for atom in molecule.nodes():
        if not is_leaving_group(atom):
            continue
        neighbors = list(molecule.neighbors(atom))
        if len(neighbors) == 1 and neighbors[0].symbol == 'C':
            return atom, neighbors[0]
    return None, None

def is_leaving_group(atom: Atom) -> bool:
    # TODO: add logic to find more complex leaving groups such as OTs
    if atom.symbol in consts.LEAVING_GROUPS:
        return True
    return False

def find_beta_carbon(molecule: nx.Graph, carbon: Atom, zaitsev: bool) -> Optional[Atom]:
    candidate_carbons = []
    for neighbor in molecule.neighbors(carbon):
        if neighbor.symbol != 'C':
            continue
        if count_available_hydrogens(molecule, neighbor) == 0:
            continue
        carbon_bonds = count_carbon_bonds(molecule, neighbor)
        candidate_carbons.append({"atom":neighbor, "carbon_bonds":carbon_bonds})
    if not candidate_carbons:
        return None
    candidate_carbons.sort(key=lambda c: c["carbon_bonds"], reverse=zaitsev)
    return candidate_carbons[0]["atom"]



def e2_reaction(molecule: nx.Graph, zaitsev: bool = True) -> nx.Graph:
    leaving_group, alpha_carbon = find_leaving_group(molecule)
    if not leaving_group:
        return molecule

    beta_carbon = find_beta_carbon(molecule, alpha_carbon, zaitsev)
    if not beta_carbon:
        return molecule

    molecule = copy.deepcopy(molecule)
    molecule.remove_node(leaving_group)
    elevate_bond(molecule, alpha_carbon, beta_carbon)
    return molecule

def elevate_bond(molecule: nx.Graph, alpha_carbon: Atom, beta_carbon: Atom):
    current_bond = molecule.edges[alpha_carbon, beta_carbon]['bond_type']
    if current_bond == "SINGLE":
        new_bond = "DOUBLE"
    elif current_bond == "DOUBLE":
        new_bond = "TRIPLE"
    elif current_bond == "TRIPLE":
        raise ValueError("fully bonded connection")
    else:
        raise ValueError(f"Unknown bond type: {current_bond}")
    molecule.remove_edge(alpha_carbon, beta_carbon)
    molecule.add_edge(alpha_carbon, beta_carbon, bond_type=new_bond)


