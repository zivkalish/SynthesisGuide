import networkx as nx
from objects import Atom

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