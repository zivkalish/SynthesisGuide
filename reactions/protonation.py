import networkx as nx
from typing import Optional
from objects import Atom
import copy
from reactions.utils import has_lone_pair
import consts

def find_protonation_site(molecule: nx.Graph) -> Optional[Atom]:
    candidates = []
    for atom in molecule.nodes():
        if atom.symbol not in consts.SORTED_PROTANBALE_ATOMS:
            continue
        if not has_lone_pair(molecule, atom):
            continue
        candidates.append(atom)
    candidates.sort(key= lambda a: consts.SORTED_PROTANBALE_ATOMS.index(a.symbol))
    return candidates[0] if candidates else None

# def protonation(molectule)


