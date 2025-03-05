import networkx as nx
from typing import List
from objects import Atom


def safe_copy_graph(graph):
    """Safely copy a molecular graph without issues with RDKit objects."""
    new_graph = nx.Graph()

    # Copy nodes
    for node in graph.nodes():
        new_graph.add_node(Atom(
            idx=node.idx,
            atomic_num=node.atomic_num,
            symbol=node.symbol,
            charge=node.charge
        ))

    # Build a mapping from old nodes to new nodes
    old_to_new = {old: new for old, new in zip(graph.nodes(), new_graph.nodes())}

    # Copy edges
    for u, v, data in graph.edges(data=True):
        new_graph.add_edge(old_to_new[u], old_to_new[v], **data)

    return new_graph


def get_connected_components(graph: nx.Graph) -> List[nx.Graph]:
    connected_components = list(nx.connected_components(graph))
    subgraphs = [graph.subgraph(component).copy() for component in connected_components]
    return subgraphs