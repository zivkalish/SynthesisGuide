import networkx as nx
from typing import List



def get_connected_components(graph: nx.Graph) -> List[nx.Graph]:
    connected_components = list(nx.connected_components(graph))
    subgraphs = [graph.subgraph(component).copy() for component in connected_components]
    return subgraphs