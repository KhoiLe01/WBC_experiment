import maxflow
from typing import List, Tuple

def Find_Densest_Subgraph(nodes: List[int], edges: List[Tuple[int, int]], c: float):
    ''' This function performs the binary search of the density of subgraph and finds the densest subgraph.
        Nodes can be arbitrary hashable identifiers (not necessarily 0..n-1).
    '''
    num_nodes = len(nodes) # n
    num_edges = len(edges) # m
    min_degree = 0.0

    max_degree = float(num_edges) 
    # the smallest difference between two possible densities is 1 / (n * (n-1))
    if num_nodes > 1:
        difference = 1.0 / (num_nodes * (num_nodes - 1))
    else:
        difference = 1.0
    
    subgraph = []
    
    # Binary search on the density
    while(max_degree - min_degree >= difference):
        least_density = (max_degree + min_degree) / 2.0
        
        # Solve the max-flow min-cut problem for this density
        flow, source_segment = solve(nodes, edges, least_density)
        
        if flow < 2 * num_edges - 2 * c * least_density: 
            min_degree = least_density
            subgraph = source_segment
        else:
            max_degree = least_density

    return subgraph

def solve(nodes: List[int], edges: List[Tuple[int, int]], g: float) -> Tuple[float, List[int]]:
    ''' Constructs the network as per the specifications given by Goldberg algorithm'''
    num_nodes = len(nodes)
    num_edges = len(edges)
    
    # Create the graph with estimated number of nodes and edges
    graph = maxflow.Graph[float](num_nodes, num_edges + 2*num_nodes)
    nodes_map = graph.add_nodes(num_nodes)
    
    # Map arbitrary node labels to internal 0..n-1 indices
    node_to_idx = {node: i for i, node in enumerate(nodes)}
    
    degrees = {node: 0 for node in nodes}
    
    # Edges between v and u have capacity 1
    for u, v in edges:
        if u not in degrees or v not in degrees:
            continue # Skip edges with nodes not in the node list
            
        degrees[u] += 1
        degrees[v] += 1
        
        # Capacity of each edge (u,v) is 1 in both directions (undirected)
        idx_u = node_to_idx[u]
        idx_v = node_to_idx[v]
        graph.add_edge(nodes_map[idx_u], nodes_map[idx_v], 1.0, 1.0)
    
    for node in nodes:
        idx = node_to_idx[node]
        # Edges from source S to node i with capacity = degree of i
        # Edges from node i to sink T with capacity = 2 * g
        graph.add_tedge(nodes_map[idx], float(degrees[node]), 2.0 * g)
        
    flow = graph.maxflow()
    
    source_segment = []
    for node in nodes:
        idx = node_to_idx[node]
        if graph.get_segment(nodes_map[idx]) == 0: # 0 means Source side
            source_segment.append(node)
            
    return flow, source_segment

