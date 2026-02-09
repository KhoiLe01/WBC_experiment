
from typing import List, Tuple
from densest_subgraph.c_densest_subgraph import Find_Densest_Subgraph

def two_approx_densest_subgraph(nodes: List[int], edges: List[Tuple[int, int]]) -> List[int]:
    ''' Returns a 2-approximation of the densest subgraph using the greedy peeling algorithm. Optimized O(N+E) implementation. '''
    if not nodes:
        return []

    n = len(nodes)
    # Map nodes to 0..n-1 for efficient array indexing
    node_to_idx = {node: i for i, node in enumerate(nodes)}

    # Build Adjacency List
    adj_sets = [set() for _ in range(n)]
    degree = [0] * n

    for u, v in edges:
        if u not in node_to_idx or v not in node_to_idx:
            continue
        idx_u = node_to_idx[u]
        idx_v = node_to_idx[v]

        if idx_u != idx_v:
            if idx_v not in adj_sets[idx_u]:
                adj_sets[idx_u].add(idx_v)
                adj_sets[idx_v].add(idx_u)
                degree[idx_u] += 1
                degree[idx_v] += 1

    # Convert to list for faster iteration
    adj = [list(s) for s in adj_sets]

    # Initialize Buckets
    max_deg = 0
    if degree:
        max_deg = max(degree)

    buckets = [set() for _ in range(max_deg + 1)]
    for i in range(n):
        buckets[degree[i]].add(i)

    min_deg = 0
    while min_deg <= max_deg and not buckets[min_deg]:
        min_deg += 1

    # Peeling Process
    curr_n = n
    curr_e = sum(degree) // 2

    best_density = curr_e / curr_n if curr_n > 0 else 0
    best_iter = 0

    removed = [False] * n
    peeling_order = []

    for i in range(n):
        # Check density of current state (before removal)
        density = curr_e / curr_n if curr_n > 0 else 0
        if density > best_density:
            best_density = density
            best_iter = i

        # Find next node to remove
        while min_deg <= max_deg and not buckets[min_deg]:
            min_deg += 1

        if min_deg > max_deg:
            break

        v = buckets[min_deg].pop()
        removed[v] = True
        peeling_order.append(v)

        # Update current stats
        deg_v = degree[v]  # This is degree in remaining graph
        curr_n -= 1
        curr_e -= deg_v

        # Update neighbors
        for u in adj[v]:
            if not removed[u]:
                deg_u = degree[u]
                buckets[deg_u].remove(u)

                deg_u -= 1
                degree[u] = deg_u
                buckets[deg_u].add(u)

                if deg_u < min_deg:
                    min_deg = deg_u

    # Reconstruct result based on best_iter
    removed_set = set(peeling_order[:best_iter])
    result = [nodes[i] for i in range(n) if i not in removed_set]
    return result