import copy
import random
import time
import numpy as np
import heapq
from subprocess import check_call

def l_infinity(a, b):
    return max(abs(a_i - b_i) for a_i, b_i in zip(a, b))

def baseline1(vertices, edges, maximum_time=0):
    start_time = time.perf_counter()
    
    if not vertices:
        return [], {}

    # Build adjacency list
    adj_list = {}
    for u, v in edges:
        try:
            adj_list[u].add(v)
        except KeyError:
            adj_list[u] = {v}
        try:
            adj_list[v].add(u)
        except KeyError:
            adj_list[v] = {u}

    vertices_set = set(vertices)
    neighbors_map = {}
    c_initial = set()

    for v in vertices:
        nbrs = adj_list.get(v, set())
        neighbors_map[v] = nbrs
        if nbrs:
            c_initial.add(frozenset(nbrs))

    candidates = set()
    
    processed_hashes = set()
    
    def process_neighborhood(s):
        h = hash(s)
        if h in processed_hashes:
            return False
        processed_hashes.add(h)
        
        s_list = list(s)
        if not s_list:
            return False
            
        try:
            first_u = s_list[0]
            if first_u in adj_list:
                common_nbrs = adj_list[first_u].intersection(vertices_set)
            else:
                common_nbrs = set()
            
            if not common_nbrs:
                return False

            for i in range(1, len(s_list)):
                u = s_list[i]
                if u in adj_list:
                    common_nbrs.intersection_update(adj_list[u])
                else:
                    common_nbrs = set()
                    
                if not common_nbrs:
                    break
        except Exception:
            common_nbrs = set()
            
        if common_nbrs:
            candidates.add(frozenset(common_nbrs | s))
            # Why? Originally we wanted to limit the number of candidates.
            # return len(candidates) > threshold
        return False

    # Shuffle the list to get a different order of processing each time, which can lead to different candidates being generated and potentially better coverage. We can also consider a more sophisticated ordering strategy in the future.
    c_list = list(c_initial)
    random.shuffle(c_list)
    c_list_copy = c_list.copy()
    random.shuffle(c_list_copy)

    n = len(c_list)
    for i in range(n):
        s1 = c_list[i]
        for j in range(i + 1, n):
            s2 = c_list_copy[j]
            if s1 == s2:
                continue
            if not s1.isdisjoint(s2):
                intersection = s1.intersection(s2)
                process_neighborhood(intersection)
                checkpoint_time = time.perf_counter()
                if maximum_time > 0 and (checkpoint_time - start_time > maximum_time):
                    return list(candidates), neighbors_map
                
    return list(candidates), neighbors_map

# State of the art algortithm for MBE: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=10633882
def ambea_baseline(input_edges, len_v, timeout=0):
    edges = [(v, u-len_v) for v, u in input_edges]
    max_v = max(u for u, v in edges) + 1
    d = {i: [] for i in range(max_v)}
    for v, u in edges:
        d[v].append(u)
    print("number of edges:", len(edges))
    print(sum(len(d[v]) for v in d))
    l = [" ".join([str(u) for u in d[v]]) for v in d]
    with open('edges.txt', 'w') as f:
        f.write("\n".join(l))
    
    # example
    # ./MBE_ALL -i "test_new" -s 2 > out.txt 
    with open('bicliques.txt', 'w') as f:
        check_call(['MBE_ALL.exe', '-i', 'edges.txt', '-s', '5', '-T', str(timeout), '-S'], stdout=f, cwd='.')
        
    with open('bicliques.txt', 'r') as f:
        content = f.readlines()
    
    candidates = []
    
    for line in content:
        edges_list = line.split("|")
        if len(edges_list) > 1:
            candidate = ()
            v_set = [int(i) for i in edges_list[0].strip().split()]
            u_set = [int(i) + len_v for i in edges_list[1].strip().split()]
            candidate = frozenset(set(v_set) | set(u_set))
            candidates.append(candidate)
    
    return candidates

def greedy(vertices, edges, len_v, baseline: int, c: float, maximum_time=0):
    fraction = 0.1
    start_time = time.perf_counter()
    if baseline == 1:
        candidates, neighbors_map = baseline1(vertices, edges, fraction*maximum_time)
    elif baseline == 2:
        candidates = ambea_baseline(edges, len_v, (fraction/1.5)*maximum_time if maximum_time < 2000 else (fraction*1.25)*maximum_time)
        neighbors_map = {}
        for u, v in edges:
            if u not in neighbors_map:
                neighbors_map[u] = set()
            neighbors_map[u].add(v)
            if v not in neighbors_map:
                neighbors_map[v] = set()
            neighbors_map[v].add(u)
    print(f"Total candidates generated: {len(candidates)}")
    
    # edges to IDs
    edge_list = list(edges)
    num_edges = len(edge_list)
    remaining_edges_count = num_edges
    
    is_edge_covered = bytearray(num_edges) 
    
    # map (u, v) pairs to their edge ID(s)
    # Storing tuple(sorted((u,v))) -> list of IDs
    pair_to_edge_indices = {}
    for idx, (u, v) in enumerate(edge_list):
        key = tuple(sorted((u, v)))
        if key not in pair_to_edge_indices:
            pair_to_edge_indices[key] = []
        pair_to_edge_indices[key].append(idx)
        
    cand_edges_indices = [] # cand_index -> list of edge_indices
    cand_sizes = []
    
    # edge_index -> list of cand_indices
    edge_to_cand_indices = [[] for _ in range(num_edges)]
    
    for c_idx, cand in enumerate(candidates):
        cand_list = sorted(list(cand))
        current_edges = []
        cand_size = len(cand_list)
        cand_sizes.append(cand_size)
        
        # Identify edges within this candidate
        for i in range(cand_size):
            u = cand_list[i]
            if u not in neighbors_map: continue
            
            nbrs = neighbors_map[u]
            for j in range(i + 1, cand_size):
                v = cand_list[j]
                if v in nbrs:
                    # Edge exists, retrieve ID
                    key = tuple(sorted((u, v)))
                    if key in pair_to_edge_indices:
                        for eid in pair_to_edge_indices[key]:
                            current_edges.append(eid)
                            edge_to_cand_indices[eid].append(c_idx)
                            
        cand_edges_indices.append(current_edges)
    
    # track uncovered edges count for each candidate
    cand_uncovered_count = [len(ce) for ce in cand_edges_indices]
    
    pq = [] 
    # Heap stores (ratio, c_idx)
    
    for c_idx in range(len(candidates)):
        cnt = cand_uncovered_count[c_idx]
        if cnt > 0:
            ratio = (cand_sizes[c_idx] + c) / cnt
            heapq.heappush(pq, (ratio, c_idx))
            
    selected_candidates = []
    
    # phase 1: lazy greedy selection
    while remaining_edges_count > 0 and pq:
        
        # uncomment this to enforce a hard timeout on the greedy phase.
        # current_time = time.perf_counter()
        # if maximum_time > 0 and (current_time - start_time > 1.2*maximum_time):
        #     # return (sum(len(c) for c in selected_candidates) + 2 * remaining_edges_count, False)
        #     break
            
        ratio, c_idx = heapq.heappop(pq)
        
        cnt = cand_uncovered_count[c_idx]
        if cnt == 0:
            continue
            
        # update check
        current_ratio = (cand_sizes[c_idx] + c) / cnt
        if ratio < current_ratio - 1e-9:
            heapq.heappush(pq, (current_ratio, c_idx))
            continue
            
        # select candidate
        selected_candidates.append(candidates[c_idx])
        cand_uncovered_count[c_idx] = 0 # Mark as selected/processed
        
        # update data structures for covered edges
        edges_to_remove = []
        for eid in cand_edges_indices[c_idx]:
            if is_edge_covered[eid] == 0:
                is_edge_covered[eid] = 1
                edges_to_remove.append(eid)
                remaining_edges_count -= 1
        
        if not edges_to_remove:
            continue
            
        # update state for covered edges
        for eid in edges_to_remove:
            u, v = edge_list[eid]
            if u in neighbors_map: neighbors_map[u].discard(v)
            if v in neighbors_map: neighbors_map[v].discard(u)
            
            for impacted_c_idx in edge_to_cand_indices[eid]:
                if cand_uncovered_count[impacted_c_idx] > 0:
                    cand_uncovered_count[impacted_c_idx] -= 1
    
    print(f"Remaining edges after Phase 1: {remaining_edges_count}")
    
    # phase 2: we greedily select vertex in v that covers the most remaining edges until we cover all edges.
    if remaining_edges_count > 0:
        vertex_pq = []
        for v in neighbors_map:
            d = len(neighbors_map[v])
            if d > 0:
                cost = (d + 1 + c) / d
                heapq.heappush(vertex_pq, (cost, v))
                
        while remaining_edges_count > 0 and vertex_pq:
             
            # uncomment this to enforce a hard timeout on the greedy phase. Not recommend here
            # current_time = time.perf_counter()
            # if maximum_time > 0 and (current_time - start_time > maximum_time):
            #     return (sum(len(c) for c in selected_candidates) + 2 * remaining_edges_count, False)

            cost, v = heapq.heappop(vertex_pq)
            
            if v not in neighbors_map:
                continue
                
            current_nbrs = neighbors_map[v]
            d = len(current_nbrs)
            
            if d == 0:
                del neighbors_map[v]
                continue
                
            real_cost = (d + 1 + c) / d
            if cost < real_cost - 1e-9:
                heapq.heappush(vertex_pq, (real_cost, v))
                continue
                
            star_cand = set(current_nbrs)
            star_cand.add(v)
            selected_candidates.append(frozenset(star_cand))
            
            # remove edges
            remaining_edges_count -= d
            
            # update neighbors
            for n in list(current_nbrs):
                if n in neighbors_map:
                    neighbors_map[n].discard(v)
                    # We rely on lazy update for n's degree in PQ
            
            del neighbors_map[v]
    
    return sum(len(i) for i in selected_candidates), True
