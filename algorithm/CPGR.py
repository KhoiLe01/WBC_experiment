from copy import deepcopy
import os
from subprocess import run
import time
import math

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "datasets")
CPGR_BINARY_PATH = os.path.join(SCRIPT_DIR, "..", "cpgr.exe")
DEFAULT_DENSITY = 100
DEFAULT_EXPERIMENT_NO = 1
DELTA_STEP = 0.05

def _graph_file_path(left_size, right_size, density=DEFAULT_DENSITY, experiment_no=DEFAULT_EXPERIMENT_NO):
    return os.path.join(
        DATA_DIR,
        f"bipartite_graph_{max(left_size, right_size)}_{density}_{experiment_no}.mtx",
    )


def _write_bipartite_graph(left_size, right_size, edges, density=DEFAULT_DENSITY, experiment_no=DEFAULT_EXPERIMENT_NO):
    os.makedirs(DATA_DIR, exist_ok=True)
    graph_path = _graph_file_path(left_size, right_size, density, experiment_no)
    sorted_edges = sorted(edges)

    with open(graph_path, "w", buffering=1024 * 1024) as f:
        f.write("%%MatrixMarket matrix coordinate pattern general\n")
        f.write(f"{left_size} {right_size} {len(edges)}\n")

        if sorted_edges:
            chunk = []
            append_line = chunk.append
            for row, col in sorted_edges:
                append_line(f"{row} {col}\n")
                if len(chunk) >= 8192:
                    f.write("".join(chunk))
                    chunk.clear()

            if chunk:
                f.write("".join(chunk))

    return graph_path

def _parse_tricliques(temp_file_path, original_node_count):
    cliques = {}

    if not os.path.exists(temp_file_path):
        return []

    with open(temp_file_path, "r") as f:
        for raw_line in f:
            parts = raw_line.split()
            if len(parts) < 2:
                continue

            left_value = int(parts[0])
            right_value = int(parts[1])

            if left_value == -1 or right_value == -1:
                break

            if left_value > original_node_count:
                clique = left_value
                clique_data = cliques.setdefault(clique, {"left": set(), "right": set()})
                clique_data["right"].add(right_value)
            elif right_value > original_node_count:
                clique = right_value
                clique_data = cliques.setdefault(clique, {"left": set(), "right": set()})
                clique_data["left"].add(left_value)

    tricliques = []
    for clique_id in sorted(cliques):
        clique_data = cliques[clique_id]
        if clique_data["left"] and clique_data["right"]:
            tricliques.append((clique_data["left"], clique_data["right"]))

    return tricliques

def _covered_edges_from_tricliques(tricliques, uncovered_edges=None):
    covered_edges = set()
    filtered_tricliques = []
    seen_edges = set()

    for left_vertices, right_vertices in tricliques:
        triclique_edges = {
            (left_vertex, right_vertex)
            for left_vertex in left_vertices
            for right_vertex in right_vertices
        }

        if uncovered_edges is not None:
            triclique_edges &= uncovered_edges

        triclique_edges -= seen_edges
        if not triclique_edges:
            continue

        seen_edges.update(triclique_edges)
        covered_edges.update(triclique_edges)
        filtered_tricliques.append((left_vertices, right_vertices))

    return filtered_tricliques, covered_edges

def _run_cpgr(left_size, right_size, delta, density=DEFAULT_DENSITY, experiment_no=DEFAULT_EXPERIMENT_NO):
    temp_clique_path = os.path.join(DATA_DIR, "tempCliqueEdges.mtx")
    if os.path.exists(temp_clique_path):
        os.remove(temp_clique_path)

    command = [
        CPGR_BINARY_PATH,
        str(max(left_size, right_size)),
        str(density),
        str(experiment_no),
        f"{delta:.12f}",
    ]
    result = run(command, cwd=SCRIPT_DIR, capture_output=True, text=True)

    tricliques = _parse_tricliques(temp_clique_path, max(left_size, right_size))
    return tricliques, result

def run_cpgr(left_size, right_size, edges, eps=DELTA_STEP):
    if eps <= 0:
        raise ValueError("eps must be positive")

    start_time = time.perf_counter()
    
    n = max(left_size, right_size)
    chosen_tricliques = []
    current_edges = deepcopy(edges)
    
    
    while current_edges:
        delta_lb = math.log((2 * n**2)/len(current_edges), 2) / math.log(n, 2)
        print(f"Current edges: {len(current_edges)}, Delta lower bound: {delta_lb:.4f}")
        if delta_lb > 1:
            print(f"Delta lower bound {delta_lb:.4f} exceeds 1. No tricliques can be found. Ending extraction.")
            break
        _write_bipartite_graph(left_size, right_size, current_edges)

        best_delta = None
        best_covered_edges = set()
        best_ratio = 0
        best_tricliques = []

        candidate_delta = max(delta_lb, 1.0)
        while candidate_delta > delta_lb:
            tricliques, _ = _run_cpgr(left_size, right_size, candidate_delta)
            tricliques, covered_edges = _covered_edges_from_tricliques(tricliques, current_edges)
            nodes_used = sum(len(left) + len(right) for left, right in tricliques)
            
            if covered_edges and len(covered_edges)/nodes_used > best_ratio:
                best_delta = candidate_delta
                best_covered_edges = covered_edges
                best_tricliques = tricliques
                best_ratio = len(covered_edges)/nodes_used

            candidate_delta = round(candidate_delta - eps, 12)

        if not best_covered_edges:
            break

        current_edges.difference_update(best_covered_edges)
        chosen_tricliques.append(
            {
                "delta": best_delta,
                "tricliques": best_tricliques,
                "triclique_node_counts": [
                    len(left_vertices) + len(right_vertices)
                    for left_vertices, right_vertices in best_tricliques
                ],
                "covered_edges": len(best_covered_edges),
                "remaining_edges": len(current_edges),
            }
        )

    _write_bipartite_graph(left_size, right_size, current_edges)
    print(
        f"Tripartite extraction complete. Selected {len(chosen_tricliques)} triclique iteration(s); "
        f"remaining edges: {len(current_edges)}"
    )
    end_time = time.perf_counter()
    print(f"Total execution time: {end_time - start_time:.2f} seconds")

    all_triclique_node_counts = []
    for iteration in chosen_tricliques:
        all_triclique_node_counts.extend(iteration["triclique_node_counts"])

    return {
        "triclique_node_counts": all_triclique_node_counts,
        "remaining_edges": len(current_edges),
        "remaining_edge_list": sorted(current_edges),
        "iterations": len(chosen_tricliques),
        "details": chosen_tricliques,
    }