import argparse
import matplotlib.pyplot as plt
import numpy as np
import random
import time

# SSL for fetching dataset
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

# Algorithms
from algorithm.r_tree import RTree
from algorithm.baseline import greedy

# Load data
from data.data_loader import load_data_adults, load_data_credits, load_data_gamma, load_data_popsim

# Tree for processing edges
from scipy.spatial import KDTree as SciPyKDTree

EPSILON = 1e-9

def run_one_iteration(delta, sample_size, algo, baseline_mode, c, distance_metric, max_time, v, u):
    """
    Runs a single iteration of the experiment.
    """
    len_v = len(v)
    len_u = len(u)
    vertices = [i for i in range(len_v + len_u)]
    edges = []

    v_np = np.array(v)
    u_np = np.array(u)

    tree = SciPyKDTree(u_np)

    indices = tree.query_ball_point(v_np, r=delta + EPSILON, p=float(distance_metric))

    edges = []
    v_offset = len(v_np)
    for i, neighbors in enumerate(indices):
        for neighbor_idx in neighbors:
            edges.append((i, neighbor_idx + v_offset))
    
    print(
        f"Delta: {delta}, Data Size: {sample_size}, Generated graph with {len(vertices)} vertices, length v {len_v} and {len(edges)} edges."
    )
    result = {
        "delta": delta,
        "sample_size": sample_size,
        "r_approx_res": 0,
        "r_approx_time": 0,
        "baseline_res": 0,
        "baseline_time": 0,
    }

    # Approx r-Tree Approach
    if algo == "ouralgo":
        start_r_approx = time.perf_counter()
        r_tree_approx = RTree(v, u, delta, c, distance_metric)
        end_build = time.perf_counter()
        print(f"R-Tree built in {(end_build - start_r_approx):.4f}s")
        res, _ = r_tree_approx.solve(edges)
        end_r_approx = time.perf_counter()
        result["r_approx_res"] = res
        result["r_approx_time"] = end_r_approx - start_r_approx
        print(
            f"Delta: {delta}, Data Size: {sample_size}, R-Tree Result: {res}, Time: {(end_r_approx - start_r_approx):.4f}s"
        )

    # --- Baseline Approach ---
    if algo == "baseline":
        start_baseline = time.perf_counter()
        if not edges:
            baseline_res = 0
        else:
            baseline_res, completed = greedy(
                vertices, edges, len_v, baseline_mode, c, max_time
            )

        end_baseline = time.perf_counter()
        result["baseline_res"] = baseline_res
        result["baseline_time"] = end_baseline - start_baseline
        print(
            f"Delta: {delta}, Data Size: {sample_size}, Baseline Result: {result['baseline_res']}, Time: {result['baseline_time']:.4f}s"
        )

    return result

def main(dataset, c, delta_list, algo, baseline_mode, seed, v_min, v_max, u_min, u_max, distance_metric, max_time):
    # Parameters

    match dataset:
        case "adults":
            data_list, v, u = load_data_adults()
        case "credits":
            data_list, v, u = load_data_credits()
        case "gamma":
            data_list, v, u = load_data_gamma()
        case "popsim":
            data_list, v, u = load_data_popsim(seed, v_min, v_max, u_min, u_max)
    random.seed(seed)
    sample_size = len(data_list)
    # Prepare tasks
    tasks = []
    for delta in delta_list:
        tasks.append((delta, sample_size, algo, baseline_mode, c, distance_metric, max_time, v, u))

    print(f"Starting execution of {len(tasks)} tasks sequentially...")

    results_aggregated = {}

    # structure to hold results for each delta and sample size
    # can be used for debugging and later for plotting
    for delta in delta_list:
        results_aggregated[delta] = {
            "sample_sizes": [sample_size],
            "r_times_approx": {sample_size: []},
            "baseline_times": {sample_size: []},
            "r_results_approx": {sample_size: []},
            "baseline_results": {sample_size: []},
        }

    start_total = time.perf_counter()

    # Run sequentially
    # Technically can run in parallel, but recorded memory usage will be different.
    results_raw = []
    for task in tasks:
        results_raw.append(run_one_iteration(*task))
        time.sleep(10)  # Small delay to track logs more easily
    end_total = time.perf_counter()
    print(f"Total execution time: {end_total - start_total:.2f}s")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to run the experiment.")
    parser.add_argument(
        "--dataset",
        type=str,
        choices=["adults", "credits", "gamma", "popsim"],
        default="adults",
        help="Dataset to use for the experiment.",
    )
    parser.add_argument(
        "--distance_metric",
        type=str,
        default="inf",
        help="Distance metric to use for graph construction (default: inf). Use 1 for L1, 2 for L2, etc., and inf for L-infinity.",
    )
    parser.add_argument(
        "--c",
        type=int,
        default=0,
        help="0 for WBC, above 0 for GWBC",
    )
    parser.add_argument(
        "--delta_list",
        type=float,
        nargs="+",
        default=[0.13, 0.14, 0.15, 0.16],
        help="List of delta values to test."
    )
    parser.add_argument(
        "--algo",
        type=str,
        choices=["baseline", "ouralgo"],
        default="baseline",
        help="Algorithm to run: baseline or ouralgo.",
    )
    parser.add_argument(
        "--baseline_mode",
        type=int,
        choices=[1, 2],
        default=1,
        help="Baseline mode to use (1-2). 1 is for baseline 1, 2 is for AMBEA baseline. No need to specify if algo is not baseline.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=2026,
        help="Random seed for reproducibility.",
    )
    parser.add_argument(
        "--v_min",
        type=int,
        default=4000,
        help="Minimum number of vertices in set V to use (for sampling). Only applicable for the popsim dataset.",
    )
    parser.add_argument(
        "--v_max",
        type=int,
        default=6000,
        help="Maximum number of vertices in set V to use (for sampling). Only applicable for the popsim dataset.",
    )
    parser.add_argument(
        "--u_min",
        type=int,
        default=4000,
        help="Minimum number of vertices in set U to use (for sampling). Only applicable for the popsim dataset.",
    )
    parser.add_argument(
        "--u_max",
        type=int,
        default=6000,
        help="Maximum number of vertices in set U to use (for sampling). Only applicable for the popsim dataset.",
    )
    parser.add_argument(
        "--max_time",
        type=int,
        default=3600,
        help="Maximum time (in seconds) to allow for the baseline algorithm to run before terminating it.",
    )
    args = parser.parse_args()
    main(args.dataset, args.c, args.delta_list, args.algo, args.baseline_mode, args.seed, args.v_min, args.v_max, args.u_min, args.u_max, args.distance_metric, args.max_time)