from typing import List, Tuple, Optional
import numpy as np
import heapq
from densest_subgraph.two_approx_densest_subgraph import two_approx_densest_subgraph
from densest_subgraph.c_densest_subgraph import Find_Densest_Subgraph
from decimal import Decimal

EPSILON = 1e-9

class RTreeNode:
    """Node in an R tree."""

    def __init__(
        self,
        points_v: List[Tuple],
        points_u: List[Tuple],
        delta: float,
    ):
        """
        Initialize a node.
        """
        self.points_v = points_v  # points in V within this node region
        self.points_u = (
            points_u  # points in U in this region
        )
        self.delta = delta
        self.left = None
        self.right = None

        # uncovered edges
        self.previous_uncovered_edges = [(v, u) for v in points_v for u in points_u]
        
        # cache result
        self.previous_density = 0
        self.previous_covered_edge_list = []
        self.initial_call = True

    def is_leaf(self) -> bool:
        """Check if this is a leaf node (contains a point)."""
        return self.points_v is not None and len(self.points_v) == 1


class RTree:
    """
    R Tree
    """

    def __init__(
        self, points_v: List[Tuple], points_u: List[Tuple], delta: float, c: float = 0, distance: str = "inf"
    ):
        self.points_v = points_v
        self.points_u = points_u
        self.delta = delta
        self.distance = distance
        self.c = c
        self.n = len(points_v)
        self.dim = len(points_v[0]) if self.n > 0 else 0
        self.all_nodes = []

        # build a balanced R-tree (using MBRs and widest axis split)
        if len(points_v) > 0:
            pv_arr = np.array(points_v, dtype=np.float64)
        else:
            pv_arr = np.empty((0, 0))

        pu_arr = (
            np.array(points_u, dtype=np.float64)
            if len(points_u) > 0
            else np.empty((0, self.dim))
        )

        self.root = self._build_tree(pv_arr, pu_arr, [], delta)
        print(len(self.all_nodes), "nodes in R-tree built.")

    def _build_tree(
        self,
        points_v_arr: np.ndarray,
        points_u_arr: np.ndarray,
        qualified_u_list: List[Tuple],
        delta: float,
    ):
        m = len(points_v_arr)
        if m == 0:
            return None

        # compute mbrs
        mbr_min = points_v_arr.min(axis=0)
        mbr_max = points_v_arr.max(axis=0)

        # process candidate points_u against MBR of points_v
        new_qualified_list = []
        still_candidates_arr = points_u_arr  # Default if no candidates

        if len(points_u_arr) > 0:
            if self.distance == "inf":
                # condition: u covers MBR if (u <= min + delta) AND (u >= max - delta)
                upper_bound = mbr_min + delta
                lower_bound = mbr_max - delta

                # creating boolean mask for all candidates
                mask = np.all(
                    (points_u_arr >= lower_bound - EPSILON) & (points_u_arr <= upper_bound + EPSILON), axis=1
                )
            else:
                # l_p distance
                p = float(self.distance)
                delta_p = delta**p + EPSILON
                threshold = 1e7  # set a limit to avoid memory issues

                if len(points_u_arr) * len(points_v_arr) * points_v_arr.shape[1] < threshold:
                    diff = np.abs(points_u_arr[:, np.newaxis, :] - points_v_arr[np.newaxis, :, :]) ** p
                    max_dist_p = np.max(np.sum(diff, axis=2), axis=1)
                    mask = max_dist_p <= delta_p
                else:
                    # memory efficient loop over candidates
                    mask = np.zeros(len(points_u_arr), dtype=bool)
                    for i in range(len(points_u_arr)):
                        d_p = np.max(np.sum(np.abs(points_v_arr - points_u_arr[i]) ** p, axis=1))
                        if d_p <= delta_p:
                            mask[i] = True

            newly_qualified_arr = points_u_arr[mask]

            if len(newly_qualified_arr) > 0:
                # convert only the newly qualified points to list of tuples
                new_qualified_list = [tuple(x) for x in newly_qualified_arr]

            still_candidates_arr = points_u_arr[~mask].copy()
        else:
            # maintain correct shape even if empty
            if self.dim > 0:
                still_candidates_arr = np.empty((0, self.dim))
            else:
                still_candidates_arr = np.array([])

        # combine inherited qualified points with newly qualified ones
        current_qualified_list = qualified_u_list + new_qualified_list

        # convert points_v back to list of tuples for Node storage
        points_v_list = [tuple(x) for x in points_v_arr]

        node = RTreeNode(
            points_v=points_v_list,
            points_u=current_qualified_list,
            delta=delta,
        )
        self.all_nodes.append(node)
        if m == 1:
            # leaf
            return node

        # split along widest dimension
        widths = mbr_max - mbr_min
        rounded_widths = np.round(widths, 10) 
        split_axis = np.argmax(rounded_widths)

        # IMPORTANT: python sort is stable, np use quicksort by default. Changing this will provide slightly different results.
        sort_idxs = np.argsort(points_v_arr[:, split_axis], kind='stable')
        sorted_v = points_v_arr[sort_idxs]

        mid = m // 2

        # recurse
        node.left = self._build_tree(
            sorted_v[:mid], still_candidates_arr, current_qualified_list, delta
        )
        node.right = self._build_tree(
            sorted_v[mid:], still_candidates_arr, current_qualified_list, delta
        )

        return node

    def solve(self, edges: List[Tuple]) -> Tuple[int, List[Tuple[set, set]]]:
        """
        Search for all points within the given orthogonal ranges.
        """
        if self.root is None:
            return 0, []

        result = 0
        result_set = []

        # build set of uncovered edges for quick lookup
        all_points = self.points_v + self.points_u
        points_index = {point: idx for idx, point in enumerate(all_points)}
        
        B_int = set(edges)
        
        # priority queue for greedy selection (min-heap based on weight = 1/density)
        pq = []
        
        # initial population of PQ
        for node in self.all_nodes:
            weight, covered_edges = self._calculate_node_weight(
                node, points_index, B_int, self.n
            )
            if weight != float("inf"):
                heapq.heappush(pq, (weight, id(node), node))

        while len(B_int) > 0:
            if not pq:
                # fail safe break
                result += 2 * len(B_int)
                break

            # get best candidate
            weight, _, node = heapq.heappop(pq)
            
            current_weight, current_covered_edges = self._calculate_node_weight(
                node, points_index, B_int, self.n
            )

            # if node quality degraded (higher weight), push back and try next best
            if current_weight > weight:
                if current_weight != float("inf"):
                    heapq.heappush(pq, (current_weight, id(node), node))
                continue
            
            # process valid best node
            if not current_covered_edges:
                continue

            set_u = set([i[0] for i in current_covered_edges])
            set_v = set([i[1] for i in current_covered_edges])
            result_set.append((set_u, set_v))
            result += len(set_u) + len(set_v)
            
            # remove covered edges from B
            edges_to_remove = set(current_covered_edges)
            B_int = B_int - edges_to_remove
            
            # re-evaluate node
            if len(B_int) > 0:
                 new_weight, _ = self._calculate_node_weight(node, points_index, B_int, self.n)
                 if new_weight != float("inf"):
                     heapq.heappush(pq, (new_weight, id(node), node))

        return result, result_set

    def _calculate_node_weight(
        self,
        node: RTreeNode,
        points_index: dict,
        B: set[Tuple],
        num_points_v: int,
    ) -> Tuple[float, List[Tuple[int, int]]]:

        if not hasattr(node, "int_points"):
            node.int_points = [points_index[pt] for pt in node.points_v] + [
                points_index[pt] for pt in node.points_u
            ]

            node.active_int_edges = []
            for u, v in node.previous_uncovered_edges:
                node.active_int_edges.append((points_index[u], points_index[v]))

        points = node.int_points

        # filter edges
        active_edges = []
        for edge in node.active_int_edges:
            if edge in B:
                active_edges.append(edge)

        # use cache if no change in active edges
        if len(active_edges) == len(node.active_int_edges) and not node.initial_call:
            return node.previous_density, node.previous_covered_edge_list

        # update the node's state
        node.active_int_edges = active_edges
        
        # GWBC: use exact algo        
        if self.c > 0:
            subgraph = Find_Densest_Subgraph(nodes=points, edges=active_edges, c=self.c)
        # WBC: use 2-approximation
        else:
            subgraph = two_approx_densest_subgraph(nodes=points, edges=active_edges)
        
        covered_edges = 0
        covered_edges_list = []

        subgraph_set = set(subgraph) if not isinstance(subgraph, set) else subgraph

        # get covered edges count and list
        for u, v in active_edges:
            if u in subgraph_set and v in subgraph_set:
                covered_edges += 1
                covered_edges_list.append((u, v))
        
        node.initial_call = False

        if covered_edges == 0:
            node.previous_density = float("inf")
            node.previous_covered_edge_list = []
            return float("inf"), covered_edges_list
        else:
            density = (len(subgraph) + self.c) / covered_edges
            node.previous_density = density
            node.previous_covered_edge_list = covered_edges_list
            return density, covered_edges_list
