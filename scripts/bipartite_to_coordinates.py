import pandas as pd
import numpy as np
from scipy.io import mmread
from scipy.spatial.distance import pdist, squareform
import torch
import torch.nn as nn
import torch.optim as optim

def get_edges_from_mtx(input_path, output_path):
    # Read the matrix as a (preferably sparse) SciPy matrix
    matrix = mmread(input_path)
    print(f"Matrix type: {type(matrix)}, shape: {matrix.shape}")

    # Use sparse COO representation to get non-zero indices efficiently
    # Fall back to numpy nonzero if mmread returned a dense array-like
    try:
        coo = matrix.tocoo()
        rows_idx = coo.row
        cols_idx = coo.col
    except Exception:
        arr = np.asarray(matrix)
        rows_idx, cols_idx = np.nonzero(arr)

    n_rows = matrix.shape[0]
    num_u = n_rows.shape[0]
    num_v = matrix.shape[1]

    # Preserve original mapping: source = column_index + n_rows, target = row_index
    sources = (cols_idx + n_rows).astype(int)
    targets = rows_idx.astype(int)

    df = pd.DataFrame({'source': sources, 'target': targets})

    # Remove duplicate edges if any (fast, vectorized)
    before = len(df)
    df = df.drop_duplicates()
    after = len(df)

    df.to_csv(output_path, index=False, header=False)
    print(f"Nonzeros (raw): {before}, Unique edges written: {after}")
    
    return num_u, num_v

def embed_coordinates(num_u, num_v, output_csv_path, output_embeddings_path, embedding_dim=2, margin=0.01, lambda_rep=10.0, epochs=1000, learning_rate=0.001):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")
    
    total_nodes = num_u + num_v
    R = 1.0


    # Data Loading & Index Shifting
    print(f"Loading edges from u.data via read_data()...")
    # df = read_data()
    if output_csv_path.endswith('.csv'):
        df = pd.read_csv(output_csv_path, header=None, names=['source', 'target'])
    else:
        df = pd.read_csv(output_csv_path, sep=',', header=None, names=['source', 'target', 'rating', 'timestamp'])[['source', 'target']]
        df['source'] = df['source'] - 1
        df['target'] = df['target'] - 1

    # Move edge tensors directly to the GPU
    edges_u = (torch.tensor(df['source'].values, dtype=torch.long)).to(device)
    edges_v = (torch.tensor(df['target'].values, dtype=torch.long)).to(device)

    # Initialize the mask directly on the GPU to save transfer time
    adj_mask = torch.zeros((total_nodes, total_nodes), dtype=torch.bool, device=device)
    adj_mask[edges_u, edges_v] = True
    adj_mask[edges_v, edges_u] = True
    adj_mask.fill_diagonal_(True)

    repulsion_mask = ~adj_mask

    print(f"Loaded {len(edges_u)} edges. Beginning optimization...")
    print(f"Number of non zero entries: {adj_mask.sum()}")
    print(f"Total Nodes: {total_nodes}")


    # Optimizer Setup
    embeddings = nn.Embedding(total_nodes, embedding_dim).to(device)
    nn.init.normal_(embeddings.weight, std=0.1)
    optimizer = optim.AdamW(embeddings.parameters(), lr=learning_rate, weight_decay=1e-4)

    scheduler = optim.lr_scheduler.StepLR(optimizer, step_size=500, gamma=0.5)

    # Train
    for epoch in range(epochs):
        optimizer.zero_grad()

        coords = embeddings.weight

        # torch.cdist will now execute entirely on the GPU
        dist_matrix = torch.cdist(coords, coords, p=2.0)

        edge_dists = dist_matrix[edges_u, edges_v]
        loss_attr = torch.clamp(edge_dists - R, min=0.0).mean()

        non_edge_dists = dist_matrix[repulsion_mask]
        loss_rep = torch.clamp((R + margin) - non_edge_dists, min=0.0).mean()

        loss = loss_attr + (lambda_rep * loss_rep)

        loss.backward()
        optimizer.step()

        scheduler.step()

        if epoch % 50 == 0:
            print(f"Epoch {epoch:3d} | Total Loss: {loss.item():.4f} | "
                f"Attr: {loss_attr.item():.4f} | Rep: {loss_rep.item():.4f}")

    # Save embeddings
    final_coordinates = embeddings.weight.detach().cpu().numpy()

    out_df = pd.DataFrame(final_coordinates, columns=[f'dim_{i+1}' for i in range(embedding_dim)])
    out_df['partition'] = ['U'] * num_u + ['V'] * num_v
    out_df['original_id'] = list(range(num_u)) + list(range(num_v))

    out_df.to_csv(output_embeddings_path, index=False)
    print(f"\nOptimization complete. Coordinates saved to '{output_embeddings_path}'.")

def verify(edges_path, embedding_path):
    global device
    global num_u # Need num_u from training script for correct offset

    if edges_path.endswith('.csv'):
        edges_df_pd = pd.read_csv(edges_path, header=None, names=['source', 'target'])
    else:
        edges_df_pd = pd.read_csv(edges_path, sep=',', header=None, names=['source', 'target', 'rating', 'timestamp'])[['source', 'target']]
        edges_df_pd['source'] = edges_df_pd['source'] - 1
        edges_df_pd['target'] = edges_df_pd['target'] - 1

    coords_df_pd = pd.read_csv(embedding_path)

    # Convert to PyTorch tensors and move to device
    edges_u_tens = torch.tensor(edges_df_pd['source'].values, dtype=torch.long, device=device)
    edges_v_tens = torch.tensor(edges_df_pd['target'].values, dtype=torch.long, device=device)

    # Rebuild coordinates with the same ordering used in test(): U first, then V.
    u_df = coords_df_pd[coords_df_pd['partition'] == 'U'].sort_values('original_id')
    v_df = coords_df_pd[coords_df_pd['partition'] == 'V'].sort_values('original_id')

    num_u = len(u_df)

    print(f"Number of U nodes: {num_u}, Number of V nodes: {len(v_df)}")
    dim_cols = [col for col in coords_df_pd.columns if col.startswith('dim_')]

    coords_u_tens = torch.tensor(u_df[dim_cols].values, dtype=torch.float, device=device)
    coords_v_tens = torch.tensor(v_df[dim_cols].values, dtype=torch.float, device=device)

    coords_all_tens = torch.cat((coords_u_tens, coords_v_tens), dim=0)

    edge_distances = torch.norm(coords_all_tens[edges_u_tens] - coords_all_tens[edges_v_tens], p=2, dim=1)
    max_r = edge_distances.max().item()

    print(f"Total edges: {len(edges_df_pd)}")
    print(f"Maximum distance: {max_r}")

    # Count all U-V pairs whose distance is within max_r.
    all_uv_distances = torch.cdist(coords_u_tens, coords_v_tens, p=2)

    num_actual_edges = (all_uv_distances <= max_r).sum().item()
    print(f"Number of actual edges: {num_actual_edges}")

if __name__ == "__main__":
    input_path = ""
    output_csv_path = ""
    output_embeddings_path = ""
    num_u, num_v = get_edges_from_mtx(input_path, output_csv_path)
    embed_coordinates(num_u, num_v, output_csv_path, output_embeddings_path)
    verify(output_csv_path, output_embeddings_path)