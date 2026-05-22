# Algorithm instruction and usage guide

## Overview

This repository compares our algorithm performance with the baselines with regard to the WBC and GWBC problems

## Set up:

This repository uses Python [3.12.9](https://www.python.org/downloads/release/python-3129/). Please also install [pip](https://packaging.python.org/en/latest/tutorials/installing-packages/) for package installation. [Virtual environment](https://docs.python.org/3/library/venv.html) is also highly suggested.

To install and activate virtual environment:

```
py -m venv venv
.\venv\Scripts\activate
```

After that, install required packages:

```
pip install -r requirements.txt
```

Torch with cuda installation for parallel optimization for points' embeddings:

```
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu126
```

Test cuda device by running

```python
import torch
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"Using device: {device}")
```

## AMBEA

AMBEA source code can be found [here](https://github.com/ISCS-ZJU/AMBEA) or under the folder [AMBEA](AMBEA). The sake of usage, we have made some minor changes to the source code so that the input and output are more friendly for our Python code. Details of the changes and differences can be found under the folder AMBEA. Recompilation might be tricky due to inconsistency in architechture and g++/gcc version. We have precompiled and included the executable [file](MBE_ALL.exe).

In high level, in the modified version, the input file, edges.txt, for the executable file contains $n_v$ rows, corresponding to vertices $0, ..., n_v-1$ of set $V$. Each row $i$ is a space-separated list of vertex index $j$ of set U if and only if $(i,j)$ is an edge. Then, the executable file ingests the file and outputs the output file bicliques.txt. Essentially, each row in the output file is a maximal biclique in the format $V'|U'$ where $V'$ is a subset of vertices in $V$, $U'$ is a subset of vertices in $U$, and $V', U'$ forms a valid maximal biclique. The output file is finally read by our Python program to get the set of candidates.

The modified version also allows timeout argument to return the current set of maximal bicliques if exceed that time. In the file [baseline.py](/algorithm/baseline.py), row 118-119:

```python
with open('bicliques.txt', 'w') as f:
        check_call(['MBE_ALL.exe', '-i', 'edges.txt', '-s', '5', '-T', str(timeout), '-S'], stdout=f, cwd='.')
```

The -i is the flag for input file, -s is the flag for algorithm (5 for AMBEA), -T is the flag for time for timeout, and -S is the flag for randomly shuffling the vertices for diverse output.

## CPGR

CPGR source code can be found [here](https://github.com/SanazRabinia1/CPGR-Paper/tree/master) or under folder [CPGR](CPGR).

In high level, the CPGR try to create tricliques based on a given value of delta between 0 and 1 (this is different from the 𝛿 in our paper). To fit the greedy fashion of the baselines and our algorithm, we calculate the lower bound of delta $LB_{\text{delta}}$ and iterate through every value of delta between $LB_{\text{delta}}$ and 1 with a default step of 0.05 to see which value returns the tricliques with the largest ratio $\frac{\text{number of edges cover}}{\text{number of nodes used by the tricliques}}$. The process is repeated until CPGR cannot find any tricliques. The result is calculated as the number of vertices used by the tricliques plus two times the number of remaining uncovered edges.

The size of the step can be modified in line 12 of [this file](algorithm/CPGR.py)

## Converting any bipartite graphs to geometric graphs

To translate an arbitrarily bipartite graph to 𝛿-disk bipartite graph in $\ell_2$, run [this file](scripts/bipartite_to_coordinates.py). The theoretical inspiration is the [sphericity of a graph](<https://en.wikipedia.org/wiki/Sphericity_(graph_theory)>), while finding the embedding is completed by running pytorch optimzation on the distance matrix.

Currently support any Matrix Market (mtx) file from this [website](https://sparse.tamu.edu/), [movieslen dataset 100K](https://grouplens.org/datasets/movielens/100k/) (or use directly [this file](data/u.data)), and [movieslen dataset 1M](https://grouplens.org/datasets/movielens/1m/) (or use directly [this file](data/ratings.dat)). Any additional dataset can be modified to work. Please refer to the [this file](scripts/bipartite_to_coordinates.py) logic.

In the file, please fill in the path of the original data file (input_path), the path of the to-be-converted CSV edge file (output_csv_path), and the path to the output embedding file (output_embeddings_path). In addition, to find a good embedding, please tune the hyperparameter of the function embed_coordinates.

To run the file:

```
cd scripts
python bipartite_to_coordinates.py
```

## Run the code

Our Python code can be run by executing

```
python -u main.py --dataset popsim --distance_metric inf --c 0 --delta_list 0.04 --algo ouralgo --v_min 4000 --v_max 6000 --u_min 4000 --u_max 6000
```

Or on real dataset:

```
python -u main.py --dataset real_bipartite --distance_metric inf --c 0 --delta_list 0.04 --algo ouralgo --real_bipartite_csv data/embeddings.csv
```

The flags in the arguments are:

- --dataset: Dataset to use for the experiment. Can currently only be adults, credits, gamma, popsim, or real bipartite dataset (contain embeddings of a real bipartite graph generated by [this file](scripts/bipartite_to_coordinates.py)).
- --distance_metric: Distance metric to use for graph construction (default: inf). Use 1 for L1, 2 for L2, etc., and inf for L-infinity.
- --c: 0 for WBC, above 0 for GWBC.
- --delta_list: List of delta values to test.
- --algo: Algorithm to run: baseline or ouralgo.
- --baseline_mode: Baseline mode to use (1-2). 1 is for baseline 1, 2 is for AMBEA baseline. No need to specify if algo is not baseline.
- --seed: Random seed for reproducibility. Default is 2026.
- --v_min: Minimum number of vertices in set $V$ to use (for sampling). Only applicable for the popsim dataset (default: 4000).
- --v_max: Maximum number of vertices in set $V$ to use (for sampling). Only applicable for the popsim dataset (default: 6000).
- --u_min: Minimum number of vertices in set $U$ to use (for sampling). Only applicable for the popsim dataset (default: 4000).
- --u_max: Maximum number of vertices in set $U$ to use (for sampling). Only applicable for the popsim dataset (default: 6000).
- --max_time: Maximum time (in seconds) to allow for the baseline algorithm to run before terminating it (default: 3600).
- --real_bipartite_csv: path to CSV of the embeddings of the real bipartite graph. Only apply if --dataset is set to real_bipartite.

Memory performance can also be recorded using psrecord:

```
psrecord "python -u main.py --dataset popsim --distance_metric inf --c 0 --delta_list 0.04 --algo ouralgo --v_min 4000 --v_max 6000 --u_min 4000 --u_max 6000" --interval 3 --include-children --include-io --log log_performance.txt > log.txt
```

where

- log_performance.txt contains the memory activity of running the code
- log.txt containts the output of the script in the terminal

Note:

1. [this file](data/popsim_1M.csv) might have different newline characters in different environment. Please make sure to use correct encoding to read this dataset.
