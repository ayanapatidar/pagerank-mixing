from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

DATA_PATH = Path(__file__).resolve().parent.parent/"data"/"web-Stanford.txt.gz"


def load_graph(path=DATA_PATH):
    # returns (A, node_ids, info) for the directed graph stored at `path`.

    # A: n x n CSR matrix with A[i, j] = 1 when page i links to page j.
    # node_ids[i]: the original SNAP id of the page in row i.
    # info: dict of counts recorded while loading.
    
    edges = pd.read_csv(
        path, sep="\t", comment="#", header=None,
        names=["src", "dst"], dtype=np.int64,
    )
    n_raw = len(edges)
    edges = edges.drop_duplicates()

    # Remap the original ids to 0..n-1 so that row i of A is one page.
    node_ids, inverse = np.unique(edges.to_numpy().ravel(), return_inverse=True)
    idx = inverse.reshape(-1, 2)
    src, dst = idx[:, 0], idx[:, 1]

    n = len(node_ids)
    A = sparse.csr_matrix(
        (np.ones(len(src), dtype=np.float64), (src, dst)), shape=(n, n)
    )

    out_degree = np.asarray(A.sum(axis=1)).ravel()
    info = {
        "nodes": n,
        "edges": int(A.nnz),
        "duplicate_edges_dropped": n_raw - len(edges),
        "self_loops": int((src == dst).sum()),
        "dangling": int((out_degree == 0).sum()),
    }
    return A, node_ids, info