import networkx as nx
import numpy as np
import pytest
from scipy import sparse

from src.chain import build_chain, distances, stationary

ALPHAS = [0.5, 0.85, 0.99]


def random_graph(seed, n=50, p=0.08, n_dangling=5):
    rng = np.random.default_rng(seed)
    M = (rng.random((n, n)) < p).astype(float)
    np.fill_diagonal(M, 0)
    M[rng.choice(n, n_dangling, replace=False)] = 0
    return sparse.csr_matrix(M)


@pytest.mark.parametrize("alpha", ALPHAS)
def test_matches_networkx(alpha):
    A = random_graph(seed=0)
    PT, dangling = build_chain(A)
    pi, _ = stationary(PT, dangling, alpha)

    G = nx.from_scipy_sparse_array(A, create_using=nx.DiGraph)
    ranks = nx.pagerank(G, alpha=alpha, tol=1e-14, max_iter=100_000)
    reference = np.array([ranks[i] for i in range(A.shape[0])])

    assert pi == pytest.approx(reference, abs=1e-8)


@pytest.mark.parametrize("alpha", ALPHAS)
def test_contraction_bound(alpha):
    A = random_graph(seed=0)
    PT, dangling = build_chain(A)
    pi, _ = stationary(PT, dangling, alpha)

    start = np.zeros(A.shape[0])
    start[0] = 1.0
    d = distances(start, pi, PT, dangling, alpha)
    d = d[d >= 1e-8]  # stay well above rounding error
    ratio = d / (alpha ** np.arange(len(d)) * d[0])

    assert (ratio <= 1 + 1e-6).all()
    assert (np.diff(ratio) <= 1e-6).all()