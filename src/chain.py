import numpy as np
from scipy import sparse
import math
from .metrics import tv


def build_chain(A):
    # Return (PT, dangling) for the adjacency matrix A

    # PT: the transpose of the row-normalised link matrix (CSR)
    # dangling: a boolean array marking pages with no out-links
    
    out_degree = np.asarray(A.sum(axis=1)).ravel()
    dangling = out_degree == 0
    inverse_degree = np.zeros_like(out_degree)
    inverse_degree[~dangling] = 1.0 / out_degree[~dangling]
    P = sparse.diags(inverse_degree) @ A
    return P.T.tocsr(), dangling


def step(x, PT, dangling, alpha):
    n = len(x)
    dangling_mass = x[dangling].sum()
    x_new = alpha * (PT @ x + dangling_mass / n) + (1.0 - alpha) / n
    return x_new / x_new.sum()

def stationary(PT, dangling, alpha, precision=1e-15):
    # Return the stationary distribution and the steps used.
    n = PT.shape[0]
    steps = math.ceil(math.log(precision) / math.log(alpha))
    x = np.full(n, 1.0 / n)
    for _ in range(steps):
        x = step(x, PT, dangling, alpha)
    return x, steps


def distances(x0, pi, PT, dangling, alpha, floor=1e-12, max_steps=10_000):
    # Return the TV distances from pi along the trajectory started at x0.
    # Entry k is the distance after k steps. 
    x = x0.copy()
    out = [tv(x, pi)]
    while out[-1] >= floor and len(out) <= max_steps:
        x = step(x, PT, dangling, alpha)
        out.append(tv(x, pi))
    return np.array(out)