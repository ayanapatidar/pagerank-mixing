import numpy as np


def tv(x, y):
    # Total variation distance between two distributions
    return 0.5 * np.abs(x - y).sum()

