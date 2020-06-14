# RA, 2020-06-13

"""
Simulated annealing for finding the minimum of a simple function
"""

import time
import numpy as np
import matplotlib.pyplot as plt

from itertools import count


def f(x):
    return (x ** 4) + 2 * (x ** 3) - (x ** 2) - x + 2


def main():
    rs = np.random.RandomState(1)

    ax1: plt.Axes
    (fig, ax1) = plt.subplots()

    (a, b) = (-2, 2)

    xx = np.linspace(a, b, 1001)
    ax1.plot(xx, f(xx))
    plt.ion()
    plt.pause(2)
    plt.show()

    n = 101  # Number of states
    T1 = 1e+2  # Initial temperature
    T0 = 1e-2  # Stop at this temperature
    # x-coordinate of states
    X = dict(enumerate(np.linspace(-2, 2, n)))
    # States
    J = set(X.keys())
    # Neighbors
    N = {i: [j for j in (i - 4, i - 2, i - 1, i + 1, i + 2, i + 4) if (j in J)] for i in J}
    # N = {i: list(rs.choice(list(J), 3)) for i in J}

    i = max(J)
    T = T1

    for k in count():
        # Progress [0, 1]
        p = (np.log(T0 / T) / np.log(T0 / T1))

        # Gibbs likelihood of state
        L = (lambda j: np.exp(-f(X[j]) / T))

        T *= 0.99

        if (T < T0):
            break

        ax1.scatter(X[i], f(X[i]), marker='o', s=32, c=p, vmin=0, vmax=1, cmap=plt.cm.get_cmap("cool"), zorder=k)
        if not (k % 1):
            ax1.set_title(F"Step {k}, temperature {T:.3}")
            plt.pause(0.01)

        # Proposal
        j = rs.choice(N[i])
        if (L(j) >= L(i)) or (rs.rand() <= (L(j) / L(i))):
            # Accept
            i = j

    plt.ioff()
    plt.show()


if __name__ == '__main__':
    main()


