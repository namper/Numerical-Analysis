# heat equation approximation
# using crank nicholson scheme for parabolic equations
import numpy as np
from heat_equation import u_zero, miu_1, miu_2, f, u_exact


def factorization_method(tau: float, h: float, t_0: int = 0, t_1: int = 2, x_0: int = -2,
                         x_f: int = 2):
    assert t_1 > t_0
    assert x_f > x_0

    # grid params
    m = int((t_1 - t_0) // tau)
    n = int((x_f - x_0) // h)

    print(f'shape >> m: {m}, n:{n}')
    # initiating factorization constants
    gamma = tau / h ** 2
    c = 1 + 2 * gamma

    # boundary conditions
    y = np.zeros(shape=(m, n))
    y.T[0] = miu_1(x_0, m)
    y.T[n - 1] = miu_2(x_f, m)
    y[0] = u_zero(t_0, n)

    # computing alpha vector
    alpha = np.zeros(shape=n)
    alpha[0] = 0
    for i in range(1, n):
        alpha[i] = gamma / (c - gamma * alpha[i - 1])

    beta = np.zeros(shape=n)

    for j in range(1, m):
        # computing beta given j-th layer
        beta[0] = y[j][n - 1]
        for i in range(1, n):
            _F = -y[j][i - 1] - tau * f(x_0 + (i - 1) * h, t_0 + (j + 1) * h)
            beta[i] = (gamma * beta[i - 1] - _F) / (c - gamma * alpha[i - 1])

        for i in range(n - 2, 0, -1):
            y[j][i] = alpha[i + 1] * y[j][i + 1] + beta[i + 1]

    return y


if __name__ == '__main__':
    approx = factorization_method(tau=10 ** (-2), h=10 ** (-3))
    from utils import check_error
    print(check_error(approx, u_exact))
