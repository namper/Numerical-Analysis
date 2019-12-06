"""
    Heat Equation with given parameters
    x: {-2, 2} and t: {0, 2}
    u(x, t) := 1/(x^2 + t^2)
    Attempting to solve given heat eq:
    du/dt = du^2/dx^2 + 2(-t^2 + t^3 + 3x^2 +tx^2)/(t^2 + x^2+1)^3
    assumptions:
        f(x_i, t_{j+1}) = y_i^{j+1} ~ y(x_i, t_{j+1}) |   t < 1/2 * h^2
"""

import numpy as np
from heat_equation import u_zero, miu_1, miu_2, f, u_exact


def approx_differential(tau: float, h: float, t_0: int = 0, t_1: int = 2, x_0: int = -2, x_f: int = 2):
    assert t_1 > t_0
    assert x_f > x_0
    assert tau / h ** 2 <= 0.5

    m = int((t_1 - t_0) // tau)
    n = int((x_f - x_0) // h)

    alpha = tau / h ** 2
    beta = 1 - 2 * alpha

    y = np.zeros(shape=(m, n))

    # boundary conditions
    y.T[0] = miu_1(x_0, m)
    y.T[n - 1] = miu_2(x_f, m)
    y[0] = u_zero(t_0, n)

    for i in range(1, n - 1):
        for j in range(m - 2):
            y[j + 1][i] = alpha * y[j + 1][i + 1] + beta * y[j][i - 1] + tau * f(i, j)
    return y


if __name__ == '__main__':
    approximated = approx_differential(tau=0.001, h=0.1)
    from utils import check_error
    print(check_error(approximated, u_exact))
