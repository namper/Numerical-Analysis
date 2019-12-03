# heat equation approximation
# using crank nicolson's scheme for parabolic equations
# given du/dx = d^2u/dt^2 + f(t)
#t/2h^2u_i+1^k - (t/h^2 + 1)u_i^k + t/2h^2^u_{i-1}^k = phi_i^{k-1}
# phi_i^{k-1} = - tau/(2h^2)u_{i+1}^{k-1} - (tau/h^2 - 1)u_i^k - t/2h^2 u_{i-1}^{k-1}- tau*f(x_i, t_{k-1/2})
# copyright (m.okropiridze)

import numpy as np
from heat import u_exact, u_zero, miu_1, miu_2, f
from typing import Union, Callable


def factorization_method(tau: float = 0.02, h: float = 0.25, t_0: float = 0, t_1: float = 2, x_0: float = -2, x_f: float =2): 
    assert t_1 > t_0 
    assert x_f > x_0

    # grid params
    m = int((t_1 - t_0)//tau)
    n = int((x_f - x_0)//h)

    # initiating factorization constants
    gamma = tau/h**2
    C = 1 + 2*gamma

    # boundary conditions
    y = np.zeros(shape=(m, n))
    y.T[0] = miu_1(x_0, m)
    y.T[n-1] = miu_2(x_f, m)
    y[0] = u_zero(t_0, n)

    # computing alpha vector
    alpha = np.zeros(shape=n)
    alpha[0] = 0
    for i in range(1, n):
        alpha[i] = gamma/(C-gamma*alpha[i-1])


    beta = np.zeros(shape=n)

    for j in range(1, m):
        # computing beta given j-th layer
        beta[0] = y[j][n-1]
        for i in range(1, n):
            F = -y[j][i-1] - tau*f(x_0+(i-1)*h, t_0 + (j+1)*h)
            beta[i] = (gamma*beta[i-1]-F)/(C-gamma*alpha[i-1])

        for i in range(n-2, 0, -1):
            y[j][i] = alpha[i+1]*y[j][i+1]+beta[i+1]

    return y

if __name__ == '__main__':
    approx = factorization_method()
