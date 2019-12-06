import numpy as np


# building the program from exact solution artificially
def f(x: float, t: float) -> float:
    #  NOTE: CHANGE THIS GIVEN YOUR FUNCTION
    return -2 * (-t ** 2 + t ** 3 + 3 * x ** 2 + t * x ** 2) / (t ** 2 + x ** 2 + 1) ** 3


def u_exact(x: float, t: float) -> float:
    # NOTE: CHANGE THIS GIVEN YOUR FUNCTION
    return 1 / (x ** 2 + t ** 2 + 1)


def miu_1(x_initial: float, m: int):
    res = np.zeros(shape=m)
    for t in range(m):
        res[t] = u_exact(x_initial, t)
    return res


def miu_2(x_final: float, m: int):
    res = np.zeros(shape=m)
    for t in range(m):
        res[t] = u_exact(x_final, t)
    return res


def u_zero(t_initial: int, n: int):
    res = np.zeros(shape=n)
    for x in range(n):
        res[x] = u_exact(x, t_initial)
    return res
