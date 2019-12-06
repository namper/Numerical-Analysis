from typing import Callable
from numpy import ndarray


def check_error(approximated: ndarray, exact_fun: Callable):
    average_er = 0
    n, m = approximated.shape
    for i, column in enumerate(approximated):
        for j, val in enumerate(column):
            average_er += abs(val - exact_fun(j, i))
    return average_er / (n * m)
