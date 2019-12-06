import numpy as np
import matplotlib.pyplot as plt
from typing import Callable
from heat_equation import u_exact
from mpl_toolkits import mplot3d


def plot_exact(x_lin, t_lin, plt_num, exact_func: Callable):
    x = np.linspace(*x_lin, plt_num)
    t = np.linspace(*t_lin, plt_num)
    plt_x, plt_t = np.meshgrid(x, t)
    plt_z = exact_func(plt_x, plt_t)
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_surface(plt_x, plt_t, plt_z, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
    plt.show()


if __name__ == '__main__':
    plot_exact(
        x_lin=(-4, 4),
        t_lin=(-2, 2),
        plt_num=30,
        exact_func=u_exact)
