import heat 
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


def plot_exact(x_lin, t_lin, plt_num, exact_func):
	x = np.linspace(*x_lin, plt_num)
	t = np.linspace(*t_lin, plt_num)
	X, T = np.meshgrid(x, t)
	Z = exact_func(X, T)
	print(Z.shape)
	fig = plt.figure()
	ax = plt.axes(projection='3d')
	ax.plot_surface(X, T, Z, rstride=1, cstride=1,cmap='viridis', edgecolor='none')
	plt.show()


if __name__ == '__main__':
	plot_exact(
		x_lin= (-2, 2),
		t_lin = (0, 2),
		plt_num = 30,
		exact_func = heat.u_exact)
