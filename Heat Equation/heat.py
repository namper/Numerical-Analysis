"""
Heat Equation with given parameters
	x: {-2, 2} and t: {0, 2} 
	u(x, t) := 1/(x^2 + t^2)
Attempting to solve given heat eq:
	du/dt = du^2/dx^2 + 2(-t^2 + t^3 + 3x^2 +tx^2)/(t^2 + x^2+1)^3
	assumptions:
		f(x_i, t_{j+1}) = y_i^{j+1} ~ y(x_i, t_{j+1}) |   t < 1/2 * h^2 

Forensic Analysis: Mishiko Okropiridze
"""

import numpy as np

# building the program from exact solution artifically


def f(x: float, t: float) -> float:
	# NOTE: CHANGE THIS GIVEN YOUR FUNCTION
	return -2*(-t**2 + t**3 + 3*x**2 + t*x**2)/(t**2 + x**2+1)**3

def u_exact(x: float, t: float) -> float:
	# NOTE: CHANGE THIS GIVEN YOUR FUNCTION	
	return 1/(x**2 + t**2 + 1)

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



def approx_differential(tau: float, h: float, t_0: float =0, t_1: float = 2, x_0: float = -2, x_f: float = 2):
	assert t_1 > t_0
	assert x_f > x_0
	assert tau/h**2 <= 0.5

	alpha = tau/h**2
	m = int((t_1 - t_0)//tau)

	n = int((x_f - x_0)//h)

	beta = 1 - 2*alpha

	y = np.zeros(shape=(m, n))
	

	# boundary conditions
	y.T[0] = miu_1(x_0, m)
	y.T[n-1] = miu_2(x_f, m)
	y[0] = u_zero(t_0, n)


	for i in range(1, n-1):
		for j in range(m-2):
			y[j+1][i] = alpha * y[j+1][i+1] + beta*y[j][i-1] + tau*f(i, j)

	return y

if __name__ == '__main__':
	approximated = approx_differential(tau=0.02, h=0.25)

	average_er = 0 
	n,m = approximated.shape
	for i, column in enumerate(approximated):
		for j, val in enumerate(column):
			average_er += abs(val - u_exact(j, i))

	print(average_er/(n*m))
