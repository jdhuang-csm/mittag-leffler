import numpy as np
from scipy.special import gamma
from scipy.signal import residue

# Global Pade approximation of Miffag-Leffler function as described in https://arxiv.org/abs/1912.10996
# Valid for z < 0, 0 < alpha < 1, beta >= alpha, alpha != beta != 1
# Implemented by Jake Huang

def solve_poly_coefs(alpha, beta, m=7, n=2):
	"""
	Solve for polynomial coefficients for rational expression.
	Only implemented for (m,n) = (7,2), which was shown to yield best results.
	:param float alpha: alpha parameter of ML function
	:param float beta: beta parameter of ML function
	:param int m: m value for approximation (7)
	:param int n: n value for approximation (2)
	"""
	# Perform checks
	check_ab(alpha, beta)
	check_mn(m, n)

	if beta > alpha:
		A = np.zeros((7, 7))
		np.fill_diagonal(A[:3, :3], [1, 1, 1])
		A[6, 2] = 1

		np.fill_diagonal(A[:4, 3:], -gamma(beta - alpha) / gamma(beta))
		np.fill_diagonal(A[1:5, 3:], gamma(beta - alpha) / gamma(beta + alpha))
		np.fill_diagonal(A[2:6, 3:], -gamma(beta - alpha) / gamma(beta + 2 * alpha))
		np.fill_diagonal(A[3:6, 3:6], gamma(beta - alpha) / gamma(beta + 3 * alpha))
		np.fill_diagonal(A[4:6, 3:5], -gamma(beta - alpha) / gamma(beta + 4 * alpha))
		A[5, 3] = gamma(beta - alpha) / gamma(beta + 5 * alpha)
		A[6, 6] = -1

		y = np.array([0, 0, 0, -1,
					  gamma(beta - alpha) / gamma(beta),
					  -gamma(beta - alpha) / gamma(beta + alpha),
					  -gamma(beta - alpha) / gamma(beta - 2 * alpha)
					  ])

		x = np.linalg.solve(A, y)
	
	else:
		raise ValueError('Not implemented for alpha==beta')

	return x
	
def get_partial_frac(coefs, alpha, beta):
	"""
	Get partial fraction decomposition of rational expression
	:param coefs array: polynomial coefficients
	:param float alpha: alpha parameter of the ML function
	:param float beta: beta parameter of the ML function
	"""
	if beta > alpha:
		p1, p2, p3, q0, q1, q2, q3 = coefs
		r, p, k = residue([1, p3, p2, p1], [1, q3, q2, q1, q0])
	else:
		raise ValueError('Not implemented for alpha==beta')
		
	return r, p, k

def ml_pade_approx(z, alpha, beta, m=7, n=2, decompose=True):
	"""
	Evaluate the Pade approximation of the ML function
	:param float z: alpha parameter of ML function
	:param float alpha: alpha parameter of ML function
	:param float beta: beta parameter of ML function
	:param int m: m value for approximation (7)
	:param int n: n value for approximation (2)
	:param decompose bool: if True, use the partial fraction decomposition (exactly equal with lower computation time). 
	If False, use the rational expression
	"""
	coefs = solve_poly_coefs(alpha, beta, m, n)

	if beta > alpha:
		func = create_approx_func(alpha, beta, m, n, decompose)
		out = func(z)
	else:
		raise ValueError('Not implemented for alpha==beta')

	return out
	
def create_approx_func(alpha, beta, m=7, n=2, decompose=True):
	"""
	Create a function to evaluate the ML approximation for fixed alpha, beta, m, n
	:param float alpha: alpha parameter of ML function
	:param float beta: beta parameter of ML function
	:param int m: m value for approximation (7)
	:param int n: 2 value for approximation (2)
	:param decompose bool: if True, use the partial fraction decomposition (exactly equal with lower computation time). 
	If False, use the rational expression
	"""
	coefs = solve_poly_coefs(alpha, beta, m, n)
	if beta > alpha:
		if decompose:
			r, p, k = get_partial_frac(coefs, alpha, beta)
			def approx_func(z):
				check_z(z)
				return 2 * (np.real(r[0] / (-z - p[0])) + np.real(r[2] / (-z - p[2])))
		else:
			p1, p2, p3, q0, q1, q2, q3 = coefs
			def approx_func(z):
				check_z(z)
				return (1 / gamma(beta - alpha)) * (p1 + p2*(-z) + p3*(-z)**2 + (-z)**3) / (q0 + q1*(-z) + q2*(-z)**2 + q3*(-z)**3 + (-z)**4)
			
	return approx_func
	

# Checks

def check_z(z):
	if np.max(z) >= 0:
		raise ValueError('Approximation is only valid for z < 0')
		

def check_ab(alpha, beta):
	if beta < alpha:
		raise ValueError('Approximation is only valid for beta >= alpha')
	elif (0 < alpha < 1) == False:
		raise ValueError('Approximation is only valid for 0 < alpha < 1')
	elif alpha == 1 or beta == 1:
		raise ValueError('Approximation is not valid if alpha = 1 or beta = 1')


def check_mn(m, n):
	if m != 7 or n != 2:
		raise ValueError('Only implemented for (m,n) = (7,2)')
	