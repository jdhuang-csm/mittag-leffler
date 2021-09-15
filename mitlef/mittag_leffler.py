import numpy as np
from scipy.special import gamma
from .ml_internal import LTInversion

# Python port of MATLAB implementation of the Mittag-Leffler function
# Written by Konrad Hinsen

def ml(z, alpha, beta=1., gama=1.):
	"""
	Calculate Mittag-Leffler function
	
	:param array z: real or complex input
	:param array times: array of measurement times
	"""
	eps = np.finfo(np.float64).eps
	if np.real(alpha) <= 0 or np.real(gama) <= 0 or np.imag(alpha) != 0. \
	   or np.imag(beta) != 0. or np.imag(gama) != 0.:
		raise ValueError('ALPHA and GAMA must be real and positive. BETA must be real.')
	if np.abs(gama-1.) > eps:
		if alpha > 1.:
			raise ValueError('GAMMA != 1 requires 0 < ALPHA < 1')
		if (np.abs(np.angle(np.repeat(z, np.abs(z) > eps))) <= alpha*np.pi).any():
			raise ValueError('|Arg(z)| <= alpha*pi')
	
	if(type(z)==complex):
		return np.vectorize(ml_, [np.complex128])(z, alpha, beta, gama)
	else:
		return np.vectorize(ml_, [np.float64])(z, alpha, beta, gama)

def ml_(z, alpha, beta, gama):
	# Target precision 
	log_epsilon = np.log(1.e-15)
	# Inversion of the LT
	if np.abs(z) < 1.e-15:
		return 1/gamma(beta)
	else:
		return LTInversion(1, z, alpha, beta, gama, log_epsilon)
		



