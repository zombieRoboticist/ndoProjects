import numpy as np


def exp(func,c,offset)
	"""exponential inequality constraint/barrier function

	Input:
		func: a function that inputs x and returns a constant
		c: constant multiplier to tune barrier
		offset: a numpy array of the x offset
	Output:
		barrier function for lower inequality constraint"""

	return lambda x: np.e**(-c*func(x-offset))

def ln(func,c,offset)
	"""logarithmic inequality constraint/barrier function

	Input:
		func: a function that inputs x and returns a constant
		c: constant multiplier to tune barrier = 1/mu
		offset: a numpy array of the x offset
	Output:
		barrier function for lower inequality constraint"""

	return lambda x:-1*c*np.ln(func(x-offset))

def quad_inequal(func,c,offset):
	"""Quadratic inequality constraint/barrier function

	Input:
		func: a function that inputs x and returns a constant
		c: constant multiplier to tune barrier =mu/2
		offset: a numpy array of the x offset
	Output:
		barrier function for lower inequality constraint"""

	return lambda x:c*min(func(x-offset),0)**2


def lowerBound(func,c,offset,type=0):
	"""lower boundry barrier function

	Input:
		func: a function that inputs x and returns a constant
		c: constant multiplier to tune barrier =mu/2
		offset: a numpy array of the x offset
		type: the type of barrier function to use (0=exponential, 1=natural log, 2= quadratic, else= func is barrier function)
	Output:
		barrier function for lower inequality constraint"""

	if type==0:
		return exp(func,c,offset)
	elif type==1:
		return ln(func,c,offset)
	elif type==2:
		return quad_inequal(func,c,offset)
	else:
		return lambda x:c*func(x-offset)

def upperBound(func,c,offset,type=0):
	"""upper boundry barrier function

	Input:
		func: a function that inputs x and returns a constant
		c: constant multiplier to tune barrier =mu/2
		offset: a numpy array of the x offset
		type: the type of barrier function to use (0=exponential, 1=natural log, 2= quadratic, else= func is barrier function)
	Output:
		barrier function for upper inequality constraint"""

	return lowerBound(lambda x:-1*func(x),c,offset,type=type)


def quad_equal(func,c,offset):
	"""Quadratic equality constraint function

	Input:
		func: a function that inputs x and returns a constant (trys to find (func(x)=0))
		c: constant multiplier to tune constraint =mu/2
		offset: a numpy array of the x offset
	Output:
		quadratic equality constraint"""

	return lambda x:c*func(x-offset)**2