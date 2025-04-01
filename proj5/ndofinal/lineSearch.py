from zoom import zoom;
import numpy as np;

def lineSearch(amax,phi,phiprime,c1=0.0001,c2=.9,maxIterations=100000):
	"""BFGS line search algorithm to find step length

	Inputs:
		amax: maximum value for step size
		phi: scalar function for function vallue in a direction
		phiprime: scalar derivative of phi
		c1,c2: (opt) tuning values for line search (default: c1=.0001 (sufficient flattness) c2=.9 (sufficient decrease))
		maxIterations: (opt) max number of iterations to run (default=100000)
	Output: 
		step length"""
	ai=1.0
	# print(amax)
	# initialize variables
	ai2= 0.0;
	count=0;
	# print("eval0: "+str(phi(0)))
	# print("deval0: "+str(phiprime(0)))
	# print("eval1: "+str(phi(ai)))
	# print("deval1: "+str(phiprime(ai)))
	while count<maxIterations:
		
		# print(str(c1*phiprime(0)));
		# check sufficient decrease
		if (phi(ai)>(phi(0)+c1*ai*phiprime(0))) or ((phi(ai)>=phi(ai2)) and (count>0)):
			# print("zoooominggg"+str(count))
				# call zoom and return value
			return zoom(ai2,ai,phi,phiprime,c1,c2);
			#check curavture condition
		# print(abs(phiprime(ai)))
		# print(-1*c2*phiprime(0))
		# break
		if (abs(phiprime(ai))<=(-1*c2*phiprime(0))):
			# print("succssesss")
			#return current value
			return ai;
			#if phi increasing
		if phiprime(ai)>=0:
			#call zoom and return
			return zoom(ai, ai2,phi,phiprime,c1,c2);
		#reset ai2
		ai2=ai;
		count +=1;

		# interpolate ai
		ai*=1.1;

		# if viable point not found print warning and return Null
	print("Warning: max number of itterations reached"+ai);
	return None;