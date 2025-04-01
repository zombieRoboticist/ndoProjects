import numpy as np;


def zoom(a0,ahi,phi, phiPrime, c1, c2, maxIterations=1000000):
	"""BFGS line search algorithm to find step length

	Inputs:
		a0: initial value for step size
		ahi: upper bound value for step size
		phi: scalar function for function vallue in a direction
		phiprime: scalar derivative of phi
		c1,c2: (opt) tuning values for line search (default: c1=.0001 (sufficient flattness) c2=.9 (sufficient decrease))
		maxIterations: (opt) max number of iterations to run (default=100000)
	Output: 
		step length"""
	# print(a0)
	# print(ahi)
	alo=a0;
	count=0;
	aj=0;
	# print("eval0: "+str(phi(0)))
	# print("deval0: "+str(phiPrime(0)))
	# print("eval1: "+str(phi(ahi)))
	# print("deval1: "+str(phiPrime(ahi)))
	# while number of itterations < some max number (to prevent infinite loops)
	while count<maxIterations:
		count+=1;
		# interpolate aj
		aj=(ahi+alo)/2
		# print(ahi)
		# print(aj)
		# print(alo)
		# print("evalj: "+str(phi(aj)))
		# print("devalj: "+str(phiPrime(aj)))
		# check sufficient decrease condition
		if(( phi(aj)>phi(0)+c1*(aj)*phiPrime(0) )or (phi(aj)>= phi(alo))):
			# if not viable move upper bound
			ahi=aj;
			
		else:
			# evaluate flatness condition
			if abs(phiPrime(aj))<=-c2*phiPrime(0):
				# if flat enough return aj
				return aj;
			if phiPrime(aj)*(ahi-alo)>=0:
				# enforce phi(ahi)>phi(alo)
				ahi=alo;
				# print(ahi)
			# if sufficient decrease satisfied but not flatness, update alo
			alo=aj;
			# print(alo)
		# if count>10:
		# 	return
	# if viable point not found print warning and return Null
	print("Warning: max number of itterations reached -z"+str(aj));
	return a0;