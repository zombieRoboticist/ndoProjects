from lineSearch import lineSearch 
import numpy as np

def bfgs (funcs, x0,H0,e, maxIterations = 1000):
	"""runs the Broyden–Fletcher–Goldfarb–Shanno (BFGS) quazi-Newton nonlinear unconstrained n-dimensional optimzation algorithm

	Inputs:
		funcs: the function to optimize (scalar output, vector input)
		x0: the initial guess for optimal inputs
		H0: initial guess for hessian
		e: the error tollerance in gradient magnitude
		maxIterations: (optional) the maximum number of iterations to run (default = 1000)
	Outputs:
		xk: the optimized inputs"""
	# initialize variables
	h=np.array(H0);
	xk=np.array(x0);
	dim=max(np.shape(H0));
	fk=grad(funcs,xk, isfk=1);
	count=0;

	# print out headers
	print("Iterations:\t|grad|")

	# loop through until the norm of the gradient is less than e or the maximum number of iterations is reached
	while count<maxIterations and np.linalg.norm(fk)>e:
		#calculate search direction
		pk=np.array(1);
		if dim==1:
			pk=np.array(-1*h*fk);
			pk=pk[:,0];

		else:
			pk=np.zeros(dim);
			c=0;
			for i in -1*np.matmul(h,fk):
				pk[c]=i[0];
				c+=1;

		#create scalar function in direction pk
		phi=lambda a:funcs(xk+a*pk);

		#create anonomous function for getting the derivative
		der=lambda a: phi(complex(a,0.000000000000000000000000000001)).imag/0.000000000000000000000000000001
		
		#compute line search
		a=lineSearch(1,phi,der,maxIterations=1000000);
		
		#calc xk+1
		xk1=xk+a*pk;

		#calc sk
		sk=a*pk;
		#calc yk
		yk=grad(funcs,xk1)-grad(funcs,xk);
		#calc rho k
		rhok=1/((yk)@cols(sk))[0]

		#calc Hk+1
		#eqn 6.17
		if dim==1:
			h=((np.eye(dim)-rhok*(sk*np.transpose(yk))))*h*(np.eye(dim)-rhok*(yk*np.transpose(sk)))+rhok*(sk*np.transpose(sk));
		else:
			a1=rhok*mult(cols(sk),sk)
			# print(a1)
			a21=np.eye(dim)-rhok*mult(cols(sk),yk)
			a22=np.eye(dim)-rhok*mult(cols(yk),(sk))
			a2=a21@(h@a22)
			h=a2+a1;

		# iterate count
		count +=1;

		# update gradient
		fk=grad(funcs,xk1, isfk=1);

		# update position
		xk=xk1;

		# print step info
		print(str(count)+'\t\t'+str(np.linalg.norm(fk)))

	# print warning if maximum number of iterations are reached
	if count>=maxIterations:
		print("warning max maxIterations reached");

	# return final position
	return xk



def grad(func, x,step=0.000000000000000000000000000001,isfk=0):
	"""use complex step to calculate function gradient

	Inputs:
		func - the function to calculate gradient of
		x - the point to evaluate at
		step - (optional) the imaginary step size to use
	Output: 
		the complex step approximation of the gradient""" 

	# initialize complex inputs
	y=np.zeros(np.shape(x), dtype=complex)
	for i in range(np.size(y)):
		y[i]=x[i];

	# 	initialize outputs
	z=np.zeros(np.shape(y))

# iterate through the input directions and get derivative in each direction using complex step
	for i in range(np.size(z)):
		# keep indicies in correct bounds
		if i==np.size(x):
			break;
		# add complex step to correct input direction
		y[i] =complex(x[i],step)
		# calculate derivative in that direction
		z[i]=func(y).imag/step;
		# remove the complex step from the variable
		y[i] =x[i]
	# convert to column vector if needed
	if isfk:
		z=cols(z)
# return output
	return z

def cols(x):
	"""define a row numpy array as a column matrix. Note this was implamented due to native python methods/python libraries acting oddly

	Inputs:
		x - the array to transpose
	Output: 
		z - the transposed array"""
	y=np.zeros((np.size(x),1));
	cont=0;
	for i in x:
		y[cont,0]=i;
		cont+=1;
	return y;
def mult(x,y):
	"""define a row numpy array as a column matrix. Note this was implamented due to native python methods/python libraries acting oddly

	Inputs:
		x - the first array to multiply 
		y - the second array to multiply
	Output: 
		z - the matrix resulting from the multiplication"""
	z=np.zeros((np.shape(x)[0],np.shape(x)[0]))
	for i in range(np.size(x)):
		for j in range(np.size(x)):
			z[i,j]=x[i][0]*y[j]
		# print(z)
	return z

#------------- test casses ----------------------
if __name__=='__main__':
	# one dimentional quadratic
	func=lambda x: x[0]*x[0];
	h0=np.array(np.eye(1))*.1
	x0=10*np.ones(1);
	print("bfgs"+str(bfgs(func,x0,h0,.000001)))
	x0=-10*np.ones(1);
	print("bfgs"+str(bfgs(func,x0,h0,.000001)))

	# two dimensional polynomial
	func = lambda x:x[0]*x[0]*x[0]*x[0]+ x[1]*x[1]
	h0=np.array(np.eye(2))*.1;
	x0=10*np.ones(2);
	print("bfgs"+str(bfgs(func,x0,h0,.0000001)))
	x0=-10*np.ones(2);
	print("bfgs"+str(bfgs(func,x0,h0,.0000001)))

	#two dimensional polynomial non-zeros minimum/position
	func = lambda x: (x[0]-1)*(x[0]-1)*(x[0]-1)*(x[0]-1)+ (x[1]+2)*(x[1]+2)+1
	h0=np.array(np.eye(2))*.1;
	x0=4*np.ones(2);
	print("bfgs"+str(bfgs(func,x0,h0,.0000001)))
	x0=-7*np.ones(2);
	print("bfgs"+str(bfgs(func,x0,h0,.0000001)))

	# two dimensional rosenbrock function
	func = lambda x: (((1-x[0])*(1-x[0]))+ 100*((x[1]-(x[0]*x[0]))*(x[1]-(x[0]*x[0]))))
	h0=np.array(np.eye(2))*.05;
	x0=0*np.ones(2);
	print("bfgs"+str(bfgs(func,x0,h0,.0000001)))

	# two dimensional polynomial, multiple minima
	func = lambda x: (x[0]-1)*(x[0]-1)*(x[0]-1)*(x[0]-1)+ 45*(x[1]-3)*(x[1]-3)+(x[1]-3)*(x[1]-3)*(x[1]-3)*(x[1]-3)-17*(x[1]-3)*(x[1]-3)*(x[1]-3)
	h0=np.array(np.eye(2))*.1;
	x0=4*np.ones(2);
	print("bfgs"+str(bfgs(func,x0,h0,.0000001)))
	x0=-7*np.ones(2);
	print("bfgs"+str(bfgs(func,x0,h0,.0000001)))

	#three dimensional polynomial

	func = lambda x: (x[0]-1)*(x[0]-1)*(x[0]-1)*(x[0]-1)+ (x[1]+2)*(x[1]+2)+1+5*(x[2]-3)*(x[2]-3)+(x[2]-3)*(x[2]-3)*(x[2]-3)*(x[2]-3)
	h0=np.array(np.eye(3))*.1;
	x0=4*np.ones(3);
	print("bfgs"+str(bfgs(func,x0,h0,.0000001)))
	x0=-7*np.ones(3);
	print("bfgs"+str(bfgs(func,x0,h0,.0000001)))