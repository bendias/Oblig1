# Oblig1
Oblig1 i INF5620
import numpy as np
import matplotlib.pyplot as plt

def solve(beta,theta,epsilon,num_periods,time_steps_per_period,plot = True):


### Setting the coefficients for the numerical scheme:

	theta = theta* np.pi/float(180) # converting the angle to radians
	P = 2*np.pi 			# Scaled period = 2pi

	T = P * num_periods 				# The total time should be number of periods * time per period( defined as P )
	n = int(time_steps_per_period * num_periods) 	# total number of steps = times per period * number of periods
	dt = T/float(n)					# The length of one step is the total time divided by how many steps we take total.

	t = np.linspace(0,T,n+1) # filling in the time_array
	x = np.zeros_like(t)     # Empty array to fill in x-components
	y = np.zeros_like(t)	 # Empty array to fill in y-components
	Theta = np.zeros_like(t) # Empty array to fill in Theta-components


### Setting the initial conditions:
	
	x[0] = (1-epsilon)*np.sin(theta)	# Inital condition for x0
	y[0] = 1 - (1 + epsilon)*np.cos(theta)  # Inital condition for y0

	L = np.sqrt(x[0]**2 + (y[0] - 1)**2)	# The L_bar-value for finding x1 and y1

	x[1] = .5*(   2*x[0] - dt**2 * beta/float(1 - beta) * ( 1 - beta/float(L) ) * x[0]   )				# Aproximated x1
	y[1] = .5*(   2*y[0] - dt**2 * beta/float(1 - beta) * ( 1 - beta/float(L) ) * (y[0] - 1) - (dt)**2*beta   )	# Aproximated y1


	for i in range(1,n):
		L = np.sqrt(x[i]**2 + (y[i] - 1)**2)		# The L_bar-value for finding x_n and y_n, where N != {0,1}

		x[i+1] = 2*x[i] - x[i-1] - (dt)**2 * beta/float(1 - beta) * ( 1 - beta/float(L) ) * x[i] 			#Scheme
		y[i+1] = 2*y[i] - y[i-1] - (dt)**2 * beta/float(1 - beta) * ( 1 - beta/float(L) ) * ( y[i] - 1 ) - (dt)**2*beta #Scheme



	for i in range(n+1):
		Theta[i] = np.arctan(x[i]/float(1.0-y[i])) # Filling in the Theta-values, since we now have the x- and y-components

###
### exact solution for non-elastic pendulum, since we want to compare it with Thata-array when angle is less than 10 degrees.
###

	def scaledclassicalnonelastic(theta,t):
		scne = np.zeros_like(t)
		for i in range(len(t)):
			scne[i] = theta * np.cos(t[i])
		return scne



	if (plot == True): # If plot == True, then plot!

		plt.plot(x,y,'k')				# Plotting x vs y
		plt.title('Simulation of elastic pendulum')
		plt.xlabel('x')
		plt.ylabel('y')
		plt.legend(['beta = %g and epsilon = %g'%(beta,epsilon)])		
		#plt.gca().set_aspect('equal')
		plt.show()
		

		if (theta < 10 * np.pi / float(180)): 	          # Theta less than 10 degrees: compare it with exact solution!
			scne = scaledclassicalnonelastic(theta,t) # The exact solution
			plt.plot(t,Theta,'k')			  # Plotting Theta as a function of time
			plt.plot(t,scne,'r')			  # Comparing with exact solution
			plt.title('How theta varies with time')
			plt.xlabel('time')
			plt.ylabel('Angle')	
			plt.legend(['Theta calculated from the numerical scheme','Scaled classical non-elastic pendulum'])
			plt.show()
		else:
			plt.plot(t,Theta,'k')			# Plotting Theta as a function of time, with no comparing.
			plt.title('How Theta varies with time')
			plt.xlabel('time')
			plt.ylabel('Angle')	
			plt.legend(['Theta calculated from the numerical scheme'])
			plt.show()


		return x,y,Theta,t
	else:
		return x,y,Theta,t


if __name__ == '__main__': # Don't want to run this test function if we're importing this program.


#
# Defining exact solutions to compare with the numerical scheme in the test
#



	def exact_vibration(t,eps,beta):
		f_ = np.zeros_like(t)
		for i in range(len(t)):
			f_[i] = -epsilon*( np.cos(np.sqrt(beta/float(1-beta)) * t[i]))
		return f_

#
# Defining constants
#
	beta = .9
	epsilon = 1e-1
	tol = 1e-15
	tol2= 1e-1 # expection big error in c)
# b)
# Testing if zero angle & zero extension gives zero motion
#

	x,y,Theta,t = solve(beta,0,0,6,60,False)

	for i in range(len(t)):
		assert abs(x[i]) < tol,'bug in x-array'
		assert abs(y[i]) < tol,'bug in y-array'
		assert abs(Theta[i]) < tol,'bug in Theta-array'
#c)
# Testing if zero angle gives zero motion in x-direction
# and if the vibration holds the correct frequency and amplitude



	x,y,Theta,t = solve(beta,0,epsilon,6,60,False)

	for i in range(len(t)):
		assert abs(x[i]) < tol,'bug in x-array'
		assert abs(Theta[i]) < tol,'bug in Theta-array'


	f_ = exact_vibration(t,epsilon,beta)
	for i in range(len(t)):
		assert abs(y[i] - f_[i]) < tol2, 'bug in vibration equation'
	plt.plot(t,y,'k')
	plt.plot(t,f_,'r--')
	plt.title('Pure Vertical Motion')
	plt.xlabel('time')
	plt.ylabel('y')	
	plt.legend(['y computed from numerical scheme','exact discrete solution'])
	plt.show()


# d)

def demo(beta,theta):
	solve(beta,theta,0,3,600,True)
demo(.9999,15)

