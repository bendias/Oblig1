import sympy as sym


V, t, I, w, dt, b, a = sym.symbols('V t I w dt b a')  # global symbols
f = None  # global variable for the source term in the ODE

def ode_source_term(u):
    """Return the terms in the ODE that the source term
    must balance, here u'' + w**2*u.
    u is symbolic Python function of t."""
    return sym.diff(u(t), t, t) + w**2*u(t)

def residual_discrete_eq(u):
    """Return the residual of the discrete eq. with u inserted."""
   # u_kplus1 = 2*u(t) - u(t).subs(t,t-dt) + dt**2 * w**2*(V*t + I - u(t))
   # R = u(t).subs(t,t+dt) - u_kplus1
    R = sym.diff(u(t),t,t) - (u(t).subs(t,t+dt) - 2*u(t) + u(t).subs(t,t-dt))/dt**2
    return sym.simplify(R)

def residual_discrete_eq_step1(u):
    """Return the residual of the discrete eq. at the first
    step with u inserted."""

    u_1 = u(0) + dt*V + .5*dt**2  *( ode_source_term(u).subs(t,0) - w**2*u(0))
    R = u(t).subs(t,dt) - u_1
    return sym.simplify(R)
"""
    f = sym.diff(u(t),t,t).subs(t,0) + w**2*u(t).subs(t,0)  
    R = u(t).subs(t,dt) - (1-w**2/2)*u(0) + 1/2 * f
    return sym.simplify(R)
"""

def DtDt(u, dt):
    """Return 2nd-order finite difference for u_tt.
    u is a symbolic Python function of t.
    """
    return (u(t+dt) + u(t-dt) - 2*u(t))/dt**2

def main(u):
    """
    Given some chosen solution u (as a function of t, implemented
    as a Python function), use the method of manufactured solutions
    to compute the source term f, and check if u also solves
    the discrete equations.
    """
    print '=== Testing exact solution: %s ===' % u
    print "Initial conditions u(0)=%s, u'(0)=%s:" %(u(t).subs(t, 0), sym.diff(u(t), t).subs(t, 0))

    # Method of manufactured solution requires fitting f
    global f  # source term in the ODE
    f = sym.simplify(ode_source_term(u))

    # Residual in discrete equations (should be 0)
    print 'residual step1:', residual_discrete_eq_step1(u)
    print 'residual:', residual_discrete_eq(u)

def linear():
    main(lambda t: V*t + I)

def quadratic():
    main(lambda t: b*t**2 + V*t + I)

def cubic():
    main(lambda t: a*t**3 + b*t**2 + V*t + I)

if __name__ == '__main__':
    linear()
    quadratic()
    cubic()

def solver(n,T,I,V,w,plot = True):
	import numpy as np
	import matplotlib.pyplot as plt

	dt = 1/float(n)
	t = np.linspace(0,T,n+1)
	u = np.zeros_like(t)
	u[0] = I

	u[1] = u[0] + dt*V + .5 * dt**2 * w**2 * (I + V*t[0] - u[0])

	for i in range(1,n):
		u[i+1] = 2*u[i] - u[i-1] + dt**2 * w**2 *(I + V*t[i] - u[i])
	if (plot == True):
		plt.plot(t,u,'k')
		plt.title('numerical scheme to solve an ODE')
		plt.xlabel('t')
		plt.ylabel('u')
		plt.show()
	return t,u
solver(1000,1,1,1,1,True)


""" 

RUN EXAMPLE: 

=== Testing exact solution: <function <lambda> at 0x7f786a2fc398> ===
Initial conditions u(0)=I, u'(0)=V:
residual step1: 0
residual: 0
=== Testing exact solution: <function <lambda> at 0x7f786a2fc398> ===
Initial conditions u(0)=I, u'(0)=V:
residual step1: b*dt**2
residual: b*dt**2*(t**2*w**2 + 2)
=== Testing exact solution: <function <lambda> at 0x7f786a2c0758> ===
Initial conditions u(0)=I, u'(0)=V:
residual step1: dt**2*(a*dt + b)
residual: dt**2*(a*t**3*w**2 + 6*a*t + b*t**2*w**2 + 2*b)
"""
