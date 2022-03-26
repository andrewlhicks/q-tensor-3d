from firedrake import SpatialCoordinate
from firedrake import FacetNormal
from firedrake import interpolate, assemble
from firedrake import dx, ds
from ufl.operators import *
import printoff as pr

class linesearch:
    def ls(name,*args,**kwargs):
        names = ('backtrack','exact1','exact2','none')
        if name not in names:
            raise ValueError(f'Must choose from {", ".join(names)}')
        if name == 'backtrack':
            return linesearch.backtrack(*args,**kwargs)
        elif name == 'exact1':
            return linesearch.exact1(*args,**kwargs)
        elif name == 'exact2':
            return linesearch.exact2(*args,**kwargs)
        else:
            return linesearch.none(*args,**kwargs)
    def backtrack(q_prev,time_der,alpha):
        """ Given the previous guess, the time derivative, and the time step
        alpha, returns xi computed by backtracking. """
        from firedrakeplus.eqnglobals import EqnGlobals

        weak_boundary = EqnGlobals.weak_boundary
        H1_vec = q_prev.function_space()

        xi = 8*alpha # Initial guess for xi, doesn't necessarily have to be 8 times the time step

        while xi > 1.0e-8: # Break the loop when xi becomes less than order 8 in magnitude
            q_next = interpolate(q_prev + xi*time_der,H1_vec)
            if compute_energy(q_next) < compute_energy(q_prev):
                return xi
            xi /= 2

        return xi
    def exact1(q_prev,time_der,alpha):
        """ Given the previous guess, the time derivative, and the time step
        alpha, returns xi computed by exact line search using Newton's
        method. """

        from numpy import abs

        H1_vec = q_prev.function_space()

        xi = 0 # Initial guess for xi, should be 0 because of initial guess for Newton's method

        for _ in range(100):
            xi_prev = xi
            q_next = interpolate(q_prev+xi_prev*time_der,H1_vec)
            first_der = compute_energy(q_next,time_der,der=1)
            secnd_der = compute_energy(q_next,time_der,time_der,der=2)
            xi = xi_prev - first_der/secnd_der
            if abs(xi - xi_prev) < 1.0e-8: break
        print(xi)
        return xi
    def exact2(q_prev,time_der,alpha):
        """ Given the previous time, the time derivative, and the time step
        alpha, returns xi compute by exact line search using the fact that xi is
        the root of a polynomial. """

        import numpy as np
        from numpy.polynomial import Polynomial
        import plot

        H1_vec = q_prev.function_space()

        x = np.linspace(0,1,5)
        q_next = [interpolate(q_prev+float(x[ii])*time_der,H1_vec) for ii in range(5)]
        y = np.array([compute_energy(q_next[ii]) for ii in range(5)])
        x_min = x[np.argmin(y)]

        poly = Polynomial.fit(x,y,4)
        xi = poly.deriv().roots()
        xi = xi[np.isclose(xi.imag, 0)]
        E = poly(xi)
        xi_min = xi[np.argmin(E)].real

        q_next = interpolate(q_prev+float(xi_min)*time_der,H1_vec)

        diff = compute_energy(q_next) - np.amin(y)

        if diff > 0:
            pr.warning(f'exact2 ls polynomial error {diff}')
            pr.info(f'{x_min}')
            return float(x_min)

        pr.info(f'{xi_min}')
        return float(xi_min)

    def none(q_prev,time_der,alpha):
        return alpha

def compute_energy(*function,der=0):
    from firedrakeplus.eqnglobals import EqnGlobals

    # Check for correct values

    if len(function) not in [1,2,3]:
        raise ValueError('Must choose 1, 2, or 3 functions.')
    if der not in [0,1,2]:
        raise ValueError('Must choose der=0, 1, or 2.')

    weak_boundary = EqnGlobals.weak_boundary
    energies = EqnGlobals.energies[2] if der==2 else EqnGlobals.energies[1] if der==1 else EqnGlobals.energies[0]

    q = function[0]

    if len(function) == 2:
        p = function[1]
    if len(function) == 3:
        # See documentation for secondVariationalDerivative for explanation
        r = function[1] # Should allow for customization in the future
        p = function[2]

    mesh = q.function_space().mesh()
    x0, x1, x2 = SpatialCoordinate(mesh)
    nu = FacetNormal(mesh)

    f = eval(EqnGlobals.forcing_f)
    g = eval(EqnGlobals.forcing_g)

    # Assemble domain integral

    domain_assembly = 0
    for energy in energies['domain']:
        domain_assembly += eval(energy)
    domain_integral = assemble(domain_assembly*dx)

    # Assemble boundary integral

    measure = determine_measure(weak_boundary[1])
    if measure is not None:
        boundary_assembly = 0
        for energy in energies['boundary']:
            boundary_assembly += eval(energy)
        boundary_integral = assemble(boundary_assembly*measure)
    else:
        boundary_integral = 0

    # print(float(domain_integral + boundary_integral))
    return float(domain_integral + boundary_integral)

def determine_measure(boundary_indicator):
    if boundary_indicator == 'none':
        return None
    if boundary_indicator == 'all':
        return ds
    if isinstance(boundary_indicator,int):
        if boundary_indicator >=0:
            return ds(boundary_indicator)
        raise ValueError('Boundary integer specified must be positive.')
    raise ValueError('Boundary specified must be \'all\', \'none\', or a positive integer.')