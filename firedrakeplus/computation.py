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

        H1_vec = q_prev.function_space()

        xi = 8*alpha # Initial guess for xi, doesn't necessarily have to be 8 times the time step

        while xi > 1e-13: # Break the loop when xi becomes small enough for machine error to occur
            q_next = interpolate(q_prev + xi*time_der,H1_vec)
            if compute_energy(q_next) < compute_energy(q_prev):
                return xi
            xi /= 2

        # pr.green(f'> ξb = {xi}')
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
        
        # pr.green(f'> ξ1 = {xi}')
        return xi
    def exact2(q_prev, search_dir, *, i=0):
        """ Given the previous step, the search direction, returns ideal time step
        alpha computed using the fact that the double well is the minimum of a
        polynomial of degree four. """

        import numpy as np

        pr.info('starting exact2 polynomial ls')

        H1_vec = q_prev.function_space()

        # take five evenly spaced sample points between 0 and 1 and calculate energies at these points
        alphas_lin = np.linspace(0,1,5)
        energies_lin = [compute_interpolate_energy(q_prev + float(alpha)*search_dir, H1_vec) for alpha in alphas_lin]
        alpha_lin_min = alphas_lin[np.argmin(energies_lin)]
        energy_lin_min = np.amin(energies_lin)

        # take the critical points of a polynomial of degree 4 fitted to the previous alphas/energies and calculate energy at these points
        alphas_poly = critical_pts_of_poly(alphas_lin, energies_lin, 4)
        energies_poly = [compute_interpolate_energy(q_prev + float(alpha)*search_dir, H1_vec) for alpha in alphas_poly] # as opposed to the previous, wrong approach: ---> energies_poly = poly(alphas_poly)
        alpha_poly_min = alphas_poly[np.argmin(energies_poly)]
        energy_poly_min = np.amin(energies_poly)

        # if polynomial minimum energy is greater than the linspace minimum energy, give warning and return linspace minimum energy
        if energy_poly_min - energy_lin_min > 1e-12:
            pr.warning(f'exact2 ls polynomial error, energy increase of {energy_poly_min - energy_lin_min}')
            pr.info('falling back to alpha from linear space')
            pr.info(f'alpha = {alpha_lin_min}')
            alpha_min = alpha_lin_min
        else:
            # return the argmin of the polynomial energy critical points
            pr.info(f'alpha = {alpha_poly_min}')
            alpha_min = alpha_poly_min

        return float(alpha_min)

    def none(q_prev,time_der,alpha):
        return alpha

def compute_energy(*functions, der=0):
    from firedrakeplus.eqnglobals import EqnGlobals
    from config import settings

    # Check for correct values

    if len(functions) not in [1,2,3]:
        raise ValueError('Must choose 1, 2, or 3 functions.')
    if der not in [0,1,2]:
        raise ValueError('Must choose der=0, 1, or 2.')

    weak_boundary = settings.options.weak_boundary
    energies = EqnGlobals.energies[2] if der==2 else EqnGlobals.energies[1] if der==1 else EqnGlobals.energies[0]

    q = functions[0]

    if len(functions) == 2:
        p = functions[1]
    if len(functions) == 3:
        # See documentation for secondVariationalDerivative for explanation
        r = functions[1] # Should allow for customization in the future
        p = functions[2]

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

    measure = determine_measure(weak_boundary)
    if measure is not None:
        boundary_assembly = 0
        for energy in energies['boundary']:
            boundary_assembly += eval(energy)
        boundary_integral = assemble(boundary_assembly*measure)
    else:
        boundary_integral = 0

    # print(float(domain_integral + boundary_integral))
    return float(domain_integral + boundary_integral)

def compute_interpolate_energy(*args, **kwargs):
    """ Computes the energy, but also adds interpolation over
    a function space given by the last arg. """

    H1_vec = args[-1]
    functions = args[:-1]
    functions = [interpolate(function, H1_vec) for function in functions]
    return compute_energy(*functions, **kwargs)

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

def critical_pts_of_poly(x_values, y_values, degree):
    """ Returns the critical points of a polynomial fitted to the x and y values
    given, and of degree specified """

    from numpy.polynomial import Polynomial

    # create polynomial using the x values, y values, and degree
    poly = Polynomial.fit(x_values, y_values, degree)
    # find its critical points by finding the roots of the derivative
    critical_pts = poly.deriv().roots().real

    return critical_pts