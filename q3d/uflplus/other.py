from ufl import conditional, exp, gt
from q3d.uflplus.overwrites import conditional
from q3d.uflplus.domain import spherical_compatible

@spherical_compatible
def smooth_transition(x, *, I: list[int | float]):
    """ A smooth transition function in the variable x on the
    interval I = [a,b], with the value of 0 at x=a and 1 at x=b
    Arguments:
        x: The variable for which the bump transition function is evaluated.
        I: The interval [a,b] where the function transitions from 0 to 1.
    """

    a = I[0]
    b = I[1]

    # smooth function f defined by the exponential
    def f(x):
        return conditional(gt(x,0), exp(-1/x), 0)
    # smooth function g such that g(x)=0 for x<=0 and g(x)=1 for x>=1
    def g(x):
        return f(x)/(f(x)+f(1-x))

    return g((x-a)/(b-a))

def reverse_smooth_transition(x, *, I: list[int | float]):
    """ A smooth transition function in the variable x on the
    interval I = [a,b], with the value of 1 at x=a and 0 at x=b.
    This is the reverse of the smooth_transition() function.
    Arguments:
        x: The variable for which the bump transition function is evaluated.
        I: The interval [a,b] where the function transitions from 1 to 0.
    """
    return 1 - smooth_transition(x, I=I)

def bump_transition(x, *, I1: list[int | float], I2: list[int | float]):
    """ A smooth transition function in the variable x on the
    intervals I1 = [a,b] and I2 = [c,d], with the value of 0 at x=a and 1 at x=b,
    and the value of 1 at x=c and 0 at x=d. Note that a < c <= b < d.
    Within [b,c], evaluates to 1, and outside [a,d], evaluates to 0.

    Credit to Yuzhi Liu for the idea of this function.

    Arguments:
        x: The variable for which the bump transition function is evaluated.
        I1: The first interval [a,b] where the function transitions from 0 to 1.
        I2: The second interval [c,d] where the function transitions from 1 to 0.
    """

    return 1 - reverse_smooth_transition(x, I=I1) - smooth_transition(x, I=I2)

@spherical_compatible
def gaussian(x, *, mu=0, sigma=1):
    """ Gaussian function in the variable x with mean mu and standard deviation sigma.
     Arguments:
        x: The variable for which the Gaussian function is evaluated.
        mu: The mean of the Gaussian distribution.
        sigma: The standard deviation of the Gaussian distribution."""

    return exp(-(x - mu)**2 / (2 * sigma**2))