from ufl import *
from ufl.classes import * # allows for evaluating UFL cache strings into UFL objects
from firedrake.ufl_expr import * # overwrites some of ufl.classes using Firedrake-specific UFL object defintinions. Needed for Firedrake solvers to work properly
from ufl.operators import ListTensor, Zero # used for type hints/checking for UFL objects

from q3d.uflplus.fy import qvectorfy, qtensorfy
from q3d.uflplus.other import smooth_transition