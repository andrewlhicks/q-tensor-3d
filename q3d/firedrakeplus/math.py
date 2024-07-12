from firedrake import Function
from firedrake import sqrt, assemble, inner, grad, dot
from firedrake import dx

__all__ = ('nrm','errorH1','errorL2')

class nrm:
    def inf(function): # REMOVE?
        abs_function = Function(function._function_space).interpolate(abs(function)) # CHANGE, this code seems sketchy
        with abs_function.dat.vec_ro as v:
            norm = v.max()[1]
        return norm
    def H1(function,mesh):
        return sqrt(assemble((inner(grad(function),grad(function)) + dot(function,function)) * dx(domain=mesh)))
    def L2(function,mesh):
        return sqrt(assemble(dot(function,function) * dx(domain=mesh)))

def errorH1(func_comp,func_true,mesh):
    return nrm.H1(func_comp - func_true,mesh)

def errorL2(func_comp,func_true,mesh):
    return nrm.L2(func_comp - func_true,mesh)