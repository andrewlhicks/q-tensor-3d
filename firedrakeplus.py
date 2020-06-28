from firedrake import *

def errorH1(func_comp,func_true):
    diff = func_comp - func_true
    return sqrt(assemble((inner(grad(diff),grad(diff)) + dot(diff,diff)) * dx))

def errorL2(func_comp,func_true):
    diff = func_comp - func_true
    return sqrt(assemble((dot(diff,diff)) * dx))

# END OF CODE