from firedrake import op2
from firedrake.slate.slac.compiler import PETSC_ARCH

# Kernel string for eigen decomposition

with open('eigen.cpp','r') as file:
    kernel_str = file.read()

# Kernel for eigen decomposition

kernel = op2.Kernel(kernel_str, 'get_reordered_eigendecomposition', cpp=True, include_dirs= ["%s/include/eigen3" % PETSC_ARCH])

# END OF CODE