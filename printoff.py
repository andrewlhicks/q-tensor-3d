from firedrake import *
from misc import color

def initPrintoff(L1,L2,L3,A,B,C,ep,ksp_type,pc_type,dt,end,visualize,manufactured,meshsize,init_meshsize,max_meshsize):
    print()
    print(f"{color.uline}PRELIMINARY INFO:{color.end}")
    print()
    
    if 0>=L1 or -L1>=L3 or L3>=2*L1 or -3/5*L1-1/10*L3>=L2:
        print(f"{color.warning}WARNING: L1, L2, and L3 do not satisfy the proper inequalities{color.end}")
        print()
    
    print(f"Constants: L1 = {L1},")
    print(f"           L2 = {L2},")
    print(f"           L3 = {L3},")
    print(f"            A = {A},")
    print(f"            B = {B},")
    print(f"            C = {C}")
    print(f"      epsilon = {ep}")
    print()
    print(f"KSP type: \"{ksp_type}\"")
    print(f"PC type:  \"{pc_type}\"")
    print()
    print(f"Time step: {dt}")
    print(f"End time: {end}")
    print(f"No. time steps: {end/dt:0.0f}")
    print()
    
    if visualize == 0:
        print(f"Visualize in Paraview? No")
    elif visualize == 1:
        print(f"Visualize in Paraview? Yes")
    else:
        raise ValueError("Variable 'visualize' must be 0 or 1.")
    print()
    
    if manufactured == 0:
        print(f"Manufactured solution? No")
        print()
        print(f"Mesh size: {meshsize}")
    elif manufactured == 1:
        print(f"Manufactured solution? Yes")
        print()
        print(f"Init mesh size: {init_meshsize} x {init_meshsize} x {init_meshsize}")
        print(f"Max mesh size:  {max_meshsize} x {max_meshsize} x {max_meshsize}")
    else:
        raise ValueError("Variable 'manufactured' must be 0 or 1.")
    print()
    
    print(f"{color.uline}ERROR CALCULATIONS:{color.end}")
    print()

def errorPrintoff(meshsize,q_soln,g,calctime):
    print(f"Mesh size:    {meshsize} x {meshsize} x {meshsize}")
    print(f"L2 error:     {sqrt(assemble(dot(q_soln - g, q_soln - g) * dx)):0.15f}")
    print(f"H1 error:     {sqrt(assemble( (inner(grad(q_soln-g),grad(q_soln-g)) + dot(q_soln - g, q_soln - g) ) * dx)):0.15f}")
    print(f"Time elapsed: {calctime:0.2f} seconds")
    print()