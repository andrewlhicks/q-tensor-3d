from settings import *
from misc import color

def initPrintoff():
    if omit_initial_printoff == 1:
        return
    
    print()
    print(f"{color.uline}PRELIMINARY INFO:{color.end}")
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
        print("Visualize in Paraview? No")
    elif visualize == 1:
        print("Visualize in Paraview? Yes")
    
    print()
    
    if manufactured == 0:
        print("Manufactured solution? No")
        print()
        print(f"Mesh size: {meshsize_max} x {meshsize_max} x {meshsize_max}")
    elif manufactured == 1:
        print("Manufactured solution? Yes")
        print()
        print(f"Init mesh size: {meshsize_init} x {meshsize_init} x {meshsize_init}")
        print(f"Max mesh size:  {meshsize_max} x {meshsize_max} x {meshsize_max}")

def initPrintoff2():
    print()
    print(f"{color.uline}ERROR CALCULATIONS:{color.end}")
    print()

def summaryPrintoff(meshsize,H1_error,L2_error,time_elapsed):
    print(f"Mesh size:    {meshsize} x {meshsize} x {meshsize}")
    print(f"H1 error:     {H1_error:0.15f}")
    print(f"L2 error:     {L2_error:0.15f}")
    print(f"Time elapsed: {time_elapsed:0.2f} seconds")
    print()

# END OF CODE