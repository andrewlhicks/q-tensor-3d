""" This python program is meant to print off any relevant information to the
user in a visually pleasing way. For the initial printoff, we first print a
blank line to create visual separation between the regular command line and
the information printed off here. Then the program checks to see if the
initial printoff should be omitted or not. From this point forward, if a blank
line is needed for visual separation, it will be printed after the previous
print and not before. """

from settings import const, meshdata, options, visdata, timedata, solverdata
from misc import color
from time import sleep

def prelimTitle():
    # Print a blank line to create a better visual
    
    print()
    
    # Check to see if we should omit this printoff altogether

    if options.omit_init_printoff:
        return
    
    # Print a new section
    
    print(f"{color.uline}PRELIMINARY INFO:{color.end}")
    print()

def prelimInfo():
    # Check to see if we should omit this printoff altogether

    if options.omit_init_printoff:
        return
    
    # Begin to print the preliminary information
    
    print(f"Constants: L1 = {const.L1},                   Time step: {timedata.time_step}                         KSP type: \"{solverdata.ksp_type}\"")
    print(f"           L2 = {const.L2},                   End time: {timedata.end_time}                         PC type:  \"{solverdata.pc_type}\"")
    print(f"           L3 = {const.L3},                   No. time steps: {timedata.end_time/timedata.time_step:0.0f}")
    print(f"            A = {const.A},")
    print(f"            B = {const.B},")
    print(f"            C = {const.C},")
    print(f"      epsilon = {const.ep},")
    print(f"           L0 = {const.L0}")
    print()
    
    # If we are visualizing this in Paraview, print the path to the Paraview file
    
    if options.visualize:
        print(f"Paraview file: {color.blue}{visdata.file_path}{color.end}")
        print()
    
    # Unless we are manufacturing a solution, print the path to the mesh file; otherwise print the information for the unit cube meshes we will cycle through
    
    if options.manufactured:
        print("Manufactured solution")
        print("Mesh: unit cube mesh")
        print()
        print(f"Init mesh node struc: {meshdata.numnodes_init} x {meshdata.numnodes_init} x {meshdata.numnodes_init}")
        print(f"Max mesh node struc:  {meshdata.numnodes_max} x {meshdata.numnodes_max} x {meshdata.numnodes_max}")
        print()
    else:
        print(f"Mesh: {color.blue}{meshdata.file_path}{color.end}")
        print()

def prelimCompTitle():
    # Wait for 1 second, it looks nicer
    
    sleep(1)
    
    # Print a new section
    
    print(f"{color.uline}PRELIMINARY COMPUTATIONS:{color.end}")
    print()

def prelimCompInfo(time_elapsed):
    print(f"Finished preliminary computations in {time_elapsed:0.2f} seconds.")
    print()

def pdeSolveTitle():
    # Wait for 1 second
    
    sleep(1)
    
    # Print a new section
    
    print(f"{color.uline}PDE SOLVE:{color.end}")
    print()

def pdeSolveInfo(**kwargs):
    print()
    for arg in kwargs:
        if arg == 'mesh_numnodes':
            print(f"Mesh node struc: {kwargs['mesh_numnodes']} x {kwargs['mesh_numnodes']} x {kwargs['mesh_numnodes']}")
        elif arg == 'h1_error':
            print(f"H1 error:        {kwargs['h1_error']:0.15f}")
        elif arg == 'l2_error':
            print(f"L2 error:        {kwargs['l2_error']:0.15f}")
        if arg == 'time_elapsed':
            print(f"Time elapsed:    {kwargs['time_elapsed']:0.2f} seconds")

    print()

# END OF CODE