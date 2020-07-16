# This python program is meant to print off any relevant information to the user in a visually pleasing way. For the initial printoff, we first print a blank line to create visual
# separation between the regular command line and the information printed off here. Then the program checks to see if the initial printoff should be omitted or not. From this point
# forward, if a blank line is needed for visual separation, it will be printed after the previous print and not before.

def initPrintoff():
    from settings import A, B, C, dt, end, ep, ksp_type, L1, L2, L3, manufactured, mesh_numnodes_init, mesh_numnodes_max, mesh_filepath, omit_init_printoff, paraview_filepath, pc_type, visualize
    from misc import color
    from time import sleep
    
    # Print a blank line to create a better visual
    
    print()
    
    # Check to see if we should omit this printoff altogether
    
    if omit_init_printoff:
        return
    
    # Print a new section
    
    print(f"{color.uline}PRELIMINARY INFO:{color.end}")
    print()
    
    # Begin to print the preliminary information
    
    print(f"Constants: L1 = {L1},                   Time step: {dt}                         KSP type: \"{ksp_type}\"")
    print(f"           L2 = {L2},                   End time: {end}                         PC type:  \"{pc_type}\"")
    print(f"           L3 = {L3},                   No. time steps: {end/dt:0.0f}")
    print(f"            A = {A},")
    print(f"            B = {B},")
    print(f"            C = {C}")
    print(f"      epsilon = {ep}")
    print()
    
    # If we are visualizing this in Paraview, print the path to the Paraview file
    
    if visualize:
        print(f"Paraview file: {paraview_filepath}")
        print()
    
    # Unless we are manufacturing a solution, print the path to the mesh file; otherwise print the information for the unit cube meshes we will cycle through
    
    if manufactured:
        print("Manufactured solution")
        print("Mesh: unit cube mesh")
        print()
        print(f"Init mesh node struc: {mesh_numnodes_init} x {mesh_numnodes_init} x {mesh_numnodes_init}")
        print(f"Max mesh node struc:  {mesh_numnodes_max} x {mesh_numnodes_max} x {mesh_numnodes_max}")
        print()
    else:
        print(f"Mesh: {mesh_filepath}")
        print()

def calctimePrintoff(time_elapsed):
    from misc import color
    from time import sleep
    
    # Wait for 1 second, it looks nicer
    
    sleep(1)
    
    # Print a new section
    
    print(f"{color.uline}PRELIMINARY CALCULATIONS:{color.end}")
    print()
    print(f"Finished preliminary calculations in {time_elapsed:0.2f} seconds.")
    print()
    
    # Wait for 1 second
    
    sleep(1)
    
    # Print a new section
    
    print(f"{color.uline}PDE SOLVE:{color.end}")
    print()

def summaryPrintoff(time_elapsed):
    print(f"Finished PDE solve in {time_elapsed:0.2f} seconds.")
    print()

def summaryPrintoffManufactured(mesh_numnodes,H1_error,L2_error,time_elapsed):
    print(f"Mesh node struc: {mesh_numnodes} x {mesh_numnodes} x {mesh_numnodes}")
    print(f"H1 error:        {H1_error:0.15f}")
    print(f"L2 error:        {L2_error:0.15f}")
    print(f"Time elapsed:    {time_elapsed:0.2f} seconds")
    print()

# END OF CODE