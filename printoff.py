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

def prelimInfo():
    # Print a blank line to create a better visual
    
    print()

    # Print a new section
    
    print(f"{color.uline}PRELIMINARY INFO:{color.end}")
    print()
    
    # Begin to print the preliminary information
    
    print(f"Constants: L1 = {const.L1},                   Time step: {timedata.time_step}                         KSP type: \"{solverdata.ksp_type}\"")
    print(f"           L2 = {const.L2},                   End time: {timedata.end_time}                         PC type:  \"{solverdata.pc_type}\"")
    print(f"           L3 = {const.L3},                   No. time steps: {timedata.end_time/timedata.time_step:0.0f}")
    print(f"           q0 = {const.q0}")
    print(f"            A = {const.A},")
    print(f"            B = {const.B},")
    print(f"            C = {const.C},")
    print(f"           W0 = {const.W0},")
    print(f"           W1 = {const.W1},")
    print(f"           W2 = {const.W2},")
    print(f"      epsilon = {const.ep},")
    print(f"           L0 = {const.L0}")
    print()
    
    # If we are visualizing this in Paraview, print the path to the Paraview file
    
    if options.visualize:
        print(f"Paraview file: {color.blue}{visdata.file_path}{color.end}")
        print()

def meshInfo(mesh_name,**kwargs):
    print(f"{color.uline}MESH INFO:{color.end}")

    # assemble dicts

    dicts = []

    dicts.append({'title':'Mesh','text':mesh_name})

    for kw in kwargs:
        if kw == 'numnodes_init':
            dicts.append({'title':'Init mesh node struc','text':f'{kwargs[kw]} x {kwargs[kw]} x {kwargs[kw]}'})
        elif kw == 'numnodes_final':
            dicts.append({'title':'Final mesh node struc','text':f'{kwargs[kw]} x {kwargs[kw]} x {kwargs[kw]}'})
        elif kw == 'no_refinements':
            dicts.append({'title':'No. refinements','text':kwargs[kw]})
        elif kw == 'file_path':
            dicts.append({'title':'Path','text':f'{color.blue}kwargs[kw]{color.end}'})

    # for kw in kwargs:
    #     if kw == 'numnodes_init':
    #         print(f"Init mesh node struc: {kwargs['numnodes_init']} x {kwargs['numnodes_init']} x {kwargs['numnodes_init']}")
    #     elif kw == 'numnodes_max':
    #         print(f"Init mesh node struc: {kwargs['numnodes_max']} x {kwargs['numnodes_max']} x {kwargs['numnodes_max']}")
    #     elif kw == 'no_refinements':
    #         print(f"No. refinements: {kwargs['no_refinements']}")
    #     elif kw == 'file_path':
    #         print(f"Path: {color.blue}{kwargs['file_path']}{color.end}")
    # print()

    print_lines(*dicts)

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

def print_lines(*args):
    """ Args are dictionaries with 'title' and 'text' """

    def print_line(title,text):
        if not isinstance(title,str):
            raise TypeError('Title must be a string.')

        title_len = len(title)
        spaces_len = int(indent_len - 1 - title_len)
        if spaces_len < 0: spaces_len = 0

        spaces = ''.join([' ' for _ in range(spaces_len)])
        print(f"{title}:{spaces}{text}")

    indent_len = int(0)

    for arg in args:
        if len(arg['title']) + 2 > indent_len:
            indent_len = len(arg['title']) + 2

    print()
    for arg in args:
        print_line(arg['title'],arg['text'])
    print()
    
def pdeSolveInfo(**kwargs):
    # assemble dicts

    dicts = []

    for kw in kwargs:
        if kw == 'mesh_numnodes':
            dicts.append({'title':'Mesh node struc','text':f'{kwargs[kw]} x {kwargs[kw]} x {kwargs[kw]}'})
        elif kw == 'h1_error':
            dicts.append({'title':'H1 error','text':f'{kwargs[kw]:0.15f}'})
        elif kw == 'l2_error':
            dicts.append({'title':'L2 error','text':f'{kwargs[kw]:0.15f}'})
        elif kw == 'time_elapsed':
            dicts.append({'title':'Time elapsed','text':f'{kwargs[kw]:0.2f} seconds'})
        elif kw == 'refinement_level':
            dicts.append({'title':'Refinement level','text':kwargs[kw]})
        elif kw == 'energy':
            dicts.append({'title':'Energy','text':kwargs[kw]})
        elif kw == 'custom':
            dicts.append(kwargs[kw])

    print_lines(*dicts)

def warning(text):
    if not isinstance(text,str):
        raise TypeError('Warnings must be composed of a string.')
    print(f'{color.warning}Warning: {text}{color.end}')

# END OF CODE