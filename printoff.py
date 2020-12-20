""" Here, most of the printing is controlled. Most of it is done through
'plogging', i.e. print-logging, where the information is printed to the
console and then put into a log file. """

from settings import const, meshdata, options, visdata, timedata, solverdata
from misc import colors
from time import sleep
import functools

# Decorators

def plogger(func):
    @functools.wraps(func)
    def wrapper_plogger(*args, **kwargs):
        with open(file_path,'a') as file:
            global plog
            def plog(string='',color=None):
                file.write(string+'\n')
                if color is not None:
                    string = colors[color] + string + colors['end']
                print(string)
            value = func(*args, **kwargs)
        return value
    return wrapper_plogger

# Functions that deal with the file

def _set_file_path(new_file_path):
    global file_path
    file_path = new_file_path

def _clear_file():
    from datetime import datetime
    with open(file_path,'w') as file:
        now = datetime.now()
        file.write(now.strftime('%c') + '\n')

# Functions that print lines

def print_lines(*args):
    """ Prints a lines, each line containing a tile and text. Args are
    dictionaries with 'title' and 'text' attributes. """

    def print_line(title,text):
        if not isinstance(title,str):
            raise TypeError('Title must be a string.')

        title_len = len(title)
        spaces_len = int(indent_len - 1 - title_len)
        if spaces_len < 0: spaces_len = 0

        spaces = ''.join([' ' for _ in range(spaces_len)])
        plog(f'{title}:{spaces}{text}')

    indent_len = int(0)

    for arg in args:
        if len(arg['title']) + 2 > indent_len:
            indent_len = len(arg['title']) + 2

    plog()
    for arg in args:
        print_line(arg['title'],arg['text'])
    plog()

# Plogger functions

@plogger
def prelimInfo():
    """ Plogs the preliminary information for the PDE solve. """

    plog()
    plog(f'PRELIMINARY INFO:',color='uline')

    plog()
    plog(f'Constants: L1 = {const.L1},                   Time step: {timedata.time_step}                         KSP type: \'{solverdata.ksp_type}\'')
    plog(f'           L2 = {const.L2},                   End time: {timedata.end_time}                         PC type:  \'{solverdata.pc_type}\'')
    plog(f'           L3 = {const.L3},                   No. time steps: {timedata.end_time/timedata.time_step:0.0f}')
    plog(f'           q0 = {const.q0}')
    plog(f'            A = {const.A},')
    plog(f'            B = {const.B},')
    plog(f'            C = {const.C},')
    plog(f'           W0 = {const.W0},')
    plog(f'           W1 = {const.W1},')
    plog(f'           W2 = {const.W2},')
    plog(f'      epsilon = {const.ep},')
    plog(f'           L0 = {const.L0}')
    plog()

    if options.visualize:
        plog(f'Paraview file: {visdata.file_path}')
        plog()

@plogger
def meshInfo(mesh_name,**kwargs):
    """ Prints the mesh information, namely, the mesh name and any other desired properties. """

    plog(f'MESH INFO:',color='uline')

    # Assemble dicts

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
            dicts.append({'title':'Path','text':f'kwargs[kw]'})

    # Print lines

    print_lines(*dicts)

@plogger
def prelimCompTitle():
    """ Prints title for the preliminary computations. """
    
    sleep(1)
    
    plog(f'PRELIMINARY COMPUTATIONS:',color='uline')
    plog()

@plogger
def prelimCompInfo(time_elapsed):
    """ Prints information after the preliminary computations are finished. """

    plog(f'Finished preliminary computations in {time_elapsed:0.2f} seconds.')
    plog()

@plogger
def pdeSolveTitle():
    """ Prints title for the PDE solving. """
    
    sleep(1)
    
    plog(f'PDE SOLVE:',color='uline')
    plog()

@plogger
def pdeSolveInfo(**kwargs):
    """ Prints desired information for the PDE solving after it is complete. """

    # Assemble dicts

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

@plogger
def warning(text):
    """ Plogs a warning. """

    if not isinstance(text,str):
        raise TypeError('Warnings must be composed of a string.')
    plog(f'Warning: {text}',color='warning')

# END OF CODE