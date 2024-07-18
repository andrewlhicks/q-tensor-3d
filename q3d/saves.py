import os
import yaml
from firedrake import VTKFile, Function, CheckpointFile, DumbCheckpoint, FILE_READ, FILE_CREATE

from pathlib import Path

outfile = None
SaveMode = None

save_modes = ('r','o')
save_modes_legacy = ('resume','overwrite','OVERWRITE')

def initialize(save_mode,save_path,remote=False):
    """ Establishes the SaveMode, as well as a SavePath which depends on
    whether current or legacy SaveMode in use, then repairs the save if
    needed. """

    global SaveMode, SavePath

    SaveMode = save_mode
    
    if SaveMode in save_modes:
        SavePath = save_path
        repair_save(SavePath)
        return
    
    if SaveMode in save_modes_legacy:
        save_name = save_path
        SavePath = f'saves-remote/{save_name}' if remote else f'saves/{save_name}'
        repair_save(SavePath)
        return

    if SaveMode is None:
        raise ValueError('Must have a SaveMode other than None')
    
    raise ValueError('Must choose current ("' + '","'.join(save_modes) + '") or legacy ("' + '","'.join(save_modes_legacy) + '") save mode')

def repair_save(save_path):
    import os

    directories = (f'{save_path}/chk',f'{save_path}/energy',f'{save_path}/vis')

    for directory in directories:
        if not os.path.exists(directory):
            os.makedirs(directory)

def load_checkpoint(*names):
    """ Loads a previously saved checkpoint. Returns mesh
    if no function names are given. Otherwise returns
    function names specified. """
    
    path = f'{SavePath}/chk/checkpoint.h5'

    if not os.path.exists(path):
        raise FileNotFoundError

    with CheckpointFile(path,'r') as file:
        mesh = file.load_mesh()

        if len(names) != 0:
            functions = [file.load_function(mesh, name) for name in names]
            return (*functions,)
        
        return mesh

def checkpoint_exists():
    path = Path(f'{SavePath}/chk/checkpoint.h5')

    if path.exists():
        return True
    
    return False

def save_checkpoint(mesh,*functions):
    path = f'{SavePath}/chk/checkpoint.h5'

    with CheckpointFile(path,'w') as file:
        file.save_mesh(mesh)
        for function in functions:
            file.save_function(function)

def load_dumb_checkpoint(vector_space,name='dump'):
    """ Loads a checkpoint and outputs the function with name specified. """

    path = f'{SavePath}/chk/{name}'

    if not os.path.exists(f'{path}.h5'):
        raise FileNotFoundError

    q_dump = Function(vector_space,name=name)

    with DumbCheckpoint(path,mode=FILE_READ) as chk:
        chk.load(q_dump)

    return q_dump

def save_dumb_checkpoint(q_dump,name='dump'):
    """ Saves a checkpoint the input being the function with name specified. """

    path = f'{SavePath}/chk/{name}'

    if not isinstance(q_dump,Function):
        raise TypeError('Must be a Firedrake Function.')

    q_dump.rename(name)

    with DumbCheckpoint(path,mode=FILE_CREATE) as chk:
        chk.store(q_dump)

def load_energies():
    """ This loader for the energies ensures easy backwards compatibility with older YAML
    documents named 'energies.yml' that did not contain the !EnergyList and !TimeList
    constructors. These older documents were very much dependent on the relative import
    structure of the q-tensor-3d module, and thus imported from 'saves' directly as if
    it were a module in the . directory. A good fix to this is to make the replacement
    that I have specified in the code below. Thi is dependent on the name of the package
    for q-tensor-3d being 'q3d', but I don't plan on changing this name. So this should
    get rid of all backwards-compatibility issues. """

    with open(f'{SavePath}/energy/energies.yml') as energies_file:
        # first ensure backwards compatibility
        energies_file = energies_file.read()
        energies_file = energies_file.replace('!!python/object/new:saves','!!python/object/new:q3d.saves')
        # then load the file into yaml
        yaml_load = yaml.load(energies_file, Loader=yaml.Loader)

    times = yaml_load['times']
    energies = yaml_load['energies']

    if len(times) != len(energies):
        raise ValueError(f'Number of times {len(times)} and number of energies {len(energies)} not equal.')

    return TimeList(times), EnergyList(energies)

def save_energies(times,energies):
    if not isinstance(times,TimeList):
        raise TypeError
    if not isinstance(energies,EnergyList):
        raise TypeError

    yaml_dump = {'times':times, 'energies':energies}

    with open(f'{SavePath}/energy/energies.yml','w') as energies_file:
        energies_file.write(yaml.dump(yaml_dump))

def save_pvd(*args, time=None, path=None, mode=None):
    # default kwarg
    path = path if path is not None else 'vis/vis.pvd'
    mode = mode if mode is not None else 'a' if SaveMode in ('r','resume') else 'w'

    # check if outfile already created earlier. If not, create it
    global outfile
    if outfile is None:
        outfile = VTKFile(f'{SavePath}/{path}', mode=mode)
    
    # write to outfile
    outfile.write(*args, time=time)

# Classes for custom data types

class CustomList(list):
    name = 'Custom list'

    # Dunder methods

    def __init__(self,custom_list):
        if not isinstance(custom_list,list):
            raise TypeError('Must be list type.')
        for i in range(len(custom_list)):
            if not isinstance(custom_list[i],int) and not isinstance(custom_list[i],float):
                raise TypeError('List items must be type int or float.')
        return super().__init__(custom_list)
    def __repr__(self):
        return self.__class__.name + ': ' + super().__repr__()
    def __add__(self,other):
        if not isinstance(other,self.__class__):
            raise TypeError(f'Cannot add {self.__class__} to {other.__class__}')
        return self.__class__(super().__add__(other))
    def append(self,other):
        super().append(other)
        return self.__class__(self)

class EnergyList(CustomList):
    name = 'Energy list'

class TimeList(CustomList):
    name = 'Time list'

    # Dunder methods

    def __init__(self,times):
        self._enforce_ordered(times)
        return super().__init__(times)

    # Static methods

    @staticmethod
    def _enforce_ordered(times):
        for i in range(len(times)-1):
                if times[i] >= times[i+1]:
                    raise ValueError('List is not ordered.')

    # Class methods

    @classmethod
    def by_range(cls,t_initial,t_final=None,step=1):
        if t_final is None:
            t_final = t_initial
            t_initial = 0
        times = list(range(t_initial,t_final,step))
        return cls(times)
    @classmethod
    def by_prev(cls,t_prev,num_times=0,step=1):
        return cls.by_range(num_times).stretch(step).shift(t_prev+step)

    # Properties

    # @property
    # def initial(self):
    #     return self[0]
    @property
    def final(self):
        return self[-1]

    # Regular methods

    def shift(self,shift_factor):
        if not isinstance(shift_factor,int) and not isinstance(shift_factor,float):
            raise TypeError('Shift factor must be type int or float.')
        for i in range(len(self)):
            self[i] += shift_factor
        return self
    def stretch(self,stretch_factor):
        if not isinstance(stretch_factor,int) and not isinstance(stretch_factor,float):
            raise TypeError('Stretch factor must be type int or float.')
        for i in range(len(self)):
            self[i] *= stretch_factor
        return self
    def truncate(self,truncation_length):
        if not isinstance(truncation_length,int):
            raise TypeError('Truncation length must be type int.')
        return self.__class__(self[:truncation_length])

def list_constructor(cls):
    def constructor(loader,node):
        value = loader.construct_sequence(node)
        return cls(value)
    return constructor

def list_representer(tag):
    def representer(dumper,data):
        return dumper.represent_sequence(tag,data)
    return representer

constructors = ('EnergyList','TimeList')

for constructor in constructors:
    constructor_name = '!' + constructor
    constructor = eval(constructor)
    yaml.add_constructor(constructor_name,list_constructor(constructor))
    yaml.add_representer(constructor,list_representer(constructor_name))