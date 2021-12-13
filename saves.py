import os
from config import settings
import functools

outfile = None
SaveMode = None
SaveName = None

def initialize(save_mode,save_name,remote=False):
	global SaveMode, SaveName, current_directory

	SaveMode = save_mode
	SaveName = save_name
	current_directory = f'saves/{SaveName}'
	if remote == True:
		current_directory = f'saves-remote/{SaveName}'

def load_checkpoint(vector_space,name='dump'):
	""" Loads a checkpoint and outputs the function with name specified. """
	from firedrake import Function, DumbCheckpoint, FILE_READ

	q_dump = Function(vector_space,name=name)

	with DumbCheckpoint(f'{current_directory}/chk/{name}',mode=FILE_READ) as chk:
		chk.load(q_dump)

	return q_dump

def save_checkpoint(q_dump,name='dump'):
	""" Saves a checkpoint the input being the function with name specified. """
	from firedrake import Function, DumbCheckpoint, FILE_CREATE

	if not isinstance(q_dump,Function):
		raise TypeError('Must be a Firedrake Function.')

	q_dump.rename(name)

	with DumbCheckpoint(f'{current_directory}/chk/{name}',mode=FILE_CREATE) as chk:
	    chk.store(q_dump)

def load_energies():
	import yaml

	with open(f'{current_directory}/energy/energies.yml') as energies_file:
		yaml_load = yaml.load(energies_file, Loader=yaml.Loader)

	times = yaml_load['times']
	energies = yaml_load['energies']

	return TimeList(times), EnergyList(energies)

def save_energies(times,energies):
	import yaml
	if not isinstance(times,TimeList):
		raise TypeError
	if not isinstance(energies,EnergyList):
		raise TypeError

	yaml_dump = {'times':times, 'energies':energies}

	with open(f'{current_directory}/energy/energies.yml','w') as energies_file:
		energies_file.write(yaml.dump(yaml_dump))

def save_pvd(*args,time=None):
	from firedrake import File
	global outfile
	if outfile is None:
		mode = 'a' if SaveMode == 'resume' else 'w'
		outfile = File(f'{current_directory}/vis/vis.pvd',mode=mode)
	outfile.write(*args,time=time)

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
