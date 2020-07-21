from misc import getValues
from sympyplus import UserDefinedFunction
import yaml

settings_file = open('settings/settings.yml')
settings = yaml.load(settings_file, Loader=yaml.FullLoader)

# Set all of the settings dictionaries

const, options, meshdata, paraview, solverdata, timedata, userfunc = getValues(settings,'const options meshdata paraview solverdata timedata userfunc')

# Add 'dt' entry to 'const' dictionary

const['dt'] = timedata['time_step']

# Add the user-defined functions

userInitialGuess = UserDefinedFunction(userfunc['initialGuess'])
userBoundary = UserDefinedFunction(userfunc['boundary'])

# END OF CODE