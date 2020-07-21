from sympyplus import UserDefinedFunction
import yaml

settings_file = open('settings/settings.yml')
settings = yaml.load(settings_file, Loader=yaml.FullLoader)

# Set all of the settings dictionaries

const = settings['const']
options = settings['options']
meshdata = settings['meshdata']
paraview = settings['paraview']
solverdata = settings['solverdata']
timedata = settings['timedata']
userfunc = settings['userfunc']

# Add 'dt' entry to 'const' dictionary

const['dt'] = timedata['time_step']

# Add the user-defined functions

userInitialGuess = UserDefinedFunction(userfunc['initialGuess'])
userBoundary = UserDefinedFunction(userfunc['boundary'])

# END OF CODE