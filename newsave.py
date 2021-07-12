import os
import sys
import yaml

help_text = 'python newsave.py <name> <settings-file>.yml <constants-file>.yml'

if len(sys.argv[1:]) != 3:
	print(help_text)
	sys.exit()
if not sys.argv[2].endswith('.yml') or not sys.argv[3].endswith('.yml'):
	print(help_text)
	sys.exit()

# First, load the settings and the constants files

settings_name = sys.argv[2]
settings_path = 'settings/' + settings_name

with open(settings_path) as settings_file:
	settings_dict = yaml.load(settings_file, Loader=yaml.Loader)

constants_name = sys.argv[3]
constants_path = 'constants/' + constants_name

with open(constants_path) as constants_file:
	constants_dict = yaml.load(constants_file, Loader=yaml.Loader)

# Next, create the new save and put the settings and constants file therein

save_name = sys.argv[1]
save_path = 'saves/' + save_name

os.makedirs(save_path)
os.makedirs(save_path+'/chk') # stores checkpoints
os.makedirs(save_path+'/energy') # contains the plot of the energy decrease
os.makedirs(save_path+'/vis') # contains paraview stuff

with open(save_path + '/settings.yml','w') as settings_file:
	settings_file.write(yaml.dump(settings_dict))

with open(save_path + '/constants.yml','w') as constants_file:
	constants_file.write(yaml.dump(constants_dict))
