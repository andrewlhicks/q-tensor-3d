import json
import yaml

def load_json(path):
    with open(path) as file:
        dict = json.loads(file.read())
    return dict

def dump_json(dump,path):
    with open(path,'w') as file:
        file.write(json.dumps(dump,indent=2))

def load_yml(path):
    with open(path) as file:
        dict = yaml.load(file, Loader=yaml.Loader)
    return dict

def dump_yml(dump,path):
    with open(path,'w') as file:
        file.write(yaml.dump(dump))
