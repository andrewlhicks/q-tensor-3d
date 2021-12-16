class FromDict:
    def __init__(self,my_dict: dict) -> None:
        for key, value in my_dict.items():
            if isinstance(value,dict):
                setattr(self,key,FromDict(value))
                continue
            setattr(self,key,value)
    def __repr__(self) -> str:
        return repr(self.__dict__)

def process_settings(s: FromDict) -> None:
    """ Takes the settings FromDict (s) and provided missing data for outdated formats. """
    try:
        s.mesh.builtin
    except AttributeError:
        s.mesh.builtin = False

def process_constants(c: FromDict) -> None:
    """ Takes the constants FromDict (c) and adds L0, S0, and dt """
    from math import sqrt, ceil

    # set L0 if not explicitly set
    if c.L0 == 'auto':
        c.L0 = ceil(2*(c.A+c.B**2/c.C))

    # set S0 as minimum of double well
    c.S0 = (c.B + sqrt(c.B**2 + 24.0*c.A*c.C))/(4.0*c.C)

    # the following should be deprecated (from settings) in future
    c.dt = settings.time.step

def _config(settings_path,constants_path):
    from loaddump import load_yml

    # make settings, constants global, thus importable
    global settings
    global constants

    # load settings, constants as FromDict
    settings = FromDict(load_yml(settings_path))
    constants = FromDict(load_yml(constants_path))

    # process settings, costants
    process_settings(settings)
    process_constants(constants)

def main():
    _config('settings/settings.yml','constants/5cb.yml')
    print(constants)

if __name__ == '__main__':
    main()