class FromDict:
    def __init__(self,my_dict: dict) -> None:
        for key, value in my_dict.items():
            if isinstance(value,dict):
                setattr(self,key,FromDict(value))
                continue
            setattr(self,key,value)
    def __repr__(self) -> str:
        return repr(self.__dict__)
    def as_dict(self):
        return FromFromDict(self)
class FromFromDict(dict):
    def __init__(self,from_dict: FromDict) -> None:
        stuff = [(attr,getattr(from_dict,attr)) for attr in dir(from_dict) if not attr.startswith('__') and attr != 'as_dict']
        for key, value in stuff:
            if isinstance(value,FromDict):
                self.update({key:FromFromDict(value)})
                continue
            self.update({key:value})

def process_settings(s: FromDict) -> None:
    """ Takes the settings FromDict (s) and provided missing data for outdated formats. """
    try:
        s.mesh.builtin
    except AttributeError:
        s.mesh.builtin = False

def process_constants(c: FromDict) -> None:
    """ Takes the constants FromDict (c) and adds L0, S0, and dt.
    Note that this function may only be called AFTER calling process_settings(). """

    from math import sqrt, ceil

    # set L0 if not explicitly set
    if c.L0 == 'auto':
        c.L0 = ceil(2*(c.A+c.B**2/c.C))

    # set S0 as minimum of double well
    c.S0 = (c.B + sqrt(c.B**2 + 24.0*c.A*c.C))/(4.0*c.C)

    # the following should be deprecated (from settings) in future
    c.dt = settings.time.step

def nondimensionalize(c:FromDict,R:float):
    """ For a constants dictionary const and a length R (for a sphere-like mesh,
    the radius), returns the non-dimensionalized version of these constants """

    from math import sqrt

    # Define auxiliary constants

    AM = max(c.A,c.B,c.C)
    LM = max(c.L1,c.L2,c.L3)

    # Define non-dimensionalized constants

    c.A /= AM
    c.B /= AM
    c.C /= AM
    c.L1 /= LM
    c.L2 /= LM
    c.L3 /= LM
    c.W0 *= (R/LM)
    c.W1 *= (R/LM)
    c.W2 *= (R/LM)
    c.q0 *= R
    c.ep = sqrt(LM/AM)/R
    try:
        c.beta
    except AttributeError:
        c.beta = 1

def initialize(settings_path,constants_path=None):
    from loaddump import load_yml

    # make settings global, thus importable
    global settings

    # load settings as FromDict
    settings = FromDict(load_yml(settings_path))

    # process settings
    process_settings(settings)

    # if no constants file given, exit
    if constants_path is None: return
    
    # make constants global, thus importable
    global constants

    # load constants as FromDict
    constants = FromDict(load_yml(constants_path))

    # nondimensionalize if needed
    if settings.options.nondim:
        nondimensionalize(constants,settings.mesh.radius)

    # process constants
    process_constants(constants)

def main():
    initialize('settings/settings.yml','constants/5cb.yml')
    print(constants)

if __name__ == '__main__':
    main()