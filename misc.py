# Value checkers

def check_elastic_constants():
    import printoff as pr
    import const
    if 0>=const.L1 or -const.L1>=const.L3 or const.L3>=2*const.L1 or -3/5*const.L1-1/10*const.L3>=const.L2:
        pr.warning('L1, L2, and L3 do not satisfy the proper inequalities')

# Functions

def getValues(dictionary,keys):
    keys = keys.split(', ')
    values = []
    for key in keys:
        values.append(dictionary[key])
    if len(values) == 1:
        return values[0]
    else:
        values = tuple(values)
        return values

def get_range(range_str):
    range_str = str(range_str)  # ensures we are working with a string
    # print(range_str)
    range_str = range_str.replace(' ','') # get rid of spaces
    # print(range_str)
    range_list = range_str.split(',') # split at commas
    # print(range_list)

    crude_list = []

    for item in range_list:
        item = item.split('-')
        item = list(map(int,item))
        if len(item) == 1:
            pass
        elif len(item) == 2:
            item = list(range(item[0],item[1]+1))
        else:
            raise ValueError('Range must be between no more than two integers.')

        crude_list = crude_list + item
        # print(crude_list)

    crude_list.sort()
    # print(crude_list)

    clean_list = []
    [clean_list.append(item) for item in crude_list if item not in clean_list]
    # print(clean_list)

    return(clean_list)

# Classes

class TimeError(Exception):
    """ A custom time error """

class Timer: # adapted from 'https://realpython.com/python-timer/' on 6/24/20
    def __init__(self):
        self._start_time = None
        self.time_elapsed = None

    def start(self):
        from time import time
        if self._start_time is not None:
            raise TimeError(f"Timer is running. Use .stop() to stop it.")

        self.time_elapsed = None
        self._start_time = time()

    def stop(self):
        from time import time
        if self._start_time is None:
            raise TimeError(f"Timer is not running. Use .start() to start it")

        self.time_elapsed = time() - self._start_time
        self._start_time = None

# END OF CODE
