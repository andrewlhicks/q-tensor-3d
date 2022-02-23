import functools

# Value checkers

class check:
    def elastic_constants():
        import printoff as pr
        from config import constants as c
        if 0>=c.L1 or -c.L1>=c.L3 or c.L3>=2*c.L1 or -3/5*c.L1-1/10*c.L3>=c.L2:
            pr.warning('L1, L2, and L3 do not satisfy the proper inequalities')
    def energy_decrease(times,energies):
        import printoff as pr
        for i in range(len(energies)-1):
            change_in_energy = energies[i+1]-energies[i]
            if change_in_energy > 0:
                pr.warning(f'Energy decrease failed at time t = {times[i+1]} by {change_in_energy}')

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
        self._time_elapsed = None

    def start(self):
        from time import perf_counter
        if self._start_time is not None:
            raise TimeError(f"Timer is running. Use .stop() to stop it.")

        self._time_elapsed = None
        self._start_time = perf_counter()

    def stop(self):
        from time import perf_counter
        if self._start_time is None:
            raise TimeError(f"Timer is not running. Use .start() to start it")

        self._time_elapsed = perf_counter() - self._start_time
        self._start_time = None

    @property
    def time_elapsed(self):
        return round(self._time_elapsed)

    @property
    def str_seconds(self):
        return f'{self.time_elapsed}s'

    @property
    def str_minutes(self):
        rem = self.time_elapsed

        mns = rem // 60 # floor division
        rem = rem % 60 # remainder

        scs = rem

        return f'{mns}m {scs}s'

    @property
    def str_hours(self):
        rem = self.time_elapsed

        hrs = rem // 3600
        rem = rem % 3600

        mns = rem // 60
        rem = rem % 60

        scs = rem

        return f'{hrs}h {mns}m {scs}s'

    @property
    def str_time(self):
        if self.time_elapsed < 60:
            return self.str_seconds
        if self.time_elapsed < 3600:
            return self.str_minutes
        return self.str_hours

# Decorators

def time_this(func):
    import printoff as pr
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        timer = Timer()
        timer.start()
        value = func(*args, **kwargs)
        timer.stop()
        pr.text(f'Function \'{func.__name__}\' completed in {timer.str_time}',spaced=False)
        return value
    return wrapper

# END OF CODE
