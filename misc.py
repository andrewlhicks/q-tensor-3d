class color:
    header = '\033[95m'
    blue = '\033[94m'
    green = '\033[92m'
    warning = '\033[93m'
    fail = '\033[91m'
    end = '\033[0m'
    bold = '\033[1m'
    uline = '\033[4m'

class TimeError(Exception):
    """ A custom time error """

class Timer: # adapted from 'https://realpython.com/python-timer/' on 6/24/20
    def __init__(self):
        self._start_time = None
        self.elapsed = None
    
    def start(self):
        from time import time
        if self._start_time is not None:
            raise TimeError(f"Timer is running. Use .stop() to stop it.")
        
        self.elapsed = None
        self._start_time = time()
    
    def stop(self):
        from time import time
        if self._start_time is None:
            raise TimeError(f"Timer is not running. Use .start() to start it")
                
        self.elapsed = time() - self._start_time
        self._start_time = None

# END OF CODE