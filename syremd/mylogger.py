from mpi4py import MPI
import logging
import platform


## for now a lot of efforts should be  made to improve this part
## find time to improve it if not busy

_name = 'Logger'

colors = {
    'blue': '\033[94m',
    'yellow': '\033[93m',
    'green': '\033[92m',
    'red': '\033[91m',
    'light_blue': '\033[96m',
    'purple': '\033[95m',
    'end': '\033[0m',
}

log_levels = {
    0: '[CRITICAL]',
    1: '[WARNING]',
    2: '[INFO]',
    3: '[OUT FILES]',
    4: '[DEBUG]'
}

color_prefix = {
    0: colors['red'] + '[CRITICAL]' + colors['end'],
    1: colors['yellow'] + '[WARNING]' + colors['end'],
    2: colors['blue'] + '[INFO]' + colors['end'],
    3: colors['purple'] + '[OUT FILES]' + colors['end'],
    4: colors['green'] + '[DEBUG]' + colors['end'],
}




def myformatter():
    formatter = logging.Formatter(fmt='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(process)d %(message)s',datefmt='%Y-%m-%d,%H:%M:%S')
    comm = MPI.COMM_WORLD
    for r in range(comm.Get_size()):
        if comm.Get_rank() == r:
            handler = logging.FileHandler("/home/shiyi/yank_CHARMM_druce_delicate/test_book/run_log_myremd_{}.log".format(r),"w")
            handler.setFormatter(formatter)
            logger = logging.getLogger("logger")
            logger.setLevel(logging.DEBUG)
            logger.addHandler(handler)
    return logger

