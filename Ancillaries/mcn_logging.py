"""
Utility for logging

@author: Ben Canty
"""

import logging
from logging.handlers import TimedRotatingFileHandler
import sys
from custom_classes import pathfinder


LOGGER_LEVELS = ['notset', 'debug', 'info', 'warning', 'error', 'critical']
logging.raiseExceptions = False
# ^ If this doesn't silence the errors, use the handleError override in make_console_handler()


def level_translator(level: str):
    """ Converts logging levels (string) to logging levels (int)

    :param level: A string (or int) representing a logger level
    :return: A logging level int used by logging.log()
    """
    if isinstance(level, int):
        return level
    elif isinstance(level, str):
        level = level.lower()
    else:
        raise ValueError(f"Logger level must be int or string not: {level} ({type(level)})")

    if level == 'critical':
        return 50
    elif level == 'error':
        return 40
    elif level in ['warning', 'warn', 'exception']:
        return 30
    elif level == 'info':
        return 20
    elif level == 'debug':
        return 10
    else:
        return 0


class MyFilter(logging.Filter):
    """
    Blocks logs from filtered sources
    Some imports utilize logging and will hook to this logger, so if they are bothersome they can be blocked by path
    """
    def __init__(self, param=None):  # noqa # Subclassing is not required, only a 'filter()' attribute,
        """ A custom filter for blocking (by path element) other users of logging

        :param param: A list of blocked path elements, e.g. [r'site-packages/urllib3', ]
        """
        self.param = param

    def filter(self, record):
        """ Filters emits

        :param record: A LogRecord object
        :return: If the record was emitted from a blocked source
        """
        if self.param is None:
            allow = True
        else:
            # check the caller's path against all blocked path elements
            allow = not any(blocked in record.pathname for blocked in self.param)
        return allow


# Some modules emit logs to the top-level logger, we do not need them in our logs
muted_speakers = [r'site-packages\urllib3', r'site-packages\matplotlib', r'site-packages\pymcr', ]
f = MyFilter(muted_speakers)

# Define logging formats for the console and logfiles
console_log_format = logging.Formatter("%(asctime)s - %(module)s - %(message)s")
file_log_format = logging.Formatter("%(asctime)s - %(levelname)s - %(module)s:%(funcName)s[%(lineno)s] - %(message)s")

# Retrospective note: deploy a hierarchy of loggers rather than having members invoke the same logger.


def make_logger(level, name=None):
    """ Creates a top-level logger

    :param level: Level at which the logger operates
    :param name: A logger name
    :return: A logging.Logger object
    """
    log = logging.getLogger(name)
    log.setLevel(level)
    return log


def make_console_handler():
    """ Generates a console handler for a Logger using sys.stdout

    :return: logging.StreamHandler
    """
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(console_log_format)
    console_handler.addFilter(f)
    # console_handler.handleError = lambda *_, **__: None
    return console_handler


def make_logfile_path(name):
    """ Converts a log file name into a full path by searching for the file

    :param name: A file name (name.extension)
    :return: A path (path/name.extension)
    """
    return pathfinder(r"Logs", name)


def bind_handler(_logger, file_path):
    """ Binds a file handler to a Logger object

    If the file cannot be accessed, it will try appending "__#" to the file until it can
    or exhausts 5 attempts at which point it will give up

    :param _logger: Logger object
    :param file_path: A full path to a log file
    :return: None
    """
    trial_index = 0
    base_file_path, _ = file_path.split(".log")
    for _ in range(5):
        try:
            file_handler = TimedRotatingFileHandler(file_path, when='midnight', delay=True)
        except FileNotFoundError as fnfe:
            print(f"WARNING: {repr(fnfe)}\nPlease make sure 'venv/Logs/' exists")
            return
        except Exception as e:  # TODO: Replace with correct Exception once known (probably PermissionError)
            print(repr(e))
        else:
            file_handler.setFormatter(file_log_format)
            file_handler.addFilter(f)
            _logger.addHandler(file_handler)
            print(f"Logging file handler loaded with {file_path}")
            return
        file_path = base_file_path + f"__{trial_index}.log"
        trial_index += 1
    print("WARNING: LOGGER COULD NOT BIND A FILE (OR MODIFIED FILE)\n\tLOGGING WILL NOT BE CAPTURED BY A FILE")
