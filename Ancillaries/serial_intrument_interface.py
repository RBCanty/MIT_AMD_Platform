"""
Tools for generic serial communication with robotic devices.

SerialInstrument (class):
  * Deploys a Lock object to enforce only one caller to the instrument at a time
  * If no com port is given OR if the is_sim argument is set True, the class property is_sim will be True
  * Deploys a disconnect method

serial_protector (decorator factory):
  * Used to create decorators for subclasses of SerialInstrument that utilize the Lock and is_sim properties
  * It can be used to create a decorator for each method individually or a general decorator can be created after the
  __init__ definition (if present)
  * When simulated, it will print a message and return None
  * When real, it will execute the specified function then sleep (handled exceptions are returned, otherwise they are
  raised)

Example::

| class MyDevice(SerialInstrument):
|   def __init__(com, pos_arg1, is_sim, baudrate, **kwargs):
|     super(MyDevice, self).__init__(com, is_sim, baudrate, **kwargs)
|     ...device-specific initializations...
|
|   serial_protect = serial_protector((IOError, KeyError), 5, 0.1)
|
|   @serial_protect
|   def my_method(args, kwargs):
|     ...stuff...
|
|   @serial_protector((IOError, KeyError), 10, 3)  # Long mechanical operation
|   def long_operation(args, kwargs):
|     ...stuff...

@author: Ben Canty
"""

import time
from functools import wraps
from threading import Lock
from typing import Tuple, Type

import serial


class SerialInstrument:
    """
    A superclass for instruments that require serial communication.

    Comes with an initialization that will handle opening or simulating a serial connection and will have all required
    properties to support the serial_protector decorator.  Also includes a disconnect method

    Intended for use with the serial_protector decorator
    """
    def __init__(self, com, is_sim=False, baudrate=9600, **kwargs):
        """
        Creates a :class:`SerialInstrument` instance

        References :class:`serial.Serial`

        :param com: (str) Name of the com port (passed to serial.Serial)
        :param is_sim: (bool) True/False - If the device should be simulated (None returns)
        :param baudrate: (int) The baud rate (passed to serial.Serial)
        :param kwargs: Any arguments matching those of serial.Serial's constructor are passed into serial.Serial

        :keyword bytesize: (EIGHTBITS) Number of data bits
        :keyword parity: (PARITY_NONE) Enable parity checking
        :keyword stopbits: (STOPBITS_ONE) Number of stop bits
        :keyword timeout: (None) Set a read timeout value in seconds
        :keyword xonxoff: (False) Enable software flow control
        :keyword rtscts: (False) Enable hardware (RTS/CTS) flow control
        :keyword write_timeout: (None) Set a write timeout value in seconds
        :keyword dsrdtr: (False) Enable hardware (DSR/DTR) flow control
        :keyword inter_byte_timeout: (None) Inter-character timeout
        :keyword exclusive: (None) Set exclusive access mode (POSIX only)
        """
        self._lock = Lock()
        self.is_sim = False
        if (com is None) or is_sim:
            self.is_sim = True

        if self.is_sim:
            print(f"Simulated {type(self).__name__}:: serial opened on {com}")
            self.ser = com
        else:
            self.ser = serial.Serial(com, baudrate, **_serial_args_protector(**kwargs))

    def sub_init(self, *_, **__):
        """ Allows sub-/composed-classes to invoke their own local init method to handle their own needs """
        pass

    def disconnect_controller(self):
        """ Disconnects (releases lock and closes serial com)

        :return: None
        """
        if self.is_sim:
            print(f"Simulated {type(self).__name__}:: serial closed")
        else:
            if self._lock.locked():
                self._lock.release()
            self.ser.close()

    def close(self):
        """ Calls disconnect_controller()

        Mimics serial's close() method nomenclature

        :return:
        """
        self.disconnect_controller()

    def write(self, cmd, ignore_exceptions=()):
        """ Mimics serial's write() method, but with simulation protection.

        :param cmd: passed to serial.write()
        :param ignore_exceptions: Return Exceptions rather than raising (can be empty tuple)
        :return: None (sim), return from self.ser.write(cmd), or exception object
        """
        if self.is_sim:
            return print(f'Simulated {type(self).__name__}:: serial write "{cmd}"')
        try:
            return self.ser.write(cmd)
        except ignore_exceptions as ignored:
            return ignored

    def serial_call(self, attribute, args=None, kwargs=None, call_in_place=True, ignore_exceptions=()):
        """ Redirects call onto the serial object, but with simulation protection.

        :param attribute: Name of the attribute being called
        :param args: args for the attribute (tuple)
        :param kwargs: kwargs for the attribute (dictionary)
        :param call_in_place: True - calls the attribute with the given args ans kwargs, False - returns attribute
        :param ignore_exceptions: (Only applies to the call on the attribute), ignored exceptions are returned
        rather than raised (can be empty tuple)
        :return: (None if sim), the attribute, the return from the attribute, or an ignored exception
        :raises AttributeError:  from getattr(serial, attribute)
        """
        if kwargs is None:
            kwargs = {}
        if args is None:
            args = []
        if self.is_sim:
            return print(f'Simulated {type(self).__name__}:: serial "{attribute}(args={args}, kwargs={kwargs})"')
        ser_attr = getattr(self.ser, attribute)
        if call_in_place:
            try:
                return ser_attr(*args, **kwargs)
            except ignore_exceptions as ignored:
                return ignored
        else:
            return ser_attr


def serial_protector(handled_exceptions: Tuple[Type[Exception], ...], *,
                     timeout=-1, post_sleep=1, return_timeout=False,
                     alternate_method=None):
    """ Decorator Factory for Serial objects.

    When simulated, all calls return None; when real, calls return the wrapped function's return or an Exception object
    for all handled exceptions.

    The TimeoutError from Lock.acquire() can be returned or raised (see return_timeout).

    :param handled_exceptions: A tuple of all exceptions that are considered "handled" (can be an empty tuple)
    :param timeout: Number of seconds to wait for a Lock.acquire() call (-1 for unbounded wait)
    :param post_sleep: Number of seconds to wait after the call to the Serial object
    :param return_timeout: (Default False) If true, a timeout on the Lock acquisition returns the exception object
    otherwise the exception is raised
    :param alternate_method: Defines custom response to a function call from simulated mode, the method is given self,
    func.__name__, *args, and **kwargs (in that order).  The default behavior is to print a message that the call was
    simulated and then return None.
    :return: None on success, An exception object on
    """
    def decorator(func):
        @wraps(func)
        def wrapper_function(self, *args, **kwargs):
            if self._lock.acquire(timeout=timeout):
                try:
                    time.sleep(0.001)
                    if self.is_sim:
                        if alternate_method:
                            return alternate_method(self, func.__name__, *args, **kwargs)
                        print(f"{type(self).__name__}: Simulated Serial Call on {func.__name__}")
                        return None
                    else:
                        return func(self, *args, **kwargs)
                except handled_exceptions as r_exc:
                    return r_exc
                finally:
                    time.sleep(post_sleep)
                    self._lock.release()
                    time.sleep(0.001)
            else:
                timeout_exception = TimeoutError("Timeout when fetching MCN Status Object")
                if return_timeout:
                    return timeout_exception
                else:
                    raise timeout_exception
        return wrapper_function
    return decorator


def _serial_args_protector(**kwargs):
    """ Isolated serial arguments from the generic kwargs given to the class constructor

    References: :class:`SerialInstrument`

    :param kwargs: A set of keyword arguments given to a SerialInstrument constructor
    :return: Only keywords expected by a serial.Serial constructor
    """
    ser_kwargs = dict()
    ser_kwargs['bytesize'] = kwargs.get('bytesize', serial.EIGHTBITS)
    ser_kwargs['parity'] = kwargs.get('parity', serial.PARITY_NONE)
    ser_kwargs['stopbits'] = kwargs.get('stopbits', serial.STOPBITS_ONE)
    ser_kwargs['timeout'] = kwargs.get('timeout', None)
    ser_kwargs['xonxoff'] = kwargs.get('xonxoff', False)
    ser_kwargs['rtscts'] = kwargs.get('rtscts', False)
    ser_kwargs['write_timeout'] = kwargs.get('write_timeout', None)
    ser_kwargs['dsrdtr'] = kwargs.get('dsrdtr', False)
    ser_kwargs['inter_byte_timeout'] = kwargs.get('inter_byte_timeout', None)
    ser_kwargs['exclusive'] = kwargs.get('exclusive', None)
    return ser_kwargs
