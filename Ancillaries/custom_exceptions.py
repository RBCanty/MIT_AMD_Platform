"""
Additional Exceptions

@author: Ben Canty
"""

# A list of Exceptions to subclass
"""
Exception
 +-- StopIteration
 +-- StopAsyncIteration
 +-- ArithmeticError
 |    +-- FloatingPointError
 |    +-- OverflowError
 |    +-- ZeroDivisionError
 +-- AssertionError
 +-- AttributeError
 +-- BufferError
 +-- EOFError
 +-- ImportError
 |    +-- ModuleNotFoundError
 +-- LookupError
 |    +-- IndexError
 |    +-- KeyError
 +-- MemoryError
 +-- NameError
 |    +-- UnboundLocalError
 +-- OSError
 |    +-- BlockingIOError
 |    +-- ChildProcessError
 |    +-- ConnectionError
 |    |    +-- BrokenPipeError
 |    |    +-- ConnectionAbortedError
 |    |    +-- ConnectionRefusedError
 |    |    +-- ConnectionResetError
 |    +-- FileExistsError
 |    +-- FileNotFoundError
 |    +-- InterruptedError
 |    +-- IsADirectoryError
 |    +-- NotADirectoryError
 |    +-- PermissionError
 |    +-- ProcessLookupError
 |    +-- TimeoutError
 +-- ReferenceError
 +-- RuntimeError
 |    +-- NotImplementedError
 |    +-- RecursionError
 +-- SystemError
 +-- TypeError
 +-- ValueError
 +-- Warning
 |    +-- RuntimeWarning
 |    +-- SyntaxWarning
 |    +-- UserWarning
 |    +-- FutureWarning
 |    +-- ImportWarning
 |    +-- UnicodeWarning
 |    +-- BytesWarning
 |    +-- ResourceWarning
"""


class ExampleCustomException(Exception):
    """ Example of how to make custom exceptions

    Most of the time the Exception can literally just be a pass statement
    """
    def __init__(self, *args, **kwargs):
        super(ExampleCustomException, self).__init__(args, kwargs)
        self.custom_field: str = kwargs.get('name_of_field', 'default value')

    def get_custom_field(self) -> str:
        """ Optional (since you can access fields via the dot operator anyway) method for accessing a custom field
        using a getter.

        Nice to use if you want type hints

        :return: The value of your field
        """
        return self.custom_field


# Database
class DatabaseGeneralException(Exception):
    """ Parent class for Database-based exceptions"""
    pass


class BadQueueFormatError(DatabaseGeneralException):
    """ Exception for when inconsistencies are noticed in the queue document pulled from the database

    BadQueueFormatError.init(duration=int) allows this exception to store the last known duration (default: -1)

    BadQueueFormatError.get_duration() returns the last known duration
    """
    def __init__(self, *args, **kwargs):
        duration = kwargs.pop('duration', -1)
        super(BadQueueFormatError, self).__init__(args, kwargs)
        self.last_known_duration = duration

    def get_duration(self):
        return self.last_known_duration


class DatabaseRequestError(DatabaseGeneralException):
    """ Exception for when a call to the database returns with an error code """
    pass


# Scheduler
class TaskExecutionError(RuntimeError):
    """ When the Task issued was not completed """
    pass


class CoordinatorUserStop(RuntimeError):
    """ To signal that the user has issued a stop """
    pass


class CoordinatorSchedulingError(Exception):
    """ To signal the scheduler has failed (not used) """
    pass


class NoViableOptionsFound(Exception):
    """ To signal no operations can be scheduled without creating a conflict (or there are none to be scheduled) """
    pass
    # def __init__(self, *args, **kwargs):
    #     super(NoViableOptionsFound, self).__init__(args, kwargs)
    #     print("No Viable Options Found")


# Serial Controllers
class ControllerError(Exception):
    """ When a controller encounters an error which is not related to the (serial) com port

     Such as a missing com/controller or a safety check being failed"""
    pass


# Server
class MemberListUpdated(Exception):
    """ To signal that the members of the server have changed """
    pass


# Confirmations
class ConfirmationTimeout(ResourceWarning):
    """ Signals when a timeout has occurred when waiting for an operation confirmation """
    def __init__(self, *args, **kwargs):
        super(ConfirmationTimeout, self).__init__(args, kwargs)
        self.task_key = kwargs.get('task_key', None)

    def get_task_key(self):
        return self.task_key


class UserVoidCheckpoint(KeyError):
    """ Raised when the user removed a checkpoint and so the reference becomes undefined """
    def __init__(self, *args, **kwargs):
        super(UserVoidCheckpoint, self).__init__(args, kwargs)
        self.checkpoint_key = kwargs.get('checkpoint', None)

    def get_checkpoint_key(self):
        return self.checkpoint_key


class StatusAccessTimeout(RuntimeWarning):
    """ Raised when access to the status object times out """
    pass


# Stata
class FaultValueError(ValueError):
    """ When the reading/writing of a Fault object fails """
    def __init__(self, *args, **__):
        super(FaultValueError, self).__init__(*args)


# Loops
class BreakOuterLoop(StopIteration):
    """ Helper to control nested loops """
    def __init__(self, *args, **kwargs):
        super(BreakOuterLoop, self).__init__(args, kwargs)
        self.memory = kwargs.get('memory', None)

    def get_memory(self):
        return self.memory


# Spark
class TrayOverload(RuntimeError):
    """ When the Spark tray would be double-occupied """
    def __init__(self, *args, **__):
        super(TrayOverload, self).__init__("Plate being loaded onto occupied space!", *args)


class MissingPlate(RuntimeWarning):
    """ When the Spark tray would double-unoccupied """
    def __init__(self, *args, **__):
        super(MissingPlate, self).__init__("Plate being unloaded from an empty space!", *args)


if __name__ == '__main__':
    """
    print("Example 1____")
    example = ExampleCustomException(name_of_field='custom value')
    print(example.get_custom_field())

    print("Example 1____")
    try:
        raise ExampleCustomException(name_of_field='example value')
    except ExampleCustomException:
        pass

    print("Example 2____")
    try:
        raise ExampleCustomException(name_of_field='example value')
    except ExampleCustomException as e:
        e.print_traceback()
        # Note how it doesn't always print in order

    print("Example 3____")
    try:
        raise ExampleCustomException(name_of_field='example value')
    except ExampleCustomException as e:
        print(f"'{e.get_custom_field()}' (using getter) or '{e.custom_field}' (using dot)")

    print("Example 4____")
    try:
        raise ExampleCustomException(name_of_field='example value')
    except ExampleCustomException as e:
        print(e.__traceback__)
    """

    try:
        raise ControllerError
    except DatabaseGeneralException:
        print("Caught")
