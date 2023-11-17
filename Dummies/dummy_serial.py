# fakeSerial.py --> dummy_serial.py
# author: D. Thiebaut
# A very crude simulator for PySerial assuming it
#   is emulating an Arduino.


# a Serial class emulator
class DumSerial:

    # init():
    def __init__(self, port='COM1', baudrate=19200, timeout=1,
                 bytesize=8, parity='N', stopbits=1, xonxoff=0, rtscts=0,
                 default_read=b"1"):
        """
        The constructor.  Many of the arguments have default values and can be skipped when calling the constructor.

        :param port:
        :param baudrate:
        :param timeout:
        :param bytesize:
        :param parity:
        :param stopbits:
        :param xonxoff:
        :param rtscts:
        :param default_read: This is returned whenever a read is invoked
        """
        self.name = port
        self.port = port
        self.timeout = timeout
        self.parity = parity
        self.baudrate = baudrate
        self.bytesize = bytesize
        self.stopbits = stopbits
        self.xonxoff = xonxoff
        self.rtscts = rtscts
        self._isOpen = True
        self._receivedData = ""
        self._data = b"It was the best of times.\nIt was the worst of times.\n"
        self.in_waiting = 0
        self.default_read = default_read

    def isOpen(self):  # noqa
        """
        :return: returns True if the port to the Arduino is open.  False otherwise
        """
        return self._isOpen

    def open(self):
        """
        opens the port
        :return: None
        """
        self._isOpen = True

    def close(self):
        """
        closes the port
        :return: None
        """
        self._isOpen = False

    def write(self, bytes):  # noqa
        """
        writes a string of characters to the Arduino
        :param bytes: Ignored
        :return: None
        """
        self.in_waiting = 1

    def read(self):
        """
        reads n characters from the fake Arduino.
        :return: self.default_read
        """
        self.in_waiting = 0
        return self.default_read

    def readline(self):
        """
        reads characters from the fake Arduino until a \n is found.
        :return: self.default_read
        """
        self.in_waiting = 0
        return self.default_read

    def __str__(self):
        """
        :return: a string representation of the serial class
        """
        return "Serial<id=0xa81c10, open=%s>( port='%s', baudrate=%d," \
               % (str(self.isOpen), self.port, self.baudrate) \
               + " bytesize=%d, parity='%s', stopbits=%d, xonxoff=%d, rtscts=%d)" \
               % (self.bytesize, self.parity, self.stopbits, self.xonxoff, self.rtscts)
