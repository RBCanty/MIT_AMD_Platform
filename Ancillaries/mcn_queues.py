""" mcn_queues
Extension of the standard Queue library for logging and try/catch

@author: Ben Canty
"""

import queue
import message
from collections import deque
from mcn_logging_manager import system_log


class MCNPriorityQueue(queue.PriorityQueue):
    """ Extends :class:`queue.PriorityQueue`

    Puts queue.get() and queue.put() methods into try/except statements.
    Adds logging functionality, contains, and a __dict__ method for json
    """

    def __init__(self):
        """ Creates an MCNPriorityQueue(queue.PriorityQueue)
        """
        super(MCNPriorityQueue, self).__init__()

    def get_next(self, nowait=False):
        """ Wrapper on queue.get() and queue.get_nowait()

        :param nowait: If get (False) or get_nowait (True) should be used
        :return: The next element on the queue or None if empty or an error occurred
        """
        try:
            if nowait:
                return self.get_nowait()
            else:
                return self.get()
        except queue.Empty:
            return None
        except:  # noqa
            system_log.exception(f"Unexpected exception when pulling next message")
            return None

    def enqueue(self, msg):
        """ Wrapper on queue.put()

        :param msg: Enqueued object (in implementation: a Message object)
        :return: True/False - If the object was successfully enqueued
        """
        if not isinstance(msg, message.Message):
            system_log.info("Tried to enqueue a non-Message object")
            return False
        if (not msg[message.TYPE]) and (not msg[message.SENDER]):
            return False
        try:
            self.put(msg)
            return True
        except:  # noqa
            system_log.exception(f"Failed to enqueue: {msg}")
            return False

    def __dict__(self):
        return {"MCNPriorityQueue": list(self.queue)}

    def contains(self, item, _eq=None):
        """ Helper method for determining if the queue contains an item

        :param item: For what is being searched
        :param _eq: Default (None) will use 'item in queue', otherwise use 'any(_eq(item, i) for i in queue)'.  The
          function must take two items and return a boolean value in accordance with their equivalence.
        :return: True if present in queue, False otherwise
        """
        if _eq is None:
            return item in self.queue
        else:
            return any(_eq(item, q_item) for q_item in self.queue)


class LogQueue(queue.Queue):
    """ Mutates :class:`queue.Queue`

    Modifies __init__ and _init methods to enforce a maximum size that replaces existing elements rather than blocking
    enqueueing when full.
    Adds __contains__ functionality and put_nodup() method to control the addition of new elements
    """
    def __init__(self, size):
        """ Creates an LogQueue(queue.Queue)
        """
        super(LogQueue, self).__init__(maxsize=size)
        self.maxsize = 0

    def _init(self, maxsize):
        self.queue = deque(maxlen=maxsize)

    def __contains__(self, item):
        return item in self.queue

    def put_nodup(self, item):
        """
        Adds an item to the queue such that there are no duplicates.
        Duplicate entries are removed (most recent survives)

        :param item: new element for queue
        :return:
        """
        if item in self.queue:
            self.queue.remove(item)
        super(LogQueue, self).put(item)
