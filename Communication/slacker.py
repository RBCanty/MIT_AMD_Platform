""" slacker.py

This module posts messages to a Slack channel.

@author: Ben C
"""

import requests
from mcn_logging_manager import system_log
from mcn_status import USE_SLACK, SLACK_LEVELS, Fault
import textwrap

WEBHOOK = 'https://hooks.slack.com/services/123456798/123456798ab/123456798abcdefghijklmno'


def post_msg_to_slack(message, *,
                      headers: dict = None, additional_json_data: dict = None,
                      post_kwargs: dict = None,
                      raise_exceptions: bool = False):
    """ Posts a message to the amd_alerts Slack channel

    :param message: A string to be posted in Slack
    :param headers: __future__
    :param additional_json_data: __future__
    :param post_kwargs: __future__
    :param raise_exceptions: True - exceptions are raised, False - all exceptions are caught and returned
    :return: (response contents, response status code [str]) or (exception object, None)
    """
    if headers is None:
        headers = {}
    if post_kwargs is None:
        post_kwargs = {}
    try:
        if not isinstance(message, str):
            message = str(message)
        json_data = {'text': message}
        if additional_json_data and isinstance(additional_json_data, dict):
            json_data.update(additional_json_data)
        response = requests.post(WEBHOOK, headers=headers, json=json_data, **post_kwargs)
        return response.content.decode('utf-8'), str(response.status_code)
    except Exception as e:
        system_log.exception("Slack Bot encountered an Exception")
        if raise_exceptions:
            raise e
        return e, None


def post_fault_to_slack(new_fault: Fault, additional_message=None, **kwargs):
    """
    Reports a Fault object to the amd_alerts Slack channel in a formatted manner.

    :param new_fault: A Fault object being reported
    :param additional_message: (str) Any additional details to post to Slack
    :keyword headers: __future__
    :keyword additional_json_data: __future__
    :keyword post_kwargs: __future__
    :keyword raise_exceptions: True - exceptions are raised, False - all exceptions are caught and returned
    :return: (response contents, response status code [str]), (exception object, None), or ('config-level message', '')
    """
    # Check inputs
    if not isinstance(additional_message, (str, type(None))):
        additional_message = str(additional_message)
    if not USE_SLACK:
        return "Slack Disabled", ""
    if new_fault.level not in SLACK_LEVELS:
        return "Fault not of sufficient level", ""

    # Define constants
    newline = "\n"
    wrap_kwargs = {'width': 35, 'initial_indent': "", 'subsequent_indent': "|  "}

    # Format message
    message = f"```" \
              f"{new_fault.timestamp}\n" \
              f"{new_fault.level} Fault detected!\n" \
              f"| agent: {new_fault.location}\n"
    if new_fault.queue:
        wrapper = textwrap.wrap(f"| queue: {new_fault.queue}", **wrap_kwargs)
        message += f"{newline.join(wrapper)}\n"
    if new_fault.data:
        wrapper = textwrap.wrap(f"| data:  {new_fault.data}", **wrap_kwargs)
        message += f"{newline.join(wrapper)}\n"
    if additional_message:
        message += f"With additional comments:\n    {additional_message}\n"
    message += "```"

    return post_msg_to_slack(message, **kwargs)


if __name__ == '__main__':
    from mcn_status import V_FATAL
    from custom_exceptions import DatabaseRequestError

    USE_SLACK = True

    try:
        raise DatabaseRequestError("The fake queue step could not connect to the database.  This message will be very "
                                   "very long so as to see how Slack handles these really long messages, and if it "
                                   "would be worth tacking on the traceback to the message.  I'm not entirely sure "
                                   "that's a good idea as the 'post_fault_to_slack(Fault)' method in "
                                   "venv/Communication/sacker.py is just gonna get whatever Fault it's passed, so "
                                   "trying to add traceback would need to happen at the level of whoever generates the "
                                   "Fault object from an exception.")
    except DatabaseRequestError as dre:
        my_fault = Fault(location="_U",
                         level=V_FATAL,
                         data=repr(dre),
                         queue="Python has a built-n paragraph formatter... huh")
        resp = post_fault_to_slack(my_fault, additional_message="Hello World")
        print(resp)