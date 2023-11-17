
"""
# Contains supporting functions for the main functions to call,
# mostly here for convenience really
"""

import os


def directory_check(path_to_check):
    if not os.path.exists(path_to_check):
        os.makedirs(path_to_check)
        return ['Success', 'Made a new directory: %s' % path_to_check]
    else:
        return ['Success', 'Directory already exists: %s' % path_to_check]


def print_progress_bar(iteration, total, prefix='', suffix='', decimals=1, length=100, fill='█', print_end="\r"):
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + '-' * (length - filled_length)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end=print_end)
    if iteration == total:
        print()
