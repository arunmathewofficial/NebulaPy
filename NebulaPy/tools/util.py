import sys

def nebula_exit(errorMessage):
    """
    Custom exit function display error before exiting.
    :param message: Optional exit message.
    """
    RED = "\033[91m"
    RESET = "\033[0m"
    print(f' error: {RED}{errorMessage}{RESET}')
    sys.exit()