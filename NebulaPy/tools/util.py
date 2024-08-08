import sys
from NebulaPy import version

def nebula_exit_with_error(errorMessage):
    """
    Custom exit function display error before exiting.
    :param message: Optional exit message.
    """
    RED = "\033[91m"
    RESET = "\033[0m"
    print(f' error: {RED}{errorMessage}{RESET}')
    print(f' NebulaPy {version.__version__} exiting ...')
    sys.exit()

def nebula_warning(warnMessage):
    """
    Custom exit function display error before exiting.
    :param message: Optional exit message.
    """
    RED = "\033[91m"
    RESET = "\033[0m"
    print(f' error: {RED}{warnMessage}{RESET}')
    sys.exit()