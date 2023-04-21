import logging
from functools import partial
from logging import getLogger

from rich.console import Console
from rich.logging import RichHandler
from rich.syntax import Syntax

_console = Console()

def setup_colorlogging():
    FORMAT = "%(message)s"

    logging.basicConfig(
        level=logging.INFO, format=FORMAT, datefmt="[%X]", handlers=[RichHandler()], force=True
    )


setup_colorlogging()
log = getLogger(__name__)


def pprint_code(code: str, format: str):
    """pretty print code

    Args:
        code (str): code str
        format (str): code format, for example 'yaml'
    """
    s = Syntax(code, format)
    _console.print(s)

pprint_yaml = partial(pprint_code, format='yaml')
"""pretty print yaml code"""
