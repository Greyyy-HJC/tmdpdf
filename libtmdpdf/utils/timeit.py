import logging
import time
from typing import Optional

log = logging.getLogger(__name__)

class TimeIt:
    """
    Examples:
        >>> with TimeIt('msg'):
            ... # do_something
    """
    def __init__(self, 
                 description: str = None, 
                 logger:logging.Logger=log.info):
        """

        Args:
            description (str, optional): _description_. Defaults to None.
            logger (logging.Logger, optional): _description_. Defaults to log.info.
        """
        self.logger = logger
        self.description = description if description is not None else 'timeit'

    def __enter__(self):
        self.start = time.time()
        self.logger(f'[start] {self.description}', stacklevel=2)

    def __exit__(self, exc_type, exc_value, tb):
        self.logger(f'[end] {self.description}: {(time.time()-self.start):.2f}s', stacklevel=2)
