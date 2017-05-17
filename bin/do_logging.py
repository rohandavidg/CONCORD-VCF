#!/dlmp/sandbox/cgslIS/rohan/Python-2.7.11/python

"""
setting up logging
"""

import logging
import time
import datetime


def main(filename):
    logger = configure_logger(filename)


def configure_logger(filename):
    """
    setting up logging
    """
    logger = logging.getLogger(filename)
    logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(time.strftime(filename+"-%Y%m%d.log"))
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(asctime)s'\t'%(name)s'\t'%(levelname)s'\t'%(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger


if __name__ == "__main__":
    main(filename)
