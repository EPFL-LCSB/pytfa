# -*- coding: utf-8 -*-
"""
.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

Logging utilities

"""

import logging
import os
import time

LOGFOLDERNAME = 'logs'

def get_bistream_logger(name):
    """
    Sets up a logger that outputs INFO+ messages on stdout and DEBUG+ messages
    in the log file

    :param name: a class __name__ attribute
    :return:
    """
    logger = logging.getLogger(name)

    try:
        os.makedirs(LOGFOLDERNAME)
    except FileExistsError:
        pass

    if not logger.handlers:
        logger.setLevel(logging.DEBUG)

        # create a file handler
        filename = name+get_timestr()+'.log'
        path = os.path.join(LOGFOLDERNAME, filename)
        file_handler = logging.FileHandler(path)
        file_handler.setLevel(logging.DEBUG)

        # create a stream handler
        stream_handler = logging.StreamHandler()
        stream_handler.setLevel(logging.INFO)

        # create a logging format
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        stream_handler.setFormatter(formatter)

        # add the handlers to the logger
        logger.addHandler(file_handler)
        logger.addHandler(stream_handler)

    return logger

def get_timestr():
    timestr = time.strftime("_%Y%m%d_%H%M%S")
    return timestr

