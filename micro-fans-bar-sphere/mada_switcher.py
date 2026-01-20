import numpy as np


def switching_function(resolution, location, t, input, prev_output):
    if resolution < 1: return 1
    else: return 0
