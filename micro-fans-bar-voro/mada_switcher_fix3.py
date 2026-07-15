import numpy as np


def switching_function(resolution, location, t, input, prev_output):
    # marks a slice of width 0.15 of the meso bar mesh on the left most end (has DC BC)
    domain_inner = np.array([[-0.01, -0.01, -0.01], [0.15, 1.01, 1.01]])
    # marks a slice of width 0.3 of the meso bar mesh on the left most end (has DC BC)
    # includes the inner region
    domain_middle = np.array([[-0.01, -0.01, -0.01], [0.3, 1.01, 1.01]])
    point = np.array(location)

    def is_inside(domain, p):
        low = domain[0]
        high = domain[1]
        return np.all((p > low) * (p < high))

    in_middle = is_inside(domain_middle, point)
    in_inner = is_inside(domain_inner, point)

    # coarsest resolution
    if resolution == 0:
        if in_middle or in_inner:
            return 1
        return 0
    # medium resolution
    elif resolution == 1:
        if in_inner:
            return 1
        if in_middle:
            return 0
        # -1 allows going back to lower resolution
        # if one way path is desired return 0 here
        return -1
    # full order model
    elif resolution == 2:
        if in_inner:
            return 0
        # -1 allows going back to lower resolution
        # if one way path is desired return 0 here
        return -1
    else:
        return 0
