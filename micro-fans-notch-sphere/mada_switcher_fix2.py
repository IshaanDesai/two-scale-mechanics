import numpy as np


def switching_function(resolution, location, t, input, prev_output):
    domain_inner = np.array([[-0.01, 0.46, -0.025], [0.31, 0.54, 0.075]])
    point = np.array(location)

    def is_inside(domain, p):
        low = domain[0]
        high = domain[1]
        return np.all((p > low) * (p < high))

    in_inner = is_inside(domain_inner, point)

    if resolution == 0:
        if in_inner: return 1
        return 0
    elif resolution == 1:
        if in_inner: return 0
        return -1
    else:
        return 0

