import numpy as np


def switching_function(resolution, location, t, input, prev_output):
    domain_128 = np.array([[-0.01, 0.46, -0.025], [0.31, 0.54, 0.075]])
    domain_64 = np.array([[-0.01, 0.375, -0.025], [0.31, 0.625, 0.075]])
    point = np.array(location)

    def is_inside(domain, p):
        low = domain[0]
        high = domain[1]
        return np.all((p > low) * (p < high))

    in_d64 = is_inside(domain_64, point)
    in_d128 = is_inside(domain_128, point)

    if resolution == 0:
        if in_d64 or in_d128: return 1
        return 0
    elif resolution == 1:
        if in_d128: return 1
        if in_d64: return 0
        return -1
    elif resolution == 2:
        if in_d128: return 0
        return -1
    else:
        return 0

