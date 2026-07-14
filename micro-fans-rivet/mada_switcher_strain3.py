import numpy as np


def switching_function(resolution, location, t, input, prev_output):
    eps = np.zeros(6)
    eps[0:3] = input["strains1to3"]
    eps[3:6] = input["strains4to6"]
    e1 = 0.5 * (
        np.power(eps[0] - eps[1], 2)
        + np.power(eps[1] - eps[2], 2)
        + np.power(eps[2] - eps[0], 2)
    )
    # we div by 2 to remove sqrt 2 from sig, accounted by factor 3 -> 1.5
    e2 = 1.5 * (np.power(eps[3], 2) + np.power(eps[4], 2) + np.power(eps[5], 2))
    vm_strain = np.sqrt(e1 + e2)

    threshold32 = 5e-5
    threshold64 = 1e-4
    over_threshold32 = vm_strain > threshold32
    over_threshold64 = vm_strain > threshold64

    # coarsest resolution
    if resolution == 0:
        if over_threshold32:
            return 1
        return 0
    # medium resolution
    elif resolution == 1:
        if over_threshold64:
            return 1
        # -1 allows going back to lower resolution
        # if one way path is desired return 0 here
        return 0
    # full order model
    elif resolution == 2:
        if over_threshold64:
            return 0
        # -1 allows going back to lower resolution
        # if one way path is desired return 0 here
        return 0
    else:
        return 0
