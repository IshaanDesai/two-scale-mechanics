import numpy as np

def switching_function(resolution, location, t, input, prev_output):
    if prev_output is None:
        return 0

    sig = np.zeros(6)
    sig[0:3] = prev_output["stresses1to3"]
    sig[3:6] = prev_output["stresses4to6"]
    s1 = 0.5 * (
            np.power(sig[0] - sig[1], 2)
            + np.power(sig[1] - sig[2], 2)
            + np.power(sig[2] - sig[0], 2)
    )
    # we div by 2 to remove sqrt 2 from sig, accounted by factor 3 -> 1.5
    s2 = 1.5 * (
            np.power(sig[3], 2) + np.power(sig[4], 2) + np.power(sig[5], 2)
    )
    vm_stress = np.sqrt(s1 + s2)
    threshold = 0.015
    over_threshold = vm_stress > threshold

    # ROM
    if resolution == 0:
        if over_threshold: return 1
        return 0
    # FOM
    elif resolution == 1:
        if over_threshold: return 0
        # -1 allows going back to lower resolution
        # if one way path is desired return 0 here
        return 0
    else:
        return 0

