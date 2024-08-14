import numpy as np
from ctypes import cdll, byref, c_double, POINTER

nasmat = cdll.LoadLibrary('./np_emulator.so')

strains = np.arange(0, 6, dtype=np.float64)
strains_ = strains.ctypes.data_as(POINTER(c_double))

stresses = np.arange(0, 6, dtype=np.float64)
stresses_ = stresses.ctypes.data_as(POINTER(c_double))

nasmat.__mod_nasmat_MOD_nasmat3d_solve.argtypes = [POINTER(c_double)]
nasmat.__mod_nasmat_MOD_nasmat3d_solve.restype = c_double

nasmat.__mod_nasmat_MOD_nasmat3d_solve(strains_, stresses_)


print(stresses)
