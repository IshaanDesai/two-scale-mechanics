MKL_FOLDER=/sw/pkgs/arc/intel/2022.1.2/mkl/2022.0.2/lib/intel64

export LD_PRELOAD=$MKL_FOLDER/libmkl_def.so.2:$MKL_FOLDER/libmkl_avx2.so.2:$MKL_FOLDER/libmkl_core.so:$MKL_FOLDER/libmkl_intel_lp64.so:$MKL_FOLDER/libmkl_intel_thread.so:/sw/pkgs/arc/intel/2023.2.1/compiler/2023.2.1/linux/compiler/lib/intel64_lin/libiomp5.so

python3 ruc_nasmat.py
