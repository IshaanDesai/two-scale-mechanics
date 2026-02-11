import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent / "micro-fans-util"))
from config_io import *

if __name__ == '__main__':
    num_mm_ranks = 4
    num_workers = 4

    target_configs = [
        #(USE_ADA, NO_MADA , NO_STATELESS) ,
        #(USE_ADA, USE_MADA, NO_STATELESS) ,
        #(USE_ADA, USE_MADA, USE_STATELESS),
        #(USE_ADA, NO_MADA , USE_STATELESS),
        #(NO_ADA , USE_MADA, NO_STATELESS) ,
        #(NO_ADA , USE_MADA, USE_STATELESS),
        (NO_ADA , NO_MADA , USE_STATELESS)
    ]

    gen_config(
        num_mm_ranks,
        num_workers,
        NO_SLURM,
        MPI_INTEL,
        decomp_dim=1,
        target_configs=target_configs
    )