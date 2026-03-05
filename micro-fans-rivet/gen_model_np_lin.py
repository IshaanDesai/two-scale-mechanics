import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent / "micro-fans-util"))
from gen_lin_np_model import compute_zero_load_tangent

if __name__ == '__main__':
    compute_zero_load_tangent()