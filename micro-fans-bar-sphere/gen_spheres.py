import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent / "micro-fans-util"))
from gen_sphere_meshes import generate_spheres

generate_spheres()