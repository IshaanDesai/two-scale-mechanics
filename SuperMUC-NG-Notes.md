# Known issues and tips when installing environment on SuperMUC-NG

* Compile newest CMake from source (available versions are too old)
* Use gcc14 (Intel Compiler throws errors when compiling some libs)
* Need to install petsc and slepc from source (no compatible versions available)
* Use stack/24.xx.xx
* Load python/3.10.xx-extended (is not listed but exists)
* When using venv, create env using `python -m venv --system-site-packages --copies venv/py310/` site packages should provide a numpy and mpi4py version linked against intel-mpi. --copies is needed (fails to discover system packages otherwise)
* The following python packages can be installed via wheels:
  * `pip install --no-index --find-links=wheels --no-build-isolation --upgrade exceptiongroup joblib packaging pathspec pkgconfig poetry-core psutil "scikit-build-core[pyproject]" scikit-learn threadpoolctl tomli typing_extensions wheel`
  * when downloading wheels, make sure to have a local python install of version 3.10.xx
  * download with: e.g. `pip download --only-binary=:all: --dest wheels "scikit-build-core[pyproject]"`
  * --upgrade is required as some must be updated
  * all other packages must be installed via build from src to link against the correct MPI version
  * some python packages use an older style for pyproject.toml -> just comment out license-files and replace license with: `license = {text="LICENSE-TYPE"}`
  * some python packages have version 0.0.0 when building from src -> dependency detection will fail, remove dependency version
* As early as possible replace the MPI4PY version with 3.1.6
  * See: [DolfinX nanobind issue](https://fenicsproject.discourse.group/t/no-attribute-qualname-due-to-mpi4py-4-0-0/15375)
* Install h5py from source
  * set `HDF5_MPI=ON`
  * set `HDF5_INCLUDEDIR` to the correct path, find it using `module show hdf5/<your_version>`
* If you (for some reason) decide to compile with intel compiler: install cython with:
  * `export SETUPTOOLS_USE_DISTUTILS=1`
  * `CC=icx LINKCC=icx LDSHARED="icx -shared" pip install cython-3.2.4.zip --no-index --no-build-isolation`
* Configure petsc with:
  * `export PETSC_ARCH=arch-linux-c-opt`
  * `export PETSC_DIR=<petsc_src_dir>`
  * `./configure --prefix=<install_dir> --with-python-exec=$HOME/venv/py310/bin/python --with-shared-libraries=1 --with-memalign=64 --with-petsc4py=1 --with-eigen=1 --with-eigen-dir=$EIGEN_BASE --with-fftw=1 --with-fftw-dir=$FFTW_BASE --with-hdf5=1 --with-hdf5-dir=$HDF5_BASE --with-hypre=1 --with-hypre-dir=$HYPRE_BASE --with-mpi=1 --with-mpi4py=1 --with-nanobind=1 --with-openmp=1 --with-debugging=0`
* Install petsc4py with:
  * `PETSC_ARCH` and `PETSC_DIR` must be set as above
  * from the petsc-src dir: `pip install src/binding/petsc4py --no-index --no-build-isolation`
* Configure slepc with:
  * petsc4py must be installed
  * `PETSC_ARCH` and `PETSC_DIR` must be set as above
  * `export SLEPC_DIR=<slepc_src_dir>`
  * `./configure --prefix=<install_dir> --with-slepc4py`
* When using gcc14 and intel-mpi compiling FANS and pyFANS will result in CMAKE errors
  * comment out include(FetchContent)... and FetchContent_Declare(...), insert find_package(pybind11 REQUIRED) 
  * add FFTW3_INCLUDE_DIRS to all targets, not just FANS_FANS
  * and do the following in ./cmake/modules/FindFFTW3.cmake:

Replace:
```CMake
## MPI
if(PARALLEL STREQUAL "MPI")
  if(MPI_C_LIBRARIES)
    set_target_properties(fftw3::${kind}::mpi PROPERTIES
      IMPORTED_LOCATION "${FFTW3_${KIND}_${PARALLEL}_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIR_PARALLEL}"
      IMPORTED_LINK_INTERFACE_LIBRARIES ${MPI_C_LIBRARIES})
  endif()
endif()
```
with:
```CMake
## MPI
if(PARALLEL STREQUAL "MPI")
  set_target_properties(fftw3::${kind}::mpi PROPERTIES
    IMPORTED_LOCATION "${FFTW3_${KIND}_${PARALLEL}_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIR_PARALLEL}"
  )

  if(MPI_C_LIBRARIES AND NOT MPI_C_LIBRARIES STREQUAL "")
    set_property(TARGET fftw3::${kind}::mpi PROPERTY
      INTERFACE_LINK_LIBRARIES "${MPI_C_LIBRARIES}"
    )
  endif()
endif() 
```

     
Packages built from src:

|                         |                   |                  |                           |               |
|-------------------------|-------------------|------------------|---------------------------|---------------|
| ADIOS2-2.11.0           | ffcx-0.10.1.post0 | mpi4py-3.1.6.zip | precice-3.3.1             | slepc-3.24.2  |
| basix-0.10.0.post0      | gmsh-gmsh_4_15_0  |                  | pugixml-1.15              |               |
| cmake-4.2.3             | h5py-3.15.1.zip   | nanobind         |                           | spdlog-1.17.0 |
| cython-3.2.4.zip        | json-3.12.0       |                  | python-bindings-3.3.1.zip | ufl-2025.2.1  |
| dolfinx-0.10.0.post5    |                   | petsc-3.24.4     |                           |               |
| fenicsx-adapter-develop | mpi4py-3.1.6      |                  |                           |               |

__Remember:__ set `LD_LIBRARY_PATH CMAKE_PREFIX_PATH PATH PKG_CONFIG_PATH` when applicable

Loaded Modules:
```bash
module purge
module load stack/24.5.0
module load intel-toolkit
module load boost/1.84.0-intel25-impi
module switch intel/2025.1.1 gcc/14.3.0
module load boost/1.84.0-gcc14-impi
module load eigen/3.4.0-gcc14
module load fftw/3.3.10-gcc14-impi-openmp
module load hdf5/1.14.5-gcc14-impi
module load hypre/2.32.0-gcc14-impi
module load scotch/7.0.4-gcc14-impi-i64
module load parmetis/4.0.3-gcc14-impi-i64-r64
module load autoconf/
module load glib
module load pkg-config/
module load zlib
module load python/3.10.12-extended
export SLURM_EAR_LOAD_MPI_VERSION="intel"
```
