# Running on SuperMUC-NG

Running scripts "normally" using sbatch is impossible with the implemented task model of the micro-manager.
The provided sbatch configuration does not allow for spawning of sub-processes.
Therefore, the used workaround is by repurposing the sbatch scripts to execute using salloc.

To provide a smooth user experience, please use the provided run.sh script in the jobs directory. Usage:
```bash
bash run.sh <job-name>
```
job-name can be any file name within jobs except run or utils (without the .sh ending)

Example
```bash
bash run.sh calc-bar-sphere
```

The run script extracts all needed information from the sbatch preamble of each job-script, then allocated resources
using the found information with salloc. 
As the final step it disowns the allocation. By doing so, it seems like a regular sbatch job from the outside.

Every job copies and modifies its required config filed into a job local working environment, such that multiple jobs can work in parallel.

Prior to running any job, make sure the required configs and executables are available.

- Micro Manager requires 
  - at least PyFANS.so in each micro-fans folder
  - Model adaptive runs require PyFANS0 to PyFANS2 (compiled with input0.json to input2.json as default arg to input file loading)
  - calc and scale scripts use stateless version of configs
  - sphere meshes need to be generated, use the gen_spheres.py script for that
  - configs are generated using the gen_configs.py script