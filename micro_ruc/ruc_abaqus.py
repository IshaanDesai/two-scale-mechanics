"""
Python script to replicate the behavior of a RUC. Actually a micro simulation is not solved,
but instead a simple material model is resolved.
The material model is written by Minh Hoang Nguyen.
"""
import numpy as np
import os
from mpi4py import MPI


class MicroSimulation:

    def __init__(self, sim_id):
        """
        Constructor of MicroSimulation class.
        sim_id : int
            global ID of the micro simulation.
        """
        self._dims = 3

        # File and folder names
        case_filename = 'RUC_original'
        job_submission_filename = 'submit_job_micro_RUC.sbat'
        foldername = 'ruc_{}'.format(sim_id)

        # Create a new directory for this micro simulation
        os.system('mkdir ruc_{}'.format(sim_id))

        # Copy the generic input file into the RUC folder
        os.system('cp ' + case_filename + '.inp ' + foldername)

        # Copy the job submission file into the RUC folder
        os.system('cp ' + job_submission_filename + ' ' + foldername)

        # Change the working directory to the ruc_ folder
        os.chdir(foldername)

        # Launch simulation
        os.system('sbatch ' + job_submission_filename)


    def solve(self, strains, dt):
        assert dt != 0



        return {"stresses1to3": stresses[0:3], "stresses4to6": stresses[3:6]}
