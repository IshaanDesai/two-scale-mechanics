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
        self._sim_id = sim_id

        # File and folder names
        case_filename_initial = 'RUC_initial'
        self._case_filename_restart = 'RUC_iterate'
        foldername = 'ruc_{}'.format(sim_id)

        # Set the working directory to the micro_ruc/ folder
        os.chdir('/home/desaii/composite-multiscale/micro_ruc')

        # Create a new directory for this micro simulation
        os.system('mkdir ruc_{}'.format(sim_id))

        # Copy the input files into the RUC folder
        os.system('cp ' + case_filename_initial + '.inp ' + foldername)
        os.system('cp ' + self._case_filename_restart + '.inp ' + foldername)

        # Copy the job submission file into the RUC folder
        #os.system('cp ' + initial_job_submission_filename + ' ' + foldername)

        # Change the working directory to the ruc_ folder
        os.chdir(foldername)

        # Run the initial Abaqus simulation
        os.system('abaqus job=RUC_{0} input=RUC_initial \
                  scratch=/home/desaii/composite-multiscale/micro_ruc/ruc_{0}  \
                  interactive double=both &> log_ruc_{0}_0.log'.format(sim_id))

    def solve(self, strains, dt):
        print("Micro problem {} is being solved.".format(self._sim_id))
        assert dt != 0

        # Get strain values
        strain_values = np.zeros((6))
        for i in range(3):
            strain_values[i] = strains["strains1to3"][i]
            strain_values[i + 3] = strains["strains4to6"][i]

        inputfile = open(self._case_filename_restart, "w")
        print("Opened the input file")
        
        stresses = np.zeros(6)

        return {"stresses1to3": stresses[0:3], "stresses4to6": stresses[3:6]}
