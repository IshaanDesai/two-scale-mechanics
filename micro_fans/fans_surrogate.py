"""
Micro simulation Surrogate, require previous computation of surrogate model
"""
import os
import subprocess
from bayesvalidrox import PyLinkForwardModel, Input, PCE, ExpDesigns, Engine
import h5py
import joblib
import numpy as np
import math


class MicroSimulation:

    def __init__(self, sim_id):
        """
        Constructor of MicroSimulation class.
        """
        self._sim_id = sim_id
        self._state = None

        self._model = None
        with open('micro-fans-surrogate.pkl', 'rb') as input:
            self._model = joblib.load(input)
        if self._model is None:
            raise RuntimeError("Failed to load model.")

    def get_state(self):
        return self._state

    def set_state(self, state):
        self._state = state

    def solve(self, macro_data, dt):
        model_input = np.zeros((1, 6))
        model_input[0, 0:3] = macro_data['strains1to3']
        model_input[0, 3:6] = macro_data['strains4to6']
        model_eval, _ = self._model.eval_metamodel(model_input)
        output_data = dict()
        output_data["stresses1to3"] = model_eval["stresses1to3"][0, :]
        output_data["stresses4to6"] = model_eval["stresses4to6"][0, :]
        output_data["cmat1"] = model_eval["cmat1"][0, :]
        output_data["cmat2"] = model_eval["cmat2"][0, :]
        output_data["cmat3"] = model_eval["cmat3"][0, :]
        output_data["cmat4"] = model_eval["cmat4"][0, :]
        output_data["cmat5"] = model_eval["cmat5"][0, :]
        output_data["cmat6"] = model_eval["cmat6"][0, :]
        output_data["cmat7"] = model_eval["cmat7"][0, :]

        return output_data



