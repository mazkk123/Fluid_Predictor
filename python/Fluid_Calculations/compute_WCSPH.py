import math as m
import numpy as np
import random as rd
import re

from compute_SPH import SPH

class WCSPH(SPH):

    EXTRA_PARAMS = {
        "alpha":3,
        "beta":2,
        "SSF":1500,
        "lambda_v":0.1,
        "lambda_f":0.25,
        "lambda_vis":0.3
    }

    def __init__(self):
        
        super().__init__()

        self.compute_K()

    def compute_K(self):
        self.K = m.pow(self.EXTRA_PARAMS["SSF"], 2) * self.PARAMETERS["mass_density"]/ (self.EXTRA_PARAMS["alpha"] + self.EXTRA_PARAMS["beta"])

    def update_pressure(self):
        
        if self.particle.mass_density>2*self.PARAMETERS["mass_density"]:
            pressure_val_1 = np.power((self.particle.mass_density - self.PARAMETERS["mass_density"])/ self.PARAMETERS["mass_density"], self.EXTRA_PARAMS["alpha"])
            pressure_val_2 = np.power((self.particle.mass_density - self.PARAMETERS["mass_density"])/ self.PARAMETERS["mass_density"], self.EXTRA_PARAMS["beta"])
            self.pressure = pressure_val_1 - pressure_val_2
        if self.particle.mass_density >= (self.PARAMETERS["mass_density"]) and self.particle.mass_density < 2*self.PARAMETERS["mass_density"]:
            pressure_val_1 = np.power((self.particle.mass_density)/ self.PARAMETERS["mass_density"], self.EXTRA_PARAMS["alpha"])
            pressure_val_2 = np.power((self.particle.mass_density)/ self.PARAMETERS["mass_density"], self.EXTRA_PARAMS["beta"])
            self.pressure = pressure_val_1 - pressure_val_2
        else:
            self.particle.pressure = 0

    def CFL_condition(self, delta_time):
        
        self.velocity_max = np.linalg.norm(self.particle.velocity)
        return delta_time <= self.EXTRA_PARAMS["lambda_v"]*(self.PARAMETERS["cell_spacing"]/self.velocity_max)

    def viscous_CFL_condition(self, delta_time):
        return delta_time <= self.EXTRA_PARAMS["lambda_vis"]*(self.PARAMETERS["cell_spacing"]/(
                0.6*(self.EXTRA_PARAMS["SSF"]*self.PARAMETERS["viscosity"]) + self.EXTRA_PARAMS["SSF"]
        ))

    def force_CFL_condition(self, delta_time):
        return delta_time <= self.EXTRA_PARAMS["lambda_f"]*(np.sqrt(self.PARAMETERS["cell_spacing"] / np.linalg.norm(self.acceleration)))

    def adaptive_time_stepping(self):
        if self.CFL_condition(self.delta_time) and self.viscous_CFL_condition(self.delta_time) and self.force_CFL_condition(self.delta_time):
            pass