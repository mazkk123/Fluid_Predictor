import math as m
import numpy as np
import random as rd
import cupy as cp
import re
import time

from Particles.particles import Particle
from Fluid_Calculations.compute_sph import SPH

class WCSPH(SPH):

    EXTRA_PARAMS = {
        "alpha":3,
        "beta":2,
        "SSF":1500,
        "lambda_v":0.1,
        "lambda_f":0.25,
        "lambda_vis":0.3
    }

    def __init__(self,
                 particle: Particle=None,
                 search_method: str=None,
                 hash_table:dict=None,
                 hash_value:int=None,
                 time_stepping:str = "Euler Cromer",
                 time_schemes:dict = None,
                 collision_types:dict = None,
                 params: dict = None,
                 all_particles:list = None,
                 tank_attrs:dict = None,
                 delta_time:float = None):
        
        super().__init__(particle=particle,
                        all_particles=all_particles, 
                        time_stepping=time_stepping,
                        search_method=search_method,
                        collision_types=collision_types,
                        params=params,
                        time_schemes=time_schemes,
                        hash_table=hash_table,
                        hash_value=hash_value,
                        tank_attrs=tank_attrs,
                        delta_time=delta_time)

        self.compute_K()

    def compute_K(self):
        self.K = m.pow(self.EXTRA_PARAMS["SSF"], 2) * self.PARAMETERS["mass_density"] / (self.EXTRA_PARAMS["alpha"] + self.EXTRA_PARAMS["beta"])

    def update_pressure(self):
        
        if self.particle.mass_density>2*self.PARAMETERS["mass_density"]:
            pressure_val_1 = np.power((self.particle.mass_density - self.PARAMETERS["mass_density"])/ self.PARAMETERS["mass_density"], self.EXTRA_PARAMS["alpha"])
            pressure_val_2 = np.power((self.particle.mass_density - self.PARAMETERS["mass_density"])/ self.PARAMETERS["mass_density"], self.EXTRA_PARAMS["beta"])
            self.pressure = self.K*(pressure_val_1 - pressure_val_2)
        if self.particle.mass_density >= (self.PARAMETERS["mass_density"]) and self.particle.mass_density < 2*self.PARAMETERS["mass_density"]:
            pressure_val_1 = np.power((self.particle.mass_density)/ self.PARAMETERS["mass_density"], self.EXTRA_PARAMS["alpha"])
            pressure_val_2 = np.power((self.particle.mass_density)/ self.PARAMETERS["mass_density"], self.EXTRA_PARAMS["beta"])
            self.pressure = self.K*(pressure_val_1 - pressure_val_2)
        else:
            self.particle.pressure = 0

    def update(self):
        return super().update()