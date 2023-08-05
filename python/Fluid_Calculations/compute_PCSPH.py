import math as m
import numpy as np
import random as rd
import re

from compute_SPH import SPH
from Particles.particles import Particle

class PCSPH(SPH):

    OTHER_PARAMS = {
        "max_iterations":3,
        "min_iterations":1,
        "predictor_threshold":0.5
    }
    def __init__(self,
                 particle: Particle=None,
                 search_method: str=None,
                 hash_table: dict=None,
                 hash_value: int=None,
                 params: dict=None,):
        
        super().__init__(particle=particle,
                         search_method=search_method,
                         hash_table=hash_table,
                         hash_value=hash_value,
                         params=params)

        self.particle = particle
        self.particle.pressure_force = np.array([0, 0, 0])
        self.particle.pressure = 0

    def update_pressure_force(self):

        pressure_force = np.array([0, 0, 0])
        for id, nbr_particle in enumerate(self.neighbours_list):
            pressure_avg = self.particle.pressure/m.pow(self.particle.mass_density,2) + \
                            nbr_particle.pressure / m.pow(nbr_particle.mass_density, 2)
            kernel_grad = self.kernel_gradient(
                self.particle.initial_pos - nbr_particle.initial_pos
            )
            pressure_force += kernel_grad*pressure_avg
        self.particle.pressure_force = -m.pow(self.PARAMETERS["mass"], 2)*pressure_force
    
    def calculate_density_error(self, predicted_density:float = None):
        return predicted_density - self.PARAMETERS["mass_density"]
    
    def find_beta_const(self):
        return m.pow(self.delta_time, 2)*m.pow(self.particle.mass, 2)* \
            2 / m.pow(self.PARAMETERS["mass_density"], 2)
    
    def sum_kernel_grad_dot(self):
        kernel_grad_sum = 0
        for id, nbr_particle in enumerate(self.neighbours_list()):
            kernel_grad_sum += (
                np.dot(self.kernel_gradient(nbr_particle.initial_pos, 1), 
                       self.kernel_gradient(nbr_particle.initial_pos, 1))
            )
        return kernel_grad_sum

    def kernel_grad_sum(self):
        kernel_grad_sum = 0
        for id, nbr_particle in enumerate(self.neighbours_list()):
            kernel_grad_sum += self.kernel_gradient(nbr_particle.initial_pos, 1)
        return kernel_grad_sum

    def find_kronecker_delta(self, density_error:float = None):
        
        beta_const = self.find_beta_const()
        kernel_grad_sum = self.sum_kernel_grad_dot()
        kernel_grad_dot_sum = np.dot(-1*self.kernel_grad_sum(), self.kernel_grad_sum())

        return -1 / beta_const * (kernel_grad_dot_sum - kernel_grad_sum)

    def update(self):

        iterations = 0
        while self.calculate_density_error() > self.OTHER_PARAMS["predictor_threshold"] or \
            iterations < self.OTHER_PARAMS["max_iterations"]:

            self.velocity_prediction = self.prediction_update(particle=self.particle)[1]
            self.position_prediction = self.prediction_update(particle=self.particle)[0]

            self.predicted_density = self.density_prediction()
            self.density_error = self.calculate_density_error(self.predicted_density)

            self.pressure_correction = self.find_kronecker_delta()*self.density_error
            self.pressure += self.pressure_correction

            self.update_pressure_force()

            iterations += 1
        
        return super().update()