import math as m
import numpy as np
import random as rd
import re

from compute_SPH import SPH
from Particles.particles import Particle

class PCSPH(SPH):

    OTHER_PARAMS = {
        "max_iterations":2,
        "min_iterations":1,
        "predictor_threshold":0.5,
        "k_stiffness":7,
        "lambda_const":1
    }
    def __init__(self,
                 particle: Particle=None,
                 search_method: str=None,
                 hash_table: dict=None,
                 hash_value: int=None,
                 tank_attrs: dict=None,
                 time_stepping:str = "Euler Cromer",
                 delta_time:float = 0.02):
        
        super().__init__(particle=particle,
                         search_method=search_method,
                         hash_table=hash_table,
                         hash_value=hash_value,
                         time_stepping=time_stepping,
                         delta_time=delta_time,
                         tank_attrs=tank_attrs)

        self.particle = particle
        self.particle.pressure_force = np.array([0, 0, 0], dtype="float64")
        self.particle.pressure = 0

    def update_pressure(self):
        k_const = (self.OTHER_PARAMS["k_stiffness"]*self.PARAMETERS["mass_density"]) / self.OTHER_PARAMS["lambda_const"]
        pressure_term = ((self.particle.mass_density / self.PARAMETERS["mass_density"]) - 1)
        self.particle.pressure = k_const*pressure_term
    
    def update_pressure_force(self):
        return super().update_pressure_force()
    
    def calculate_density_error(self, predicted_density:float =None):
        if predicted_density is not None:
            return predicted_density - self.PARAMETERS["mass_density"]
    
    def find_beta_const(self):
        return m.pow(self.delta_time, 2)*m.pow(self.particle.mass, 2)* \
            2 / m.pow(self.PARAMETERS["mass_density"], 2)
    
    def sum_kernel_grad_dot(self):
        kernel_grad_sum = 0
        for id, nbr_particle in enumerate(self.neighbours_list):
            kernel_grad_sum += (
                np.dot(self.cubic_spline_kernel(nbr_position=nbr_particle.initial_pos, kernel_type=1), 
                       self.cubic_spline_kernel(nbr_position=nbr_particle.initial_pos, kernel_type=1))
            )
        return kernel_grad_sum

    def kernel_grad_sum(self):
        kernel_grad_sum = np.array([0, 0, 0], dtype="float64")
        for id, nbr_particle in enumerate(self.neighbours_list):
            kernel_grad_sum += self.cubic_spline_kernel(nbr_position=nbr_particle.initial_pos, kernel_type=1)
        return kernel_grad_sum

    def find_kronecker_delta(self, density_error:float = None):
        
        beta_const = self.find_beta_const()
        kernel_grad_sum = self.sum_kernel_grad_dot()
        kernel_grad_dot_sum = np.dot(-1*self.kernel_grad_sum(), self.kernel_grad_sum())

        """ print("Dot: " , kernel_grad_dot_sum)
        print("Kernel", kernel_grad_sum)
        print("\n\n\n\n") """

        denom = beta_const * (kernel_grad_dot_sum - kernel_grad_sum)

        """ print("Denom is: ", denom) """

        if denom == 0:
            kronecker_delta = 0
        else:
            kronecker_delta = -1 / denom

        return kronecker_delta

    def update(self):
        
        iterations = 0
        while (self.calculate_density_error(self.density_prediction()) > self.OTHER_PARAMS["predictor_threshold"]):

            if (iterations <= self.OTHER_PARAMS["max_iterations"]):
                break

            self.velocity_prediction = self.prediction_update(particle=self.particle)[1]
            self.position_prediction = self.prediction_update(particle=self.particle)[0]

            self.predicted_density = self.density_prediction()
            self.density_error = self.calculate_density_error(self.predicted_density)

            self.pressure_correction = self.find_kronecker_delta()*self.density_error

            self.particle.pressure += self.pressure_correction

            self.update_pressure_force()

            iterations += 1
        
        return super().update()