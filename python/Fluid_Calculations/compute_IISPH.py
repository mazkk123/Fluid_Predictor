import math as m
import numpy as np
import random as rd
import re

from compute_SPH import SPH
from Particles.particles import Particle

class IISPH(SPH):

    OTHER_PARAMS = {
        "lambda_stiffness":7,
        "k_stiffness":7,
        "max_threshold":0.3,
        "max_iterations":2,
        "relaxation_factor":0.5
    }

    def __init__(self, 
                 p: Particle=None,
                 search_method: str=None,
                 hash_table:dict=None,
                 hash_value:int=None,
                 params:dict = None,
                 delta_time:float = None):
        
        super().__init__(particle=p,
                        search_method=search_method,
                        hash_table=hash_table,
                        hash_value=hash_value,
                        params=params,
                        delta_time=delta_time)
        
        self.predict_advection()

    def predict_advection(self):
        self.mass_density_adv = 0
        self.predict_displacement()
        self.predict_velocity_advection()

        self.update_mass_density_advection()
        self.particle.iter_pressure = self.prev_pressure * 0.5
        self.update_acceleration_advection()

    def pressure_solve(self):
        iteration = 0
        while self.calculate_density_error(self.density_prediction()) < self.OTHER_PARAMS["max_threshold"] \
                and iteration<self.OTHER_PARAMS["max_iterations"]:
            
            self.update_displacement_iter()
            self.update_iter_pressure()
            self.particle.pressure = self.particle.iter_pressure

            iteration +=1 
        
        self.update_pressure_force()
            
    def update_mass_density(self):
        return super().update_mass_density()
    
    def predict_displacement(self):
        displacement = np.array([0, 0, 0])
        for id, nbr_particle in enumerate(self.neighbours_list):
            density_div = -1* nbr_particle.mass/ m.pow(self.particle.mass_density, 2)
            kernel_gradient = self.kernel_gradient(self.particle.initial_pos - nbr_particle.initial_pos)
            displacement += density_div*kernel_gradient

        self.particle.displacement = displacement*m.pow(self.delta_time, 2)

    def update_mass_density_advection(self):

        advection_amt = np.array([0, 0, 0])
        for id, nbr_particle in enumerate(self.neighbours_list):
            mass_advection = nbr_particle.mass * (self.particle.advected_velocity - self.particle.velocity)
            kernel_gradient = self.kernel_gradient(self.particle.initial_pos - nbr_particle.initial_pos, 1)
            advection_amt += mass_advection*kernel_gradient
        self.mass_density_adv = advection_amt*self.delta_time + self.particle.mass_density

    def predict_velocity_advection(self):

        self.update_gravity()
        self.update_buoyancy()
        self.update_viscosity()
        self.update_surface_tension()

        self.particle.advected_velocity = self.particle.velocity + self.delta_time* \
                                        (self.gravity + self.buoyancy + self.viscosity + 
                                         self.surface_tension)/self.mass

    def update_acceleration_advection(self):
        acc_adv = np.array([0, 0 , 0])
        for id, nbr_particle in enumerate(self.neighbours_list):
            acc_adv += (
                nbr_particle.mass * (self.particle.displacement - nbr_particle.displacement ) *
                self.kernel_gradient(self.particle.initial_pos - nbr_particle.initial_pos, 0)
            )
        self.particle.acceleration_adv = acc_adv
    
    def calculate_density_error(self, predicted_density:float =None):

        if predicted_density is not None:
            return predicted_density - self.PARAMETERS["mass_density"]
    
    def update_displacement_iter(self):
        displacement = np.array([0, 0, 0])
        for id, nbr_particle in enumerate(self.neighbours_list):
            density_div = -1* nbr_particle.mass/ m.pow(self.particle.mass_density, 2)
            kernel_gradient = self.kernel_gradient(self.particle.initial_pos - nbr_particle.initial_pos)
            displacement += density_div*self.particle.iter_pressure*kernel_gradient

        self.particle.displacement = displacement*m.pow(self.delta_time, 2)

    def update_pressure(self):
        pressure_const = self.OTHER_PARAMS["k_stiffness"] * self.PARAMETERS["mass_density"]
        pressure_calc = m.pow((self.particle.mass_density / self.PARAMETERS["mass_density"]), 
                             self.OTHER_PARAMS["lambda_stiffness"]) - 1
        self.particle.pressure = pressure_calc*pressure_const
        self.particle.prev_pressure = self.particle.pressure

    def neighbour_displacement(self):
        nbr_displacement = np.array([0, 0, 0])
        
        for nbr_particle in self.neighbours_list:
            nbr_displacement += nbr_particle.displacement_iter*nbr_particle.iter_pressure

        return nbr_displacement

    def displacement_dif(self):
        disp_dif = np.array([0, 0, 0])

        for nbr_particle in self.neighbours_list:
            disp_dif += (self.particle.displacement_iter*nbr_particle.iter_pressure) - \
                    (nbr_particle.displacement_iter*nbr_particle.iter_pressure)

        return disp_dif

    def update_iter_pressure(self):
        relaxation_const = (1 - self.OTHER_PARAMS["relaxation_factor"])*self.particle.iter_pressure
        acceleration_adv_const = self.OTHER_PARAMS["relaxation_factor"] * 1/self.particle.acceleration_adv
        disp_final = np.array([0, 0, 0])
        for id, nbr_particle in enumerate(self.neighbours_list):
            neighbour_displacement = self.neighbour_displacement()
            displacement_difference = self.displacement_dif()
            disp_final += (displacement_difference - neighbour_displacement) * \
                        self.kernel_gradient(self.particle.initial_pos - nbr_particle.initial_pos) * \
                        nbr_particle.mass
            
        self.particle.iter_pressure = (self.mass_density_adv - self.PARAMETERS["mass_density"] - disp_final) \
                                    *acceleration_adv_const + relaxation_const
    
    def update(self):
        
        self.pressure_solve()
        
        self.particle.velocity = self.particle.advected_velocity + self.delta_time*self.particle.pressure_force/self.particle.mass
        self.particle.initial_pos += self.particle.velocity*self.delta_time
        self.XSPH_vel_correction()

        self.choose_collision_types("Cuboid", secondary_type="Normal", particle=self.particle)