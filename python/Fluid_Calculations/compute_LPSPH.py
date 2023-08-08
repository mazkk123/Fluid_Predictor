import math as m
import numpy as np
import random as rd
import re

from Particles.particles import Particle
from compute_SPH import SPH

class LPSPH(SPH):

    OTHER_PARAMS = {
        "min_threshold":0.2,
        "max_threshold":0.5
    }

    def __init__(self, ref_radius:float=None,
                 p: Particle=None,
                 search_method: str=None,
                 hash_table:dict=None,
                 hash_value:int=None,
                 delta_time:float = None):
        
        super().__init__(particle=p,
                        search_method=search_method,
                        hash_table=hash_table,
                        hash_value=hash_value,
                        delta_time=delta_time)

        if ref_radius is not None:
            self.ref_radius = ref_radius

        self.epsilon = 2/3*self.ref_radius

        self.near_particles = []
        self.far_particles = []

        self.pressure_far = np.array([0, 0, 0])
        self.pressure_near = np.array([0, 0, 0])

        self.pressure_neighbourhood()

        self.density_error = self.particle.mass_density - self.PARAMETERS["mass_density"]

    def find_h_i(self):
        pass
    
    def pressure_neighbourhood(self):
        for id, nbr_particle in enumerate(self.neighbours_list):
            if self.is_near_particle(nbr_particle.initial_pos):
                self.near_particles.append(nbr_particle)
            elif self.is_far_particle(nbr_particle.initial_pos):
                self.far_particles.append(nbr_particle)
            else:
                self.far_particles[id] = 0
                self.near_particles[id] = 0

    def is_near_particle(self, position: np.array=None):
        if position is not None:
            return self.epsilon <= np.linalg.norm(self.particle.initial_pos - position ) and \
                np.linalg.norm(self.particle.initial_pos - position) <= self.find_h_i()

    def is_far_particle(self, position: np.array=None):
        if position is not None:
            return self.epsilon > np.linalg.norm(self.particle.initial_pos - position ) and \
                np.linalg.norm(self.particle.initial_pos - position) <= self.find_h_i()
        
    def near_pressure(self):
        for id, far_nbr in enumerate(self.near_particles):
            density_diff = (far_nbr.mass_density - self.PARAMETERS["mass_density"]) * far_nbr.mass_density * \
                            m.pow(self.ref_radius, 2)
            div_const = 2 *self.PARAMETERS["mass_density"] * m.pow(self.delta_time, 2)
            self.pressure_far += density_diff/div_const

    def far_pressure(self):
        
        for id, far_nbr in enumerate(self.far_particles):
            density_diff = far_nbr.mass_density - self.PARAMETERS["mass_density"]
            div_const = 4 * np.pi * self.PARAMETERS["mass_density"] *  \
                        m.pow(self.delta_time, 2) * np.linalg.norm(self.particle.initial_pos - 
                                                                   far_nbr.initial_pos)
            self.pressure_far += density_diff/div_const

    def calculate_density_error(self, predicted_density:float =None):

        if predicted_density is not None:
            return predicted_density - self.PARAMETERS["mass_density"]
    
    def pressure(self):
        self.particle.pressure = self.near_pressure() + self.far_pressure()

    def update_mass_density(self):
        density, density_divider = 0 
        for id, nbr_particle in enumerate(self.neighbours_list):
            density_num = self.kernel_linear(self.particle.initial_pos - nbr_particle.initial_pos, 0)*nbr_particle.mass
            density_div = nbr_particle.mass/ nbr_particle.mass_density * \
                            self.kernel_linear(self.particle.initial_pos - nbr_particle.initial_pos, 0)
            density += density_num
            density_divider += density_div

        self.particle.mass_density = self.PARAMETERS["mass_density"] + density/density_divider
    
    def update_time_step(self):
        
        self.particle.velocity += self.delta_time * self.particle.pressure_force/ self.particle.mass
        self.particle.initial_pos += m.pow(self.delta_time, 2) * self.particle.pressure_force/ self.particle.mass

    def update_all_forces(self):
        self.update_mass_density()
        self.update_gravity()
        self.update_buoyancy()
        self.update_surface_tension()

        self.all_forces = self.particle.mass_density + \
                        self.gravity + self.buoyancy + \
                        self.particle.surface_tension
        
    def local_domain(self):
        pass
    
    def update(self):

        while self.calculate_density_error(self.density_prediction()) > self.OTHER_PARAMS["max_threshold"]:

            self.velocity_prediction = self.prediction_update(particle=self.particle)[1]
            self.position_prediction = self.prediction_update(particle=self.particle)[0]

            self.pressure()

            self.predicted_density = self.density_prediction()
            self.density_error = self.calculate_density_error(self.predicted_density)

            self.update_pressure_force()
            self.update_time_step()
        
        return super().update()

