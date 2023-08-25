import math as m
import numpy as np
import random as rd
import re
import time

from Fluid_Calculations.compute_sph import SPH
from Particles.particles import Particle

class PCSPH(SPH):

    OTHER_PARAMS = {
        "max_iterations":3,
        "min_iterations":1,
    }

    def __init__(self,
                 particle: Particle=None,
                 search_method: str=None,
                 hash_table:dict=None,
                 hash_value:int=None,
                 time_stepping:str = "Euler Cromer",
                 all_particles:list = None,
                 time_schemes:dict = None,
                 params:dict = None,
                 collision_types:dict = None,
                 tank_attrs:dict = None,
                 delta_time:float = None,
                 temperature:bool = False):
        
        super().__init__(particle=particle,
                        all_particles=all_particles, 
                        time_stepping=time_stepping,
                        search_method=search_method,
                        hash_table=hash_table,
                        hash_value=hash_value,
                        time_schemes=time_schemes,
                        collision_types=collision_types,
                        params=params,
                        tank_attrs=tank_attrs,
                        temperature=temperature,
                        delta_time=delta_time)
        
        self.iterations = 0
        
        self.particle.pressure_force = np.array([0, 0, 0], dtype="float64")
        self.particle.pressure = 0

        self.update_predicted_attrs(self.particle, 4)

        self.particle.predicted_density = self.particle.mass_density
    
    def find_beta_const(self):
        return m.pow(self.delta_time, 2)*m.pow(self.particle.mass, 2)* \
            2 / m.pow(self.PARAMETERS["mass_density"], 2)

    # ---------------------------------------------------------------------- PREDICT DENSITY -----------------------------------------------------------------------------

    def update_del_x(self, particle):
        particle.delta_x = 0
        pressure_const = 2*particle.pressure_correction / m.pow(self.PARAMETERS["mass_density"], 2)
        constant = -1*m.pow(self.delta_time, 2)*particle.mass*pressure_const

        for nbr_particle in particle.neighbour_list:
            try:
                particle.delta_x += (
                    self.cubic_spline_kernel_gradient(particle.predicted_initial_pos - nbr_particle.initial_pos)
                )
            except TypeError:
                particle.delta_x += 0
        particle.delta_x *= constant

    def update_density_change(self, particle):
        sum_del_neighbour_x, sum_weights = 0, 0
        particle.density_change = 0
        for nbr_particle in particle.neighbour_list:
            try:
                sum_weights += (
                    self.cubic_spline_kernel_gradient(particle.predicted_initial_pos - 
                                                    nbr_particle.predicted_initial_pos)
                )
            except TypeError:
                sum_weights += 0
            try:
                sum_del_neighbour_x += (
                    self.cubic_spline_kernel_gradient(particle.predicted_initial_pos - 
                                                    nbr_particle.predicted_initial_pos)*nbr_particle.delta_x
                )
            except TypeError:
                sum_del_neighbour_x += 0
        particle.density_change = particle.mass*(particle.delta_x*sum_weights - sum_del_neighbour_x)

    def update_next_density(self):
        
        for nbr in self.neighbours_list:
            self.update_del_x(nbr)
            self.update_density_change(nbr)
            
        self.update_del_x(self.particle)
        self.update_density_change(self.particle)
        
        self.particle.predicted_density += self.particle.density_change

    def calculate_density_error(self, particle):
        return particle.predicted_density - self.PARAMETERS["mass_density"]

    # -------------------------------------------------------------------- COMPUTE PRESSURE -----------------------------------------------------------------------------

    def find_kronecker_delta(self, particle):
        beta_value = self.find_beta_const()
        accum_denom, denom = 0, 0
        kronecker = 0
        for nbr_particle in particle.neighbour_list:

            try:
                kernel_value = self.cubic_spline_kernel_gradient(particle.predicted_initial_pos - nbr_particle.predicted_initial_pos)
            except ValueError:
                kernel_value = 0

            try:
                accum_denom += m.pow(kernel_value, 2)
                denom += kernel_value
            except TypeError:
                accum_denom += 0
                denom += 0

        denominator = (
                    -denom*denom - accum_denom
                )
        
        if denominator==0:
            kronecker = 0
        else:
            kronecker = -1 / (
                    beta_value * denominator
            )

        return kronecker

    def update_pressure_correction(self):
        for nbr in self.neighbours_list:
           nbr.pressure_correction += self.find_kronecker_delta(nbr)*self.calculate_density_error(nbr)
            
        self.particle.pressure_correction += self.find_kronecker_delta(self.particle)*self.calculate_density_error(self.particle)
    
    def update_pressure(self):
        return super().update_pressure()
    
    def update_pressure_force(self):
        return super().update_pressure_force()

    # -------------------------------------------------------------------- UPDATE CALLS ------------------------------------------------------------------------------

    def update_predicted_attrs(self, particle, depth:int=3):
        
        if depth==0:
            return 
        
        self.find_neighbour_list(particle)
        self.update_predicted_mass_density(particle)
        self.update_predicted_gravity(particle)
        self.update_predicted_surface_tension(particle)
        self.update_predicted_viscosity(particle)
        self.update_predicted_buoyancy(particle)

        self.all_forces = particle.gravity + \
                            particle.surface_tension + \
                            particle.viscosity + \
                            particle.buoyancy

        particle.predicted_velocity = particle.velocity + self.delta_time*self.all_forces / particle.mass
        particle.predicted_initial_pos = particle.initial_pos + particle.predicted_velocity*self.delta_time

        for nbr in particle.neighbour_list:
            return self.update_predicted_attrs(nbr, depth-1)

    def update_all_forces(self):
        
        self.update_pressure()
        self.update_pressure_force()
        self.update_viscosity()
        self.update_buoyancy()
        self.update_surface_tension()
        self.update_gravity()

        self.all_forces = self.particle.pressure_force + \
                          self.particle.viscosity + \
                          self.particle.gravity + \
                          self.particle.surface_tension + \
                          self.particle.buoyancy

    def update(self):
        
        while (self.calculate_density_error(self.particle) > 0.01*self.PARAMETERS["mass_density"]) and \
            self.iterations<self.OTHER_PARAMS["max_iterations"]:

            self.update_pressure_correction()
            self.update_next_density()
            
            self.iterations += 1
        
        self.particle.mass_density = self.particle.predicted_density

        self.update_all_forces()
        self.normal_field()
        self.XSPH_vel_correction()

        self.particle.acceleration = self.all_forces / self.PARAMETERS["mass_density"]
        self.particle.next_acceleration = self.particle.acceleration

        self.choose_collision_types("Cuboid", "Normal")
        self.choose_time_stepping(self.time_stepping)

        self.adapt_to_CFL()