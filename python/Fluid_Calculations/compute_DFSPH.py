import math as m
import numpy as np
import random as rd
import re

from compute_SPH import SPH
from Particles.particles import Particle

class DFSPH(SPH):

    OTHER_PARAMS = {
        "boundary_threshold":6,
        "max_iter":1,
        "stiffness_constant":1000,
        "alpha":0.4,
        "beta_const":0.25,
        "lambda_const":0.005,
        "maximum_force":1000
    }

    def __init__(self,
                particle: Particle=None,
                search_method:str = None,
                hash_table:dict = None,
                hash_value:int = None,
                time_stepping:str = "Euler Cromer",
                tank_attrs:dict = None,
                delta_time:int = 0.02):
        
        super().__init__(particle=particle,
                         search_method=search_method,
                         hash_table=hash_table,
                         hash_value=hash_value,
                         time_stepping=time_stepping,
                         tank_attrs=tank_attrs,
                         delta_time=delta_time)
        
        self.non_pressure_f = np.array([0, 0, 0], dtype="float64")
        self.divergence = np.array([0, 0, 0], dtype="float64")
        self.divergence_force = np.array([0, 0, 0], dtype="float64")
        self.density_error = np.array([0, 0, 0], dtype="float64")

    def update(self):
        
        self.update_non_pressure_f()
        self.adapt_to_CFL()

        self.update_divergence()
        self.correct_density_error()
        self.correct_divergence_error()

        self.acceleration = self.non_pressure_f + self.particle.pressure_force

        self.XSPH_vel_correction()
        self.choose_collision_types()
        self.choose_time_stepping()
    
    def correct_density_error(self):
        iter_step = 0
        if self.mass_density - self.PARAMETERS["mass_density"] > self.OTHER_PARAMS["boundary_threshold"] or \
            iter_step < 1:

            self.update_stiffness_k()
            self.update_density_velocity()

            iter_step += 1
    
    def correct_divergence_error(self):
        iter_step = 0
        while (self.divergence > self.OTHER_PARAMS["boundary_threshold"]) or \
            iter_step < 1:

            self.update_stiffness_k_v()
            self.update_divergence_velocity()

            iter_step +=1 

    def update_non_pressure_f(self):
        
        self.update_viscosity()
        self.update_gravity()
        self.update_buoyancy()
        self.update_surface_tension()

        self.non_pressure_f += (
            self.particle.viscosity +
            self.particle.gravity +
            self.particle.buoyancy +
            self.particle.surface_tension
        )

    def update_divergence_factor(self):
        self.particle.divergence_factor = (
            self.mass_density_squared() + self.sum_mass_density_squared()
        )

    def update_stiffness_k_v(self):
        self.particle.stiffness_k_v = (
            1 / self.delta_time * 
            self.divergence * 
            1/self.particle.divergence_factor
        )

    def update_stiffness_k(self):
        self.update_density_error()
        self.particle.stiffness_k = (
            1/ m.pow(self.delta_time, 2) *
            self.density_error * self.particle.divergence_factor
        )

    def update_divergence(self):
        for nbr_particle in self.neighbours_list:
            self.divergence += (
                nbr_particle.mass * (
                self.particle.velocity - nbr_particle.velocity
                ) *
                (
                    self.cubic_spline_kernel(kernel_type=1, 
                                             nbr_position=nbr_particle.initial_pos)
                )
            )
    
    def mass_density_squared(self):
        squared_mass_density = 0
        for nbr_particle in self.neighbours_list:
            squared_mass_density += (
                nbr_particle.mass * self.cubic_spline_kernel(
                nbr_position=nbr_particle.initial_pos,
                kernel_type=1
                )
            )
        return np.linalg.norm(squared_mass_density)
    
    def sum_mass_density_squared(self):
        squared_mass_density = 0
        for nbr_particle in self.neighbours_list:
            normalized_val = np.linalg.norm((
                nbr_particle.mass * self.cubic_spline_kernel(
                nbr_position=nbr_particle.initial_pos,
                kernel_type=1
                )
            ))
            squared_mass_density += normalized_val
        return squared_mass_density
    
    def update_density_error(self):
        density_error = np.array([0, 0, 0], dtype="float64")
        for nbr_particle in self.neighbours_list:
            density_error += (
                nbr_particle.mass *
                (self.particle.pressure_force/ self.particle.mass - 
                 self.divergence_force / self.particle.mass) *
                 self.cubic_spline_kernel(
                    kernel_type=1,
                    nbr_position=nbr_particle.initial_pos
                 )
            )
        self.density_error = m.pow(self.delta_time, 2)*density_error

    
    def update_divergence_velocity(self):
        divergence_velocity = np.array([0, 0, 0], dtype="float64")
        for nbr_particle in self.neighbours_list:
            divergence_velocity += (
                nbr_particle.mass  *
                (self.particle.stiffness_k_v / self.particle.mass_density + 
                 nbr_particle.stiffness_k_v / nbr_particle.mass_density) *
                 self.cubic_spline_kernel(
                    kernel_type=1,
                    nbr_position=nbr_particle.initial_pos
                 )
            )
        self.particle.velocity -= (
            self.delta_time*divergence_velocity
        )

    def update_density_velocity(self):
        divergence_density = np.array([0, 0, 0], dtype="float64")
        for nbr_particle in self.neighbours_list:
            divergence_density += (
                nbr_particle.mass  *
                (self.particle.stiffness_k / self.particle.mass_density + 
                 nbr_particle.stiffness_k / nbr_particle.mass_density) *
                 self.cubic_spline_kernel(
                    kernel_type=1,
                    nbr_position=nbr_particle.initial_pos
                 )
            )
        self.particle.velocity -= (
            m.pow(self.delta_time, 2) * divergence_density
        )

    def update_pressure_force(self):
        for nbr_particle in self.neighbours_list:
            self.particle.pressure_force += (
                nbr_particle.mass  *
                (self.particle.stiffness_k_v / self.particle.mass_density + 
                 nbr_particle.stiffness_k_v / nbr_particle.mass_density) *
                 self.cubic_spline_kernel(
                    kernel_type=1,
                    nbr_position=nbr_particle.initial_pos
                 )
            )
        self.particle.pressure_force *= -self.particle.mass
    
    def CFL_condition(self):
        return self.OTHER_PARAMS["alpha"] * \
                (self.PARAMETERS["cell_size"] / np.sqrt(self.OTHER_PARAMS["stiffness_constant"]))
    
    def CFL_force_condition(self):
        return self.OTHER_PARAMS["beta_const"] * \
                np.sqrt((self.PARAMETERS["cell_size"]/ self.OTHER_PARAMS["maximum_force"]))
    
    def CFL_viscosity_condition(self):
        return self.OTHER_PARAMS["lambda_const"] / np.linalg.norm(self.velocity_divergence)
    
    def adapt_to_CFL(self):

        self.CFL_conditions = [self.CFL_condition(), self.CFL_force_condition(),
                               self.CFL_viscosity_condition()]
        self.delta_time = min(self.CFL_conditions)

