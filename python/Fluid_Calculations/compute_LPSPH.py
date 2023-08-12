import math as m
import numpy as np
import random as rd
import re

from Particles.particles import Particle
from Fluid_Calculations.compute_sph import SPH

class LPSPH(SPH):

    OTHER_PARAMS = {
        "min_threshold":0.2,
        "max_threshold":0.015
    }

    def __init__(self, ref_radius:float=None,
                 particle: Particle=None,
                 time_stepping:str = "Euler Cromer",
                 search_method: str=None,
                 hash_table:dict=None,
                 hash_value:int=None,
                 all_particles:list = None,
                 tank_attrs:dict = None,
                 delta_time:float = None):
        
        super().__init__(particle=particle,
                         search_method=search_method,
                         hash_table=hash_table,
                         hash_value=hash_value,
                         time_stepping=time_stepping,
                         all_particles=all_particles,
                         tank_attrs=tank_attrs,
                         delta_time=delta_time)

        if ref_radius is not None:
            self.ref_radius = ref_radius
        
        self.ref_radius = self.PARAMETERS["cell_size"] / 2
        self.epsilon = 2/3*self.ref_radius

        self.near_particles = []
        self.far_particles = []

        self.pressure_far = np.array([0, 0, 0], dtype="float64")
        self.pressure_near = np.array([0, 0, 0], dtype="float64")

        self.pressure_neighbourhood()

        self.density_error = 0

    def pressure_neighbourhood(self):
        for id, nbr_particle in enumerate(self.neighbours_list):
            if self.is_near_particle(nbr_particle.initial_pos):
                self.near_particles.append(nbr_particle)
            if self.is_far_particle(nbr_particle.initial_pos):
                self.far_particles.append(nbr_particle)

    def is_near_particle(self, position: np.array=None):
        if position is not None:
            return self.epsilon <= np.linalg.norm(self.particle.initial_pos - position )

    def is_far_particle(self, position: np.array=None):
        if position is not None:
            return self.epsilon > np.linalg.norm(self.particle.initial_pos - position )
        
    def near_pressure(self):
        pressure_near = 0
        for id, near_nbr in enumerate(self.near_particles):
            density_diff = (near_nbr.predicted_density - self.PARAMETERS["mass_density"]) * near_nbr.predicted_density * \
                            m.pow(self.ref_radius, 2)
            div_const = 2 *self.PARAMETERS["mass_density"] * m.pow(self.delta_time, 2)
            try:
                pressure_near = density_diff/div_const
            except ZeroDivisionError:
                pressure_near = 0
                
            self.pressure_near += pressure_near
        return self.pressure_near
    
    def far_pressure(self):
        density_term = 0
        for id, far_nbr in enumerate(self.far_particles):
            density_diff = far_nbr.predicted_density - self.PARAMETERS["mass_density"]
            div_const = 4 * np.pi * self.PARAMETERS["mass_density"] *  \
                        m.pow(self.delta_time, 2) * np.linalg.norm(self.particle.initial_pos - 
                                                                   far_nbr.initial_pos)
            try:
                density_term = density_diff/div_const
            except ZeroDivisionError:
                density_term = 0

            self.pressure_far += density_term
        return self.pressure_far

    def calculate_density_error(self):
        return self.particle.predicted_density - self.PARAMETERS["mass_density"]
    
    def update_pressure(self):
        self.particle.pressure = self.near_pressure() + self.far_pressure()

    def update_mass_density(self):
        mass_density, mass_density_denom = 0, 0
        for nbr_particle in self.neighbours_list:
            try:
                mass_d = nbr_particle.mass / nbr_particle.mass_density
            except ZeroDivisionError:
                mass_d = 0
            mass_density += (
                nbr_particle.mass * self.cubic_spline_kernel_pos(self.particle.initial_pos - 
                                                       nbr_particle.initial_pos)
            )
            mass_density_denom += (
                mass_d * self.cubic_spline_kernel_pos(
                    self.particle.initial_pos - nbr_particle.initial_pos
                )
            )
        try:
            final_mass_d = mass_density / mass_density_denom
        except ZeroDivisionError:
            final_mass_d = 0

        self.particle.mass_density = final_mass_d
    
    def update_predicted_density(self):
        mass_density, mass_density_denom = 0, 0
        for nbr_particle in self.neighbours_list:
            try:
                mass_d = nbr_particle.mass / nbr_particle.mass_density
            except ZeroDivisionError:
                mass_d = 0
            mass_density += (
                nbr_particle.mass * self.cubic_spline_kernel_pos(self.particle.initial_pos - 
                                                       nbr_particle.initial_pos)
            )
            mass_density_denom += (
                mass_d * self.cubic_spline_kernel_pos(
                    self.particle.initial_pos - nbr_particle.initial_pos
                )
            )
        try:
            final_mass_d = mass_density / mass_density_denom
        except ZeroDivisionError:
            final_mass_d = 0

        self.particle.predicted_density = final_mass_d
        
    def update_viscosity(self):
        return super().update_viscosity()
    
    def update_predicted_attrs(self):
        self.particle.predicted_velocity += self.delta_time * self.all_forces / self.particle.mass
        self.particle.predicted_initial_pos += m.pow(self.delta_time, 2) * (self.all_forces / self.particle.mass)

    def update_advective_forces(self):

        self.update_mass_density()
        self.update_gravity()
        self.update_buoyancy()
        self.update_surface_tension()
        self.update_viscosity()

        self.debugging_forces(0.5)

        self.all_forces = self.particle.gravity + self.particle.buoyancy + self.particle.surface_tension + \
                          self.particle.viscosity + self.particle.surface_tension

    def update_pressure_force(self):
        pressure_force = np.array([0, 0, 0], dtype="float64")
        for nbr_particle in self.neighbours_list:
            pressure_force += (
                nbr_particle.mass *
                ((nbr_particle.pressure + self.particle.pressure) / 2*nbr_particle.pressure ) *
                self.cubic_spline_kernel_gradient(self.particle.initial_pos - nbr_particle.initial_pos)
            )
        self.particle.pressure_force = pressure_force*-1

    def update(self):

        self.update_advective_forces()
        self.particle.acceleration = self.all_forces / self.particle.mass

        self.particle.predicted_velocity += self.particle.acceleration
        self.particle.predicted_initial_pos += self.delta_time*self.particle.predicted_velocity
        
        self.particle.predicted_density = self.particle.mass_density
        self.density_error = self.calculate_density_error()
        
        while self.density_error > self.OTHER_PARAMS["max_threshold"]:

            self.update_predicted_density()
            self.update_pressure()

            self.density_error = self.calculate_density_error()

            self.update_pressure_force()

            self.all_forces = self.particle.pressure_force

            self.update_predicted_attrs()
            
        
        self.particle.initial_pos = self.particle.predicted_initial_pos
        self.particle.velocity = self.particle.predicted_velocity 
        
        self.XSPH_vel_correction()
        self.choose_collision_types("Cuboid", "Normal")

