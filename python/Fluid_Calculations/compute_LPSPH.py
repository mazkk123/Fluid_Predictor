import math as m
import numpy as np
import random as rd
import re
import time

from Particles.particles import Particle
from Fluid_Calculations.compute_sph import SPH

class LPSPH(SPH):

    OTHER_PARAMS = {
        "max_iterations":3
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
        
        self.ref_radius = self.PARAMETERS["cell_size"] / 4
        self.epsilon = 2/3*self.ref_radius

        self.near_particles = []
        self.far_particles = []

        self.pressure_far = np.array([0, 0, 0], dtype="float64")
        self.pressure_near = np.array([0, 0, 0], dtype="float64")
        
        self.pressure = np.array([0, 0, 0], dtype="float64")
        self.pressure_force = np.array([0, 0, 0], dtype="float64")
        
        self.pressure_neighbourhood()
        
        self.iterations = 0

    # ------------------------------------------------------------- PRESSURE CORRECTION ------------------------------------------------------------------------
    def pressure_neighbourhood(self):
        for id, nbr_particle in enumerate(self.neigbours_list):
            if self.is_near_particle(nbr_particle.initial_pos):
                self.near_particles.append(nbr_particle)
            if self.is_far_particle(nbr_particle.initial_pos):
                self.far_particles.append(nbr_particle)
    
    def pressure_neighbourhoods(self, particle):
        for id, nbr_particle in enumerate(particle.neigbours_list):
            if self.is_near_particle_nbr(nbr_particle.predicted_initial_pos):
                particle.near_particles.append(nbr_particle)
                if self.is_far_particle_nbr(nbr_particle.predicted_initial_pos):
                particle.far_particles.append(nbr_particle)
                
    def is_near_particle(self, position: np.array=None):
        if position is not None:
            return self.epsilon <= np.linalg.norm(self.particle.initial_pos - position )

    def is_far_particle(self, position: np.array=None):
        if position is not None:
            return self.epsilon > np.linalg.norm(self.particle.initial_pos - position )
    
    def is_near_particle_nbr(self, position: np.array=None, particle):
        if position is not None:
            return self.epsilon <= np.linalg.norm(particle.prdicted_initial_pos - position )

    def is_far_particle_nbr(self, position: np.array=None, particle):
        if position is not None:
            return self.epsilon > np.linalg.norm(particle.predicted_initial_pos - position )
            
    def near_pressure(self, particle):
        pressure_near = 0
        for id, near_nbr in enumerate(particle.near_particles):
            density_diff = (near_nbr.predicted_density - self.PARAMETERS["mass_density"]) * near_nbr.predicted_density * \
                            m.pow(self.ref_radius, 2)
            div_const = 2 *self.PARAMETERS["mass_density"] * m.pow(self.delta_time, 2)
            try:
                pressure_near = density_diff/div_const
            except ZeroDivisionError:
                pressure_near = 0
                
            self.pressure_near += pressure_near
        return self.pressure_near
    
    def far_pressure(self, particle):
        density_term = 0
        for id, far_nbr in enumerate(particle.far_particles):
            density_diff = far_nbr.predicted_density - self.PARAMETERS["mass_density"]
            div_const = 4 * np.pi * self.PARAMETERS["mass_density"] *  \
                        m.pow(self.delta_time, 2) * np.linalg.norm(particle.predicted_initial_pos - 
                                                                   far_nbr.predicted_initial_pos)
            try:
                density_term = density_diff/div_const
            except ZeroDivisionError:
                density_term = 0

            self.pressure_far += far_nbr.mass*density_term
        return self.pressure_far

    def calculate_density_error(self):
        return self.particle.predicted_density - self.PARAMETERS["mass_density"]
    
    def update_pressure(self, particle):
        particle.pressure = self.near_pressure(particle) + self.far_pressure(particle)
    
    def update_pressure_force(self, particle):
        pressure_force = np.array([0, 0, 0], dtype="float64")
        for nbr_particle in particle.neighbour_list:
            pressure_force += (
                nbr_particle.mass *
                ((nbr_particle.pressure + particle.pressure) / 2*nbr_particle.pressure ) *
                self.cubic_spline_kernel_gradient(particle.predicted_initial_pos - nbr_particle.predicted_initial_pos)
            )
        self.pressure_force = pressure_force*-1

    # ----------------------------------------------------------------- PREDICTING ATTRIBS ------------------------------------------------------------------------

    def update_predicted_density(self, particle):
        density = 0 
        for nbr_particle in particle.neighbour_list:
            kernel_value = self.kernel_linear(particle.predicted_initial_pos - nbr_particle.predicted_initial_pos, 0)
            density += kernel_value*nbr_particle.mass

        particle.predicted_density = self.PARAMETERS["mass_density"] + density

    def update_predicted_mass_density(self, particle):
        """
        """
        density = 0 
        for nbr_particle in particle.neighbour_list:
            kernel_value = self.kernel_linear(particle.initial_pos - nbr_particle.initial_pos, 0)
            density += kernel_value*nbr_particle.mass

        particle.mass_density = self.PARAMETERS["mass_density"] + density

    def update_predicted_viscosity(self, particle):
        """
        """
        viscosity = np.array([0, 0, 0], dtype="float64")
        for id, nbr_particle in enumerate(particle.neighbour_list):
            vel_dif = nbr_particle.velocity - particle.velocity
            kernel_laplacian = self.kernel_laplacian(particle.initial_pos - nbr_particle.initial_pos, 2)
            try:
                mass_pressure = particle.mass/nbr_particle.mass_density
            except ZeroDivisionError:
                mass_pressure = 0
            viscosity += vel_dif*mass_pressure*kernel_laplacian

        particle.viscosity = viscosity*self.PARAMETERS["viscosity"]

    def update_predicted_gravity(self, particle):
        """
        """
        particle.gravity = self.gravity_const
    
    def update_predicted_normal_field(self, particle):
        """
        """
        normal_field = np.array([0, 0, 0],  dtype="float64")
        for id, nbr_particle in enumerate(particle.neighbour_list):
            pos_difference = particle.initial_pos - nbr_particle.initial_pos
            normal_field += (
                particle.mass * 1/particle.mass_density * self.kernel_gradient(pos_difference, 0)
            )
        return normal_field

    def update_predicted_surface_curvature(self, particle):
        """        
        """
        surface_curvature = 0
        for id, nbr_particle in enumerate(particle.neighbour_list):
            surface_curvature += (
                particle.mass * 1/particle.mass_density * self.kernel_laplacian(
                particle.initial_pos - nbr_particle.initial_pos, 0)
            )
        return surface_curvature
    
    def update_predicted_surface_tension(self, particle):
        """
        """
        normal_field = self.update_predicted_normal_field(particle)
        surface_curvature = self.update_predicted_surface_curvature(particle)
        normal_field_magnitude = np.linalg.norm(normal_field)
        
        if normal_field_magnitude >= self.PARAMETERS["tension_threshold"]:
            particle.surface_tension = (
                self.PARAMETERS["tension_coefficient"] * surface_curvature * normal_field/normal_field_magnitude 
            )

    def update_predicted_buoyancy(self, particle):
        """
        """
        buoyancy = self.PARAMETERS["buoyancy"] * (particle.mass_density - self.PARAMETERS["mass_density"])
        buoyancy *= self.gravity_const
        particle.buoyancy = buoyancy

    def find_neighbour_list(self, particle):
        particle.neighbour_list = []
        for nbr in self.hash_table[particle.hash_value]:
            if particle is not nbr:
                particle.neighbour_list.append(nbr)

    def predict_attrs(self):

        for particle in self.all_particles:

            self.find_neighbour_list(particle)
            self.update_predicted_mass_density(particle)
            self.update_predicted_gravity(particle)
            self.update_predicted_buoyancy(particle)
            self.update_predicted_surface_tension(particle)
            self.update_predicted_viscosity(particle)

            self.all_forces = particle.gravity + \
                              particle.buoyancy + \
                              particle.surface_tension + \
                              particle.viscosity
            
            particle.predicted_velocity = particle.velocity + self.all_forces*self.delta_time
            particle.predicted_initial_pos = particle.initial_pos + particle.predicted_velocity*self.delta_time
            particle.predicted_density = particle.mass_density

    def update_predicted_attrs(self):

        for particle in self.all_particles:
            
            self.update_predicted_density(particle)
            self.update_pressure_neighbourhoods(particle)
            self.update_pressure(particle)
            self.update_pressure_force(particle)
            
            particle.predicted_velocity += self.delta_time * self.pressure_force / particle.mass
            particle.predicted_initial_pos += m.pow(self.delta_time, 2) * (self.pressure_force / particle.mass)
            
    def debug_predicted_attrs(self, secs:float):

        print("Mass Density:", self.particle.predicted_density)
        print("Pressure:", self.particle.pressure)
        print("Pressure Force:", self.particle.pressure_force)
        print("Buoyancy:", self.particle.buoyancy)
        print("surface tension:", self.particle.surface_tension)
        print("viscosity:", self.particle.viscosity)
        print("\n\n")
        time.sleep(secs)
    # --------------------------------------------------------------------- UPDATE STEP ------------------------------------------------------------------------------
    def update(self):

        self.predict_attrs()
        self.debug_predicted_attrs(0)
        
        while self.calculate_density_error()> 0.05*self.PARAMETERS["mass_density"] and \
            self.iterations < self.OTHER_PARAMS["max_iterations"]:
            
            print("Entering pressure correction")
            self.debug_predicted_attrs(0.25)
            self.update_predicted_attrs()
            
            self.iterations += 1
            
        self.particle.initial_pos = self.particle.predicted_initial_pos
        self.particle.velocity = self.particle.predicted_velocity 
        
        self.XSPH_vel_correction()
        self.choose_collision_types("Cuboid", "Normal")