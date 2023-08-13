import math as m
import numpy as np
import random as rd
import re
import time

from Fluid_Calculations.compute_sph import SPH
from Particles.particles import Particle

class PCSPH(SPH):

    OTHER_PARAMS = {
        "max_iterations":2,
        "min_iterations":1,
        "stiffness_k":0.1,
        "lambda_const":7
    }

    def __init__(self,
                 particle: Particle=None,
                 search_method: str=None,
                 hash_table:dict=None,
                 hash_value:int=None,
                 time_stepping:str = "Euler Cromer",
                 all_particles:list = None,
                 tank_attrs:dict = None,
                 delta_time:float = 0.02):
        
        super().__init__(particle=particle,
                        all_particles=all_particles, 
                        time_stepping=time_stepping,
                        search_method=search_method,
                        hash_table=hash_table,
                        hash_value=hash_value,
                        tank_attrs=tank_attrs,
                        delta_time=delta_time)
    
    def find_beta_const(self):
        return m.pow(self.delta_time, 2)*m.pow(self.particle.mass, 2)* \
            2 / m.pow(self.PARAMETERS["mass_density"], 2)
    
    def cubic_spline_kernel_gradient(self, position):
        q = np.linalg.norm(position) / self.PARAMETERS["cell_size"]
        kernel_const = 1 / (np.pi*m.pow(self.PARAMETERS["cell_size"], 4))
        if q>=0 and q<=1:
            kernel_val = 9/4*m.pow(q, 2) - 3*q
        if q>=1 and q<=2:
            kernel_val = -3/4*m.pow((2 - q), 2)
        if q>=2:
            return 0
        return kernel_val*kernel_const


    # ---------------------------------------------------------------------- PREDICT DENSITY -----------------------------------------------------------------------------

    def update_del_x(self):
        pressure_const = 2*self.particle.pressure_correction / m.pow(self.PARAMETERS["mass_density"], 2)
        constant = -1*m.pow(self.delta_time, 2)*self.particle.mass*pressure_const

        for nbr_particle in self.neighbours_list:
            self.particle.delta_x += (
                self.cubic_spline_kernel_gradient(self.particle.predicted_initial_pos - nbr_particle.initial_pos)
            )
        self.particle.delta_x *= constant

    def update_density_change(self):
        sum_del_neighbour_x, sum_weights = 0, 0
        self.particle.density_change = 0

        for nbr_particle in self.neighbours_list:
            sum_weights += (
                self.cubic_spline_kernel_gradient(self.particle.predicted_initial_pos - 
                                                  nbr_particle.initial_pos)
            )
            sum_del_neighbour_x += (
                self.cubic_spline_kernel_gradient(self.particle.predicted_initial_pos - 
                                                  nbr_particle.initial_pos)*nbr_particle.delta_x
            )
        self.particle.density_change = self.particle.mass*(self.particle.delta_x*sum_weights - sum_del_neighbour_x)

    def update_next_density(self):
        self.update_del_x()
        self.update_density_change()
        self.particle.predicted_density += self.particle.density_change

    def calculate_density_error(self):
        return self.particle.predicted_density - self.PARAMETERS["mass_density"]

    # -------------------------------------------------------------------- COMPUTE PRESSURE -----------------------------------------------------------------------------

    def find_kronecker_delta(self):
        beta_value = self.find_beta_const()
        accum_denom, denom = 0, 0
        for nbr_particle in self.neighbours_list:
            accum_denom += (
                m.pow(self.cubic_spline_kernel_gradient(self.particle.predicted_initial_pos - nbr_particle.initial_pos), 2)
            )
            denom += (
                self.cubic_spline_kernel_gradient(self.particle.predicted_initial_pos - nbr_particle.initial_pos)
            )

        try:
            kronecker = -1 / (
                beta_value * (
                    -denom*denom - accum_denom
                )
            )
        except ZeroDivisionError:
            kronecker = 0

        return kronecker

    def update_pressure_correction(self):
        self.particle.pressure_correction += self.find_kronecker_delta()*self.calculate_density_error()
    
    def update_pressure(self):
        self.particle.pressure = (
            (self.OTHER_PARAMS["k_stiffness"]*self.PARAMETERS["mass_density"] /
            self.OTHER_PARAMS["lambda_const"]) *
            (m.pow(self.particle.predicted_density / self.PARAMETERS["mass_density"], self.OTHER_PARAMS["lambda_const"]) - 1)
        )
    
    def update_pressure_force(self):
        for nbr_particle in self.neighbours_list:
            self.particle.pressure_force += (
                ((self.particle.pressure / np.power(self.particle.predicted_density, 2)) + 
                 (nbr_particle.pressure / np.power(nbr_particle.predicted_density,2))) *
                  self.cubic_spline_kernel_gradient(self.particle.predicted_initial_pos -
                                                    nbr_particle.predicted_initial_pos)
            )
        self.particle.pressure_force *= m.pow(self.particle.mass, 2)

    # -------------------------------------------------------------------- UPDATE CALLS ------------------------------------------------------------------------------
    
    def update_mass_density(self, particle):
        """
        """
        density = 0 
        for id, nbr_particle in enumerate(particle.neighbour_list):
            kernel_value = self.kernel_linear(particle.initial_pos - nbr_particle.initial_pos, 0)
            density += kernel_value*nbr_particle.mass

        particle.mass_density = self.PARAMETERS["mass_density"] + density
    
    def update_viscosity(self, particle):
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

    def update_gravity(self, particle):
        """
        """
        particle.gravity = self.gravity_const
    
    def update_surface_tension(self, particle):
        """
        """
        normal_field = self.update_normal_field(particle)
        surface_curvature = self.update_surface_curvature(particle)
        normal_field_magnitude = np.linalg.norm(normal_field)
        
        if normal_field_magnitude >= self.PARAMETERS["tension_threshold"]:
            particle.surface_tension = (
                self.PARAMETERS["tension_coefficient"] * surface_curvature * normal_field/normal_field_magnitude 
            )

    def update_buoyancy(self, particle):
        """
        """
        buoyancy = self.PARAMETERS["buoyancy"] * (particle.mass_density - self.PARAMETERS["mass_density"])
        buoyancy *= self.gravity_const
        particle.buoyancy = buoyancy
    
    def update_predicted_attrs(self):
        
        for particle in self.all_particles:
            self.update_mass_density(particle)
            self.update_gravity(particle)
            self.update_surface_tension(particle)
            self.update_viscosity(particle)
            self.update_buoyancy(particle)
    
            self.all_forces = particle.gravity + \
                              particle.surface_tension + \
                              particle.viscosity + \
                              particle.buoyancy
                            
            particle.acceleration = self.all_forces / particle.mass

            particle.predicted_velocity = particle.velocity + self.delta_time*particle.acceleration
            particle.predicted_initial_pos = particle.predicted_velocity + particle.predicted_velocity*self.delta_time

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
        
        iterations = 0
        
        self.particle.pressure_force = np.array([0, 0, 0], dtype="float64")
        self.particle.pressure = 0

        self.update_predicted_attrs()

        while (self.calculate_density_error() > 0.1*self.PARAMETERS["mass_density"]) or \
            iterations < self.OTHER_PARAMS["max_iterations"]:

            self.update_pressure_correction()
            self.update_next_density()
            
            iterations += 1

        self.update_all_forces()

        self.particle.acceleration = self.all_forces / self.particle.mass

        self.choose_time_stepping()
        self.adapt_to_CFL()        
        self.XSPH_vel_correction()
        self.choose_collision_types()