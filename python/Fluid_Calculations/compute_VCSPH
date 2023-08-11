import math as m
import numpy as np
import random as rd
import re

from Fluid_Calculations.compute_DFSPH import DFSPH
from Particles.particles import Particle

class VCSPH(DFSPH):

    
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
        
        self.trapezoidal_amt = 100
        self.trapezoidal_increment = 0.05
        
    def vorticity(self):
        for nbr_particle in self.neighbour_list:
            self.particle.vorticity += (
                (nbr_particle.mass / nbr_particle.mass_density)*
                np.cross((self.particle.velocity - nbr_particle.velocity),
                        self.cubic_spline_kernel_gradient(self.particle.initial_pos - 
                                                          nbr_particle.initial_pos))
            )

    def update_vorticity(self):
        self.particle.vorticity += self.delta_time*self.update_vorticity_dt()

    def update_vorticity_del(self):
        self.update_vorticity()
        self.velocity_curl = np.array([0, 0, 0], dtype="float64")
        for nbr_particle in self.neighbour_list:
            self.velocity_curl += (
                (nbr_particle.mass / nbr_particle.mass_density)*
                -(self.particle.predicted_velocity - nbr_particle.predicted_velocity)*
                self.cubic_spline_kernel_gradient(self.particle.initial_pos - 
                                                  nbr_particle.initial_pos)
            )
        self.particle.vorticity_del = self.particle.vorticity - self.velocity_curl

    def update_vorticity_velocity_grad(self):
        self.vorticity_grad = np.array([0, 0, 0], dtype="float64")
        for nbr_particle in self.neighbour_list:
            self.vorticity_grad += (
                (nbr_particle.mass / nbr_particle.mass_density)*
                (nbr_particle.velocity - self.particle.velocity)*
                self.cubic_spline_kernel_gradient( self.particle.initial_pos - 
                                                  nbr_particle.initial_pos)
            )
        self.vorticity_grad *= self.particle.vorticity

    def update_vorticity_laplacian(self):
        dimension_const = 10
        viscosity_force = self.PARAMETERS["viscosity"]
        self.vorticity_laplacian = np.array([0, 0, 0], dtype="float64")
        for nbr_particle in self.neighbour_list:
            self.vorticity_laplacian += (
                (nbr_particle.mass / nbr_particle.mass_density)*
                ((self.particle.vorticity - nbr_particle.vorticity) / (
                    np.power(self.particle.initial_pos - nbr_particle.initial_pos,2) +
                    0.01*m.pow(self.PARAMETERS["cell_size"], 2)
                ))*
                self.cubic_spline_kernel_gradient(self.particle.initial_pos - 
                                                  nbr_particle.initial_pos)
            )
        self.vorticity_laplacian *= dimension_const*viscosity_force
    
    def update_vorticity_dt(self):

        self.update_vorticity_velocity_grad()
        self.update_vorticity_laplacian()
        return self.vorticity_grad*self.vorticity_laplacian
    
    def update_volume(self):
        trapezoidal_amt = self.trapezoidal_amt
        for nbr_particle in self.neighbour_list:

            self.volume = 0
            self.trapezoidal_amt = trapezoidal_amt
            while (self.trapezoidal_amt > 0):

                if (self.trapezoidal_amt == trapezoidal_amt):
                    self.volume += self.cubic_spline_kernel_gradient(self.particle.initial_pos - 
                                                                              nbr_particle.initial_pos)
                self.volume += 2*self.cubic_spline_kernel_gradient(self.particle.initial_pos - 
                                                                          nbr_particle.initial_pos)
                self.trapezoidal_amt -= 1

            self.particle.volume += self.volume*self.incremental_step/2

    def update_stream_function(self):
        
        self.update_vorticity_del()
        self.update_volume()

        size_const = 1 / (4*np.pi)
        for nbr_particle in self.neighbour_list:
            if np.linalg.norm(self.particle.initial_pos - nbr_particle.initial_pos) < self.PARAMETERS["cell_size"]:
                self.particle.stream_function += (
                    (nbr_particle.vorticity_del*nbr_particle.volume) / 
                    (np.linalg.norm(self.particle.initial_pos - nbr_particle.initial_pos))
                )
        self.particle.stream_function *= size_const

    def update_velocity_del(self):

        self.update_stream_function()
        self.velocity_del = np.array([0, 0, 0], dtype="float64")
        for nbr_particle in self.neighbour_list:
            self.velocity_del += (
                (nbr_particle.mass / nbr_particle.mass_density) *
                np.cross((self.particle.stream_function - nbr_particle.stream_function),
                         self.cubic_spline_kernel_gradient(self.particle.initial_pos - 
                                                           nbr_particle.initial_pos)
                )
            )

    def update_velocity(self):

        self.update_velocity_del()
        self.particle.velocity = self.particle.predicted_velocity + self.OTHER_PARAMS["alpha_vorticity"]*self.velocity_del

    def update(self):

        self.update_non_pressure_f()
        self.adapt_to_CFL()

        self.update_divergence()
        self.correct_density_error()
        self.correct_divergence_error()
        self.update_velocity()

        self.acceleration = self.non_pressure_f + self.particle.pressure_force

        self.XSPH_vel_correction()
        self.choose_collision_types("Cuboid", "Normal")
        self.choose_time_stepping(self.time_stepping)