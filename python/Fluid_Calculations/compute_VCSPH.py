import math as m
import numpy as np
import cupy as cp
import random as rd
import re
import time

from Fluid_Calculations.compute_DFSPH import DFSPH
from Particles.particles import Particle

class VCSPH(DFSPH):

    
    def __init__(self,
                particle: Particle=None,
                search_method:str = None,
                hash_table:dict = None,
                all_particles:list = None,
                hash_value:int = None,
                time_stepping:str = "Euler Cromer",
                time_schemes:dict = None,
                params:dict = None,
                collision_types:dict = None,
                tank_attrs:dict = None,
                num_particles:int = None,
                delta_time:int = 0.02):
        
        super().__init__(particle=particle,
                         search_method=search_method,
                         all_particles=all_particles,
                         num_particles=num_particles,
                         time_schemes=time_schemes,
                         params=params,
                         collision_types=collision_types,
                         hash_table=hash_table,
                         hash_value=hash_value,
                         time_stepping=time_stepping,
                         tank_attrs=tank_attrs,
                         delta_time=delta_time)
    
    # ---------------------------------------------------------------------- VORTICITY REFINEMENT --------------------------------------------------------------------------
        
    def vorticity(self, particle):
        for nbr_particle in particle.neighbour_list:
            try:
                curl_vector = cp.cross((particle.predicted_velocity - nbr_particle.predicted_velocity),
                                self.new_cubic_spline_kernel_gradient(particle.initial_pos - nbr_particle.initial_pos))
            except cp.AxisError:
                curl_vector = cp.array([0, 0, 0], dtype="float64")

            if any(curl_vector) == 0:
                particle.vorticity += 0
            else:
                particle.vorticity += (
                    (nbr_particle.mass / nbr_particle.mass_density)*curl_vector
                )

    # --------------------------------------------------------------------- VORTICITY CALCULATIONS ------------------------------------------------------------------------

    def update_vorticity(self, particle, depth:int = 3):

        if depth==0:
            return
        
        self.vorticity(particle)
        particle.vorticity_new = particle.vorticity + self.delta_time*self.update_vorticity_dt(particle)
        self.update_vorticity_del(particle)
        self.update_volume(particle)
        self.update_stream_function(particle)

        for nbr in particle.neighbour_list:
            self.update_vorticity(nbr, depth-1)

    def update_vorticity_del(self, particle):
        particle.vorticity_del = particle.vorticity_new - particle.vorticity

    def update_vorticity_velocity_grad(self, particle):
        self.vorticity_grad = cp.array([0, 0, 0], dtype="float64")
        for nbr_particle in particle.neighbour_list:
            try:
                if nbr_particle.mass_density==0:
                    mass_d = 0
                else:
                    mass_d = nbr_particle.mass / nbr_particle.mass_density
            except ZeroDivisionError:
                mass_d = 0

            self.vorticity_grad += (
                mass_d*(nbr_particle.predicted_velocity - self.particle.predicted_velocity)*
                self.new_cubic_spline_kernel_gradient( particle.initial_pos - nbr_particle.initial_pos)
            )
        self.vorticity_grad *= particle.vorticity

    def update_vorticity_laplacian(self, particle):
        
        dimension_const = 10
        viscosity_force = self.PARAMETERS["viscosity"]
        self.vorticity_laplacian = cp.array([0, 0, 0], dtype="float64")

        for nbr_particle in particle.neighbour_list:
            try:
                if nbr_particle.mass_density==0:
                    mass_d = 0
                else:
                    mass_d = nbr_particle.mass / nbr_particle.mass_density
            except ZeroDivisionError:
                mass_d = 0
            self.vorticity_laplacian += (
                mass_d * ((particle.vorticity - nbr_particle.vorticity)*
                 (particle.initial_pos - nbr_particle.initial_pos) / (
                    (cp.power(particle.initial_pos - nbr_particle.initial_pos,2) +
                    0.01*m.pow(self.PARAMETERS["cell_size"], 2))
                ))*self.new_cubic_spline_kernel_gradient(particle.initial_pos - nbr_particle.initial_pos)
            )
        self.vorticity_laplacian *= dimension_const*viscosity_force
    
    def update_vorticity_dt(self, particle):

        self.update_vorticity_velocity_grad(particle)
        self.update_vorticity_laplacian(particle)
        return self.vorticity_grad + self.vorticity_laplacian

    def update_volume(self, particle):
        try:
            if particle.mass_density == 0:
                particle.volume = 0
            else:
                particle.volume = particle.mass / particle.mass_density
        except ZeroDivisionError:
            particle.volume = 0
    # ------------------------------------------------------------------- STREAM RELATED FUNCTIONS ------------------------------------------------------------------------

    def update_stream_function(self, particle):

        size_const = 1 / (4*cp.pi)
        particle.stream_function = cp.array([0, 0, 0], dtype="float64")
        for nbr_particle in particle.neighbour_list:
            if cp.linalg.norm(particle.initial_pos - nbr_particle.initial_pos) < self.PARAMETERS["cell_size"]:
                volume_term = (nbr_particle.vorticity_del*nbr_particle.volume)
                if cp.linalg.norm(particle.initial_pos - nbr_particle.initial_pos)==0:
                    particle.stream_function += cp.array([0, 0, 0], dtype="float64")
                else:
                    particle.stream_function += (
                            volume_term / cp.linalg.norm(particle.initial_pos - nbr_particle.initial_pos)
                    )
        particle.stream_function *= size_const

    def update_velocity_del(self):

        self.update_vorticity(self.particle, 3)

        self.velocity_del = cp.array([0, 0, 0], dtype="float64")
        for nbr_particle in self.particle.neighbour_list:
            try:
                stream_curl = cp.cross((self.particle.stream_function - nbr_particle.stream_function),
                            self.new_cubic_spline_kernel_gradient(self.particle.initial_pos -  nbr_particle.initial_pos))
            except cp.AxisError:
                stream_curl = cp.array([0, 0, 0], dtype="float64")
            self.velocity_del += (
                (nbr_particle.mass / nbr_particle.mass_density) * stream_curl
                )

    def update_velocity(self):

        self.update_velocity_del()
        self.particle.velocity = self.particle.predicted_velocity + self.OTHER_PARAMS["alpha_vorticity"]*self.velocity_del

    def debugging_forces(self, secs):
        print("Volume is:", self.particle.volume)
        print("Vorticity is:", self.particle.vorticity)
        print("Vorticity del is:", self.particle.vorticity_del)
        print("Stream Function is:", self.particle.stream_function)
        print("\n\n")
        time.sleep(secs)
    
    # --------------------------------------------------------------------------- UPDATE ------------------------------------------------------------------------------------

    def update(self):

        super().update_errors()
        self.update_velocity()

        self.XSPH_vel_correction()
        self.choose_collision_types("Cuboid", "Normal")