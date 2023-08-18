import math as m
import numpy as np
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
                tank_attrs:dict = None,
                num_particles:int = None,
                delta_time:int = 0.02):
        
        super().__init__(particle=particle,
                         search_method=search_method,
                         all_particles=all_particles,
                         num_particles=num_particles,
                         hash_table=hash_table,
                         hash_value=hash_value,
                         time_stepping=time_stepping,
                         tank_attrs=tank_attrs,
                         delta_time=delta_time)
    
    def normalize(self, vector):

        if np.linalg.norm(vector) == 0:
            return 0
        
        return vector / np.linalg.norm(vector)
        
    def new_cubic_spline_kernel_gradient(self, position:np.array=None):
        if np.linalg.norm(position) == 0: q = 0
        q = np.linalg.norm(position) / self.PARAMETERS["cell_size"]
        if q>=0 and q<1:
            kernel_val = m.pow(1 - q, 2) * self.normalize(position)
            return kernel_val
        if q>=1:
            return 0
    
    # ---------------------------------------------------------------------- VORTICITY REFINEMENT --------------------------------------------------------------------------
        
    def vorticity(self, particle):
        for nbr_particle in particle.neighbour_list:
            try:
                curl_vector = np.cross((particle.predicted_velocity - nbr_particle.predicted_velocity),
                                self.new_cubic_spline_kernel_gradient(particle.initial_pos - nbr_particle.initial_pos))
            except np.AxisError:
                curl_vector = np.array([0, 0, 0], dtype="float64")

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
        self.vorticity_grad = np.array([0, 0, 0], dtype="float64")
        for nbr_particle in particle.neighbour_list:
            self.vorticity_grad += (
                (nbr_particle.mass / nbr_particle.mass_density)*
                (nbr_particle.predicted_velocity - self.particle.predicted_velocity)*
                self.new_cubic_spline_kernel_gradient( particle.initial_pos - nbr_particle.initial_pos)
            )
        self.vorticity_grad *= particle.vorticity

    def update_vorticity_laplacian(self, particle):
        
        dimension_const = 10
        viscosity_force = self.PARAMETERS["viscosity"]
        self.vorticity_laplacian = np.array([0, 0, 0], dtype="float64")

        for nbr_particle in particle.neighbour_list:
            self.vorticity_laplacian += (
                (nbr_particle.mass / nbr_particle.mass_density)*
                ((particle.vorticity - nbr_particle.vorticity)*
                 (particle.initial_pos - nbr_particle.initial_pos) / (
                    (np.power(particle.initial_pos - nbr_particle.initial_pos,2) +
                    0.01*m.pow(self.PARAMETERS["cell_size"], 2))
                ))*
                self.new_cubic_spline_kernel_gradient(particle.initial_pos - nbr_particle.initial_pos)
            )
        self.vorticity_laplacian *= dimension_const*viscosity_force
    
    def update_vorticity_dt(self, particle):

        self.update_vorticity_velocity_grad(particle)
        self.update_vorticity_laplacian(particle)
        return self.vorticity_grad + self.vorticity_laplacian

    def update_volume(self, particle):
        particle.volume = particle.mass / particle.mass_density
    
    # ------------------------------------------------------------------- STREAM RELATED FUNCTIONS ------------------------------------------------------------------------

    def update_stream_function(self, particle):

        size_const = 1 / (4*np.pi)
        particle.stream_function = np.array([0, 0, 0], dtype="float64")
        for nbr_particle in particle.neighbour_list:
            if np.linalg.norm(particle.initial_pos - nbr_particle.initial_pos) < self.PARAMETERS["cell_size"]:
                volume_term = (nbr_particle.vorticity_del*nbr_particle.volume)
                particle.stream_function += (
                        volume_term / np.linalg.norm(particle.initial_pos - nbr_particle.initial_pos)
                )
        particle.stream_function *= size_const

    def update_velocity_del(self):

        self.update_vorticity(self.particle, 3)

        self.velocity_del = np.array([0, 0, 0], dtype="float64")
        for nbr_particle in self.particle.neighbour_list:
            try:
                stream_curl = np.cross((self.particle.stream_function - nbr_particle.stream_function),
                            self.new_cubic_spline_kernel_gradient(self.particle.initial_pos -  nbr_particle.initial_pos))
            except np.AxisError:
                stream_curl = np.array([0, 0, 0], dtype="float64")
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