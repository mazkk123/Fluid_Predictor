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
        return vector / np.linalg.norm(vector)
        
    def wedland_4_kernel_gradient(self, position:np.array=None):
        q = np.linalg.norm(position)
        kernel_const = -30 / np.pi*m.pow(self.PARAMETERS["cell_size"], 5)
        kernel_term_1 = m.pow(1 - q / self.PARAMETERS["cell_size"], 4)
        kernel_term_2 = 2*q / self.PARAMETERS["cell_size"] + 1
        r_value = self.normalize(position)
        return kernel_const*kernel_term_1*kernel_term_2*r_value
    
    # ---------------------------------------------------------------------- VORTICITY REFINEMENT --------------------------------------------------------------------------
        
    def vorticity(self, particle):
        for nbr_particle in particle.neighbour_list:
            particle.vorticity += (
                (nbr_particle.mass / nbr_particle.mass_density)*
                np.cross((particle.predicted_velocity - nbr_particle.predicted_velocity),
                        self.wedland_4_kernel_gradient(particle.initial_pos - nbr_particle.initial_pos))
            )

    def update_vorticity(self):

        for nbr in self.particle.neighbour_list:
            self.find_neighbour_list(nbr)
            self.vorticity(nbr)
        
        self.vorticity(self.particle)

        for nbr in self.particle.neighbour_list:
            nbr.vorticity_new += self.delta_time*self.update_vorticity_dt(nbr)

        self.particle.vorticity_new += self.delta_time*self.update_vorticity_dt(self.particle)

    # --------------------------------------------------------------------- VORTICITY CALCULATIONS ------------------------------------------------------------------------

    def update_vorticity_del(self):

        self.update_vorticity()

        for nbr in self.particle.neighbour_list:
            nbr.vorticity_del = nbr.vorticity_new - nbr.vorticity

        self.particle.vorticity_del = self.particle.vorticity_new - self.particle.vorticity

    def update_vorticity_velocity_grad(self, particle):
        self.vorticity_grad = np.array([0, 0, 0], dtype="float64")
        for nbr_particle in particle.neighbour_list:
            self.vorticity_grad += (
                (nbr_particle.mass / nbr_particle.mass_density)*
                (nbr_particle.predicted_velocity - self.particle.predicted_velocity)*
                self.wedland_4_kernel_gradient( particle.initial_pos - nbr_particle.initial_pos)
            )
        self.vorticity_grad *= particle.vorticity

    def update_vorticity_laplacian(self, particle):
        
        dimension_const = 10
        viscosity_force = self.PARAMETERS["viscosity"]
        self.vorticity_laplacian = np.array([0, 0, 0], dtype="float64")

        for nbr_particle in particle.neighbour_list:
            self.vorticity_laplacian += (
                (nbr_particle.mass / nbr_particle.mass_density)*
                ((particle.vorticity - nbr_particle.vorticity) / (
                    np.power(particle.initial_pos - nbr_particle.initial_pos,2) +
                    0.01*m.pow(self.PARAMETERS["cell_size"], 2)
                ))*
                self.wedland_4_kernel_gradient(particle.initial_pos - nbr_particle.initial_pos)
            )
        self.vorticity_laplacian *= dimension_const*viscosity_force
    
    def update_vorticity_dt(self, particle):

        self.update_vorticity_velocity_grad(particle)
        self.update_vorticity_laplacian(particle)
        return self.vorticity_grad*self.vorticity_laplacian
    
    # ------------------------------------------------------------------- STREAM RELATED FUNCTIONS ------------------------------------------------------------------------

    def update_volume(self):
        
        for nbr in self.particle.neighbour_list:
            nbr.volume = nbr.mass / nbr.mass_density
            
        self.particle.volume = self.particle.mass / self.particle.mass_density

    def update_stream_function(self, particle):
        
        self.update_vorticity_del()
        self.update_volume()

        size_const = 1 / (4*np.pi)
        particle.stream_function = np.array([0, 0, 0], dtype="float64")
        for nbr_particle in particle.neighbour_list:
            if np.linalg.norm(particle.initial_pos - nbr_particle.initial_pos) < self.PARAMETERS["cell_size"]:
                particle.stream_function += (
                    (nbr_particle.vorticity_del*nbr_particle.volume) / 
                    (np.linalg.norm(particle.initial_pos - nbr_particle.initial_pos))
                )
        particle.stream_function *= size_const

    def update_velocity_del(self):

        for nbr in self.particle.neighbour_list:
            self.update_stream_function(nbr)
        
        self.update_stream_function(self.particle)

        self.velocity_del = np.array([0, 0, 0], dtype="float64")
        for nbr_particle in self.particle.neighbour_list:
            self.velocity_del += (
                (nbr_particle.mass / nbr_particle.mass_density) *
                np.cross((self.particle.stream_function - nbr_particle.stream_function),
                         self.wedland_4_kernel_gradient(self.particle.initial_pos -  nbr_particle.initial_pos)
                )
            )

    def update_velocity(self):

        self.update_velocity_del()
        self.particle.velocity = self.particle.predicted_velocity + self.OTHER_PARAMS["alpha_vorticity"]*self.velocity_del

    # --------------------------------------------------------------------------- UPDATE ------------------------------------------------------------------------------------

    def update(self):

        super().update_errors()
        self.update_velocity()

        self.XSPH_vel_correction()
        self.choose_collision_types("Cuboid", "Normal")