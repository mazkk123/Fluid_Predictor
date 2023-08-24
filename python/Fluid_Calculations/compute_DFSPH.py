import math as m
import numpy as np
import random as rd
import re
import time

from Fluid_Calculations.compute_sph import SPH
from Particles.particles import Particle
from  Fluid_Utilities.search_methods import SpatialHashing

class DFSPH(SPH):

    OTHER_PARAMS = {
        "max_iterations":6,
        "max_iterations_div":6,
        "divergence_error":0.6,
        "alpha_vorticity":0.5
    }

    def __init__(self, 
                 particle: Particle=None,
                 search_method: str=None,
                 all_particles:list = None,
                 time_schemes:dict = None,
                 collision_types:dict = None,
                 params:dict = None,
                 hash_table:dict=None,
                 time_stepping:str = "Euler Cromer",
                 tank_attrs:dict = None,
                 hash_value:int=None,
                 delta_time:float = None,
                 num_particles:int = None):
        
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
                         delta_time=delta_time)
        
        self.divergence = np.array([0, 0, 0], dtype="float64")
        self.divergence_force = np.array([0, 0, 0], dtype="float64")
        self.num_particles = num_particles

        self.predict_advective_forces(self.particle, 4)
        self.update_divergence_factor(self.particle)

    # ----------------------------------------------------------------- UTILITY FUNCTIONS -----------------------------------------------------------------------------

    def predict_advective_forces(self, particle, depth:int=2):
        
        if depth==0:
            return
        
        self.find_neighbour_list(particle)
        self.update_predicted_mass_density(particle)
        self.update_predicted_viscosity(particle)
        self.update_predicted_gravity(particle)
        self.update_predicted_buoyancy(particle)
        self.update_predicted_surface_tension(particle)

        self.all_forces += (
            particle.viscosity +
            particle.buoyancy +
            particle.gravity +
            particle.surface_tension
        )

        particle.predicted_velocity = particle.velocity + self.delta_time*(self.all_forces/particle.mass)
        
        for nbr in particle.neighbour_list:
            return self.predict_advective_forces( nbr , depth-1 )
        
    def update_divergence_factor(self, particle):
        particle.divergence_factor = 0
        pos_diff_sum, sum_pos_diff_norm = 0, 0 
        for nbr_particle in particle.neighbour_list:
            pos_diff_sum += (
                nbr_particle.mass*
                self.cubic_spline_kernel_gradient(particle.initial_pos - nbr_particle.initial_pos)
            )
            sum_pos_diff_norm += (
                np.linalg.norm(
                    nbr_particle.mass * self.cubic_spline_kernel_gradient(
                        particle.initial_pos - nbr_particle.initial_pos
                    )
                )
            )
        particle.divergence_factor = (
            np.linalg.norm(pos_diff_sum) +
            sum_pos_diff_norm
        )

    def calculate_density_error(self):
        return self.particle.predicted_density - self.PARAMETERS["mass_density"]

    def find_neighbour_list(self, particle):
        particle.neighbour_list = []
        for nbr in self.hash_table[particle.hash_value]:
            if particle is not nbr:
                particle.neighbour_list.append(nbr)

    def recompute_neighbour_list(self, particle):
        if self.search_method != "Neighbour":
            particle.neighbour_list = []
            hash_value = SpatialHashing(self.PARAMETERS["cell_size"], self.num_particles).find_hash_value(particle)
            try:
                for nbr in self.hash_table[hash_value]:
                    if particle is not nbr:
                        particle.neighbour_list.append(nbr)
            except KeyError:
                pass

    # ----------------------------------------------------------------- DENSITY CORRECTION -----------------------------------------------------------------------------

    def update_stiffness_k(self, particle):
        
        try:
            time_integral = 1/ m.pow(self.delta_time, 2)
            mass_diff = particle.predicted_density - self.PARAMETERS["mass_density"]
            particle.stiffness_k = (
                time_integral *
                mass_diff *
                particle.divergence_factor
            )
        except ZeroDivisionError:
            particle.stiffness_k = np.array([0, 0, 0], dtype="float64")

    def update_predicted_density(self, particle):
        for nbr_particle in particle.neighbour_list:
            particle.predicted_density += (
                nbr_particle.mass * (particle.predicted_velocity - nbr_particle.predicted_velocity) *
                self.cubic_spline_kernel_gradient(particle.initial_pos - nbr_particle.initial_pos)
            )
        particle.predicted_density *= self.delta_time
        particle.predicted_density += particle.mass_density

    def update_predicted_density_velocity(self, particle):

        for nbr in particle.neighbour_list:
            self.update_stiffness_k(nbr)

        self.update_stiffness_k(particle)

        divergence_density = np.array([0, 0, 0], dtype="float64")
        for nbr_particle in self.neighbours_list:

            try:
                if nbr_particle.mass_density == 0:
                    nbr_stiffness_density = 0
                else:
                    nbr_stiffness_density = nbr_particle.stiffness_k / nbr_particle.mass_density
            except ZeroDivisionError:
                stiffness_density = np.array([0, 0, 0], dtype="float64")

            try:
                if self.particle.mass_density == 0:
                    stiffness_density = 0
                else:
                    stiffness_density = self.particle.stiffness_k / self.particle.mass_density
            except ZeroDivisionError:
                stiffness_density = np.array([0, 0, 0], dtype="float64")

            divergence_density += (
                nbr_particle.mass*
                (
                    stiffness_density +
                    nbr_stiffness_density
                ) *
                self.cubic_spline_kernel_gradient(
                    particle.initial_pos - nbr_particle.initial_pos
                )
            )
        particle.predicted_velocity -= (
            self.delta_time * divergence_density
        )

    def correct_density_error(self):

        iter_step = 0
        self.particle.predicted_density = self.particle.mass_density

        """ self.adapt_to_CFL() """

        if self.calculate_density_error() > 0.01*self.PARAMETERS["mass_density"] and iter_step < self.OTHER_PARAMS["max_iterations"]:

            for nbr in self.neighbours_list:
                self.update_predicted_density(nbr)
                self.update_predicted_density_velocity(nbr)

            self.update_predicted_density(self.particle)
            self.update_predicted_density_velocity(self.particle)

            iter_step += 1

        if iter_step > self.OTHER_PARAMS["max_iterations"]:
            print("Corrected density is:", self.particle.predicted_density)
            time.sleep(0.1)
    
    # --------------------------------------------------------------- DIVERGENCE CORRECTION ----------------------------------------------------------------------------

    def correct_divergence_error(self):

        iterations = 0
        while any(self.particle.divergence) > self.OTHER_PARAMS["divergence_error"] and \
            iterations < self.OTHER_PARAMS["max_iterations_div"]:

            print("Entering divergence correction")

            self.update_divergence_velocity()

            iterations +=1 

    def update_stiffness_k_v(self, particle):

        divergence_factor = 1 / particle.divergence_factor

        if particle.divergence_factor == 0:
            divergence_factor = 0

        particle.stiffness_k_v = (
            1 / self.delta_time * 
            particle.divergence * 
            1 / divergence_factor
        )

    def update_divergence(self, particle):
        particle.divergence = np.array([0, 0, 0], dtype="float64")
        for nbr_particle in particle.neighbour_list:
            particle.divergence += (
                nbr_particle.mass * (
                particle.predicted_velocity - nbr_particle.predicted_velocity
                ) *
                self.cubic_spline_kernel_gradient(
                    particle.initial_pos - nbr_particle.initial_pos
                )
            )

    def update_divergence_velocity(self):

        for nbr in self.particle.neighbour_list:
            self.update_divergence(nbr)
            self.update_stiffness_k_v(nbr)

        self.update_divergence(self.particle)
        self.update_stiffness_k_v(self.particle)

        divergence_velocity = np.array([0, 0, 0], dtype="float64")
        for nbr_particle in self.particle.neighbour_list:
            try:
                divergence_velocity += (
                    nbr_particle.mass *
                    (
                        (self.particle.stiffness_k_v / self.particle.mass_density) +
                        (nbr_particle.stiffness_k_v / nbr_particle.mass_density)
                    ) *
                    self.cubic_spline_kernel_gradient(
                        self.particle.initial_pos - nbr_particle.initial_pos
                    )
                )
            except (ValueError, TypeError):
                divergence_velocity += np.array([0, 0, 0], dtype="float64")
        self.particle.predicted_velocity -= (
            self.delta_time * divergence_velocity
        )

    # -------------------------------------------------------------------- UPDATES -------------------------------------------------------------------------------------

    def debugging_forces(self, secs):

        print("Mass Density:", self.particle.mass_density)
        print("Predicted Density", self.particle.predicted_density)
        print("Density error is:", self.calculate_density_error() )
        print("Divergence:", self.particle.divergence)
        print("Divergence factor:", self.particle.divergence_factor)
        print("Buoyancy:", self.particle.buoyancy)
        print("Gravity:", self.particle.gravity)
        print("viscosity:", self.particle.viscosity)
        print("Delta time is:", self.delta_time)
        print("\n\n")
        time.sleep(secs)

    def update_attrs(self):

        self.particle.initial_pos += self.delta_time*self.particle.predicted_velocity
            
        self.recompute_neighbour_list(self.particle)

        for nbr in self.particle.neighbour_list:
            self.update_predicted_mass_density(nbr)
            self.update_divergence_factor(nbr)

        self.update_predicted_mass_density(self.particle)
        self.update_divergence_factor(self.particle)

    def update_errors(self):

        self.correct_density_error()
        self.update_attrs()

        self.correct_divergence_error()

    def update(self):
        
        self.update_errors()
        
        self.particle.velocity = self.particle.predicted_velocity

        self.XSPH_vel_correction()
        self.choose_collision_types("Cuboid", "Normal")

        """ self.debugging_forces(0.1) """