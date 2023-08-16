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
        "max_iterations":2,
        "max_iterations_div":2,
        "divergence_error":0.0003
    }

    def __init__(self,
                particle: Particle=None,
                search_method:str = None,
                hash_table:dict = None,
                hash_value:int = None,
                time_stepping:str = "Euler Cromer",
                tank_attrs:dict = None,
                num_particles:int = None,
                delta_time:int = 0.02):
        
        super().__init__(particle=particle,
                         search_method=search_method,
                         hash_table=hash_table,
                         hash_value=hash_value,
                         time_stepping=time_stepping,
                         tank_attrs=tank_attrs,
                         delta_time=delta_time)
        
        self.divergence = np.array([0, 0, 0], dtype="float64")
        self.divergence_force = np.array([0, 0, 0], dtype="float64")
        self.num_particles = num_particles

        self.predict_advective_forces()
        self.update_divergence_factor(self.particle)
        
    # ------------------------------------------------------------------- PREDICT FORCES -------------------------------------------------------------------------------

    def update_mass_density(self, particle):
        """
        """
        density = 0 
        for nbr_particle in particle.neighbour_list:
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
        particle.gravity = self.gravity_const*particle.mass
    
    def update_normal_field(self, particle):
        """
        """
        normal_field = np.array([0, 0, 0],  dtype="float64")
        for id, nbr_particle in enumerate(particle.neighbour_list):
            pos_difference = particle.initial_pos - nbr_particle.initial_pos
            normal_field += (
                particle.mass * 1/particle.mass_density * self.kernel_gradient(pos_difference, 0)
            )
        return normal_field

    def update_surface_curvature(self, particle):
        """        
        """
        surface_curvature = 0
        for id, nbr_particle in enumerate(particle.neighbour_list):
            surface_curvature += (
                particle.mass * 1/particle.mass_density * self.kernel_laplacian(
                particle.initial_pos - nbr_particle.initial_pos, 0)
            )
        return surface_curvature
    
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

    # ----------------------------------------------------------------- UTILITY FUNCTIONS -----------------------------------------------------------------------------

    def predict_advective_forces(self):
        
        for particle in self.neighbours_list:
            self.find_neighbour_list(particle)
            self.update_mass_density(particle)
            self.update_viscosity(particle)
            self.update_gravity(particle)
            self.update_buoyancy(particle)
            self.update_surface_tension(particle)

            self.all_forces += (
                particle.viscosity +
                particle.gravity +
                particle.buoyancy +
                particle.surface_tension
            )

            particle.predicted_velocity = particle.velocity + self.delta_time*(self.all_forces/particle.mass)
        
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
    
    def cubic_spline_kernel_gradient(self, position):
        q = np.linalg.norm(position) / self.PARAMETERS["cell_size"]
        kernel_const = 1 / (np.pi*m.pow(self.PARAMETERS["cell_size"], 4))
        if q>=0 and q<=1:
            kernel_val = 9/4*m.pow(q, 2) - 3*q
            return kernel_val*kernel_const
        if q>=1 and q<=2:
            kernel_val = -3/4*m.pow((2 - q), 2)
            return kernel_val*kernel_const
        if q>=2:
            return 0
    
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
        particle.stiffness_k = (
            1/ m.pow(self.delta_time, 2) *
            (particle.predicted_density - self.PARAMETERS["mass_density"]) *
            particle.divergence_factor
        )
    
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
            divergence_density += (
                nbr_particle.mass*
                (
                    particle.stiffness_k / particle.mass_density +
                    nbr_particle.stiffness_k / nbr_particle.mass_density
                )*
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

        if self.calculate_density_error() > 0.01*self.PARAMETERS["mass_density"] and \
            iter_step < self.OTHER_PARAMS["max_iterations"]:

            """ print("Entering density correction") """

            for nbr in self.neighbours_list:
                self.update_predicted_density(nbr)
                self.update_predicted_density_velocity(nbr)

            self.update_predicted_density(self.particle)
            self.update_predicted_density_velocity(self.particle)

            """ self.debugging_forces(0.1) """

            iter_step += 1
    
    # --------------------------------------------------------------- DIVERGENCE CORRECTION ----------------------------------------------------------------------------

    def correct_divergence_error(self):

        iter_step = 0
        while self.particle.divergence > self.OTHER_PARAMS["divergence_error"] and \
            iter_step < self.OTHER_PARAMS["max_iterations_div"]:

            print("Entering divergence correction")

            self.update_divergence_velocity()
            self.debugging_forces(0.1)

            iter_step +=1 
        
        self.particle.velocity = self.particle.predicted_velocity

    def update_stiffness_k_v(self, particle):
        particle.stiffness_k_v = (
            1 / self.delta_time * 
            particle.divergence * 
            1 / particle.divergence_factor
        )

    def update_divergence(self, particle):
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
        print("\n\n")
        time.sleep(secs)

    def update_attrs(self):

        self.particle.initial_pos += self.delta_time*self.particle.predicted_velocity

        self.recompute_neighbour_list(self.particle)

        for nbr in self.particle.neighbour_list:
            self.update_mass_density(nbr)
            self.update_divergence_factor(nbr)

        self.update_mass_density(self.particle)
        self.update_divergence_factor(self.particle)

    def update(self):
        
        self.correct_density_error()
        self.update_attrs()

        self.correct_divergence_error()

        self.XSPH_vel_correction()
        self.choose_collision_types()
        self.choose_time_stepping() 