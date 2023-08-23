import math as m
import numpy as np
import cupy as cp

from Particles.particles import Particle

class UtilityCalculations:

    def __init__(self, hash_table:dict = None,
                 parameters:dict = None):
        self.hash_table = hash_table
        self.PARAMETERS = parameters
    
    # --------------------------------------------------------------------- GENERAL METHODS ------------------------------------------------------------------------------------

    def kronecker_delta(self, i, j):
        if (i==j): return 1 
        else: return 0

    def reverse_kronecker_delta(self, i, j):
        if (i==j): return 0
        else: return 1

    def normalize(self, position):

        if np.linalg.norm(position)==0:
            position = np.array([0, 0, 0], dtype="float64")
            return position
        else:
            return position / np.linalg.norm(position)

    # ------------------------------------------------------------------- CALCULATION METHODS ----------------------------------------------------------------------------------

    def cubic_spline_kernel_grad_c(self, position, index_comp:int = 0):
        const_term = 15 / (np.pi*m.pow(self.PARAMETERS["cell_size"], 6))
        kernel_term = (position[index_comp] / np.linalg.norm(position)) * \
                    m.pow((self.PARAMETERS["cell_size"] - np.linalg.norm(position)), 2)
        return const_term*kernel_term
        
    def cubic_spline_kernel_pos(self, position):
        q = np.linalg.norm(position) / self.PARAMETERS["cell_size"]
        if q >= 0 and q <= 0.5:
            return 1 - 6*m.pow(q, 2) + 6*m.pow(q, 3)
        if q > 0.5 and  q<= 1:
            return 2*m.pow(1-q, 3)
        if q>1:
            return 0
        
    def cubic_spline_kernel_grad(self, position):
        const_term = 15 / (np.pi*m.pow(self.PARAMETERS["cell_size"], 6))
        kernel_term = (position / np.linalg.norm(position)) * \
                    m.pow((self.PARAMETERS["cell_size"] - np.linalg.norm(position)), 2)
        return const_term*kernel_term

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
        
    def cubic_spline_kernel_gradient_mag(self, position):
        return self.cubic_spline_kernel_gradient(np.abs(position))
    
    def new_cubic_spline_kernel_gradient(self, position:np.array=None):
        if np.linalg.norm(position) == 0: q = 0
        q = np.linalg.norm(position) / self.PARAMETERS["cell_size"]
        if q>=0 and q<1:
            kernel_val = m.pow(1 - q, 2) * self.normalize(position)
            return kernel_val
        if q>=1:
            return np.array([0, 0, 0], dtype="float64")
    
    # ---------------------------------------------------------------------- FORCE UPDATES ------------------------------------------------------------------------------------

    def update_predicted_mass_density(self, particle):
        density = 0 
        for nbr_particle in particle.neighbour_list:
            kernel_value = self.kernel_linear(particle.predicted_initial_pos - nbr_particle.predicted_initial_pos, 0)
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
        try:
            for nbr in self.hash_table[particle.hash_value]:
                if particle is not nbr:
                    particle.neighbour_list.append(nbr)
        except KeyError:
            pass
