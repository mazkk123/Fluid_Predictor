import math as m
import numpy as np
import random as rd
import re

from Fluid_Calculations.compute_sph import SPH
from Particles.particles import Particle

class PBF(SPH):
    
    OTHER_PARAMS = {
        "relaxation_factor":0.7,
        "k_const":0.1,
        "n_const":4,
        "solver_iterations":2,
        "c_const":0.01
    }

    def __init__(self,
                particle: Particle=None,
                search_method:str = None,
                hash_table:dict = None,
                hash_value:int = None,
                time_stepping:str = "Euler Cromer",
                tank_attrs:dict = None,
                all_particles:list = None,
                delta_time:int = 0.02):
        
        super().__init__(particle=particle,
                         search_method=search_method,
                         hash_table=hash_table,
                         hash_value=hash_value,
                         time_stepping=time_stepping,
                         tank_attrs=tank_attrs,
                         delta_time=delta_time)
        
        self.all_particles = all_particles
        self.constraint_grad = np.array([0, 0, 0], dtype="float64")
        self.iter = 0
        
        self.predict_positions(self.particle, 4)
    
    def predict_positions(self, particle, depth:int=4):
        
        if depth==0:
            return
            
        self.find_neighbour_list(particle)
        self.update_mass_density(particle)
        self.update_gravity(particle)
        self.update_buoyancy(particle)
        self.update_viscosity(particle)
        self.update_surface_tension(particle)
        
        self.all_forces += (
            particle.gravity +
            particle.buoyancy +
            particle.viscosity +
            particle.surface_tension
        )
        
        particle.predicted_velocity = particle.velocity + self.delta_time* \
                                        self.all_forces
        particle.predicted_initial_pos = particle.initial_pos + self.delta_time* \
                                        particle.predicted_velocity
        for nbr in particle.neighbour_list:
            return self.predict_positions(nbr, depth-1)
                
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
        particle.gravity = self.gravity_const
    
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
        
    def update_constraint_function(self, particle, depth:int=4):
        
        if depth==0:
            return 
        
        self.find_neighbour_list(particle)
        particle.constraint_function = (
            particle.mass_density / self.PARAMETERS["mass_density"] - 1
        )
        
        for nbr in particle.neighbour_list:
            return self.update_constraint_function(nbr, depth-1)

    def update_constraint_grad(self, particle):

        self.constraint_grad = np.array([0, 0, 0], dtype="float64")            
        for nbr_particle in particle.neighbour_list:
            self.constraint_grad += (
                self.cubic_spline_kernel_gradient(particle.predicted_initial_pos - 
                                            nbr_particle.predicted_initial_pos)
                )
        self.constraint_grad *= 1/self.PARAMETERS["mass_density"]
        return self.constraint_grad

    def update_constraint(self, particle, depth:int=4):
        
        if depth==0:
            return 
            
        constraint_sum = 0
        for p in particle.neighbour_list:
            constraint_sum += (
                m.pow(np.linalg.norm(self.update_constraint_grad(p)), 2)
            )
        particle.constraint = particle.constraint_function / (constraint_sum + self.OTHER_PARAMS["relaxation_factor"])
        
        for nbr in particle.neighbour_list:
            return self.update_constraint(nbr, depth-1)
        
    def update_s_corr(self, id):
        delta_q_const = np.array([0.1*self.PARAMETERS["cell_size"],
                                  0.1*self.PARAMETERS["cell_size"],
                                  0.1*self.PARAMETERS["cell_size"]], dtype="float64")
        s_corr = 0
        s_corr = m.pow((
            -self.OTHER_PARAMS["k_const"]*
            (
                self.cubic_spline_kernel_pos(particle.predicted_initial_pos - 
                                             particle.neighbour_list[id].predicted_initial_pos) / 
                self.cubic_spline_kernel_pos(delta_q_const)
            )
        ), self.OTHER_PARAMS["n_const"])
        return s_corr
    
    def update_del_position(self, particle, depth:int=4):
        
        if depth==0:
            return 
            
        constraint_term = np.array([0, 0, 0], dtype="float64")
        for id, nbr_particle in enumerate(particle.neighbour_list):
            constraint_term += (
                self.particle.constraint + nbr_particle.constraint + self.update_s_corr(id)*
                self.cubic_spline_kernel_gradient(particle.predicted_initial_pos - nbr_particle.predicted_initial_pos)
            )
        particle.del_position = 1/self.PARAMETERS["mass_density"]*constraint_term
        
        for nbr in particle.neighbour_list:
            self.update_del_position(nbr, depth-1)

    def update_position(self, particle, depth):
        
        if depth==0:
            return 
            
        particle.predicted_initial_pos += particle.del_position
        
        for nbr in particle.neighbour_list:
            return self.update_position(nbr, depth-1)
            
    def update_vorticity(self):
        for nbr_particle in particle.neighbour_list:
            particle.vorticity += (
                np.cross((nbr_particle.velocity - particle.velocity),
                          self.update_constraint_grad(nbr_particle))
            )

    def n_const(self, particle):
        return self.cubic_spline_kernel_gradient(np.linalg.norm(particle.vorticity))
    
    def update_vorticity_force(self):
        particle.vorticity_force = (
            self.OTHER_PARAMS["relaxation_factor"]*
            np.cross(self.n_const(particle)/ np.linalg.norm(self.n_const()), particle.vorticity)
        )
        
    def update_vorticity_forces(self, particle, depth:int=4):
        
        if depth==0:
            return 
            
        self.update_vorticity(particle)
        self.update_vorticity_force(particle)
        for nbr in particle.neighbour_list:
            return self.update_vorticity_forces(nbr, depth-1)
    
    def update_velocity(self, particle, depth:int=4):
        
        if depth==0:
            return
            
        particle.predicted_velocity = 1 / self.delta_time * (particle.predicted_initial_pos - particle.initial_pos)
        
        for nbr in particle.neigbour_list:
            return self.update_velocity(nbr, depth-1)
    
    def XSPH_viscosity(self):
        correction_term = np.array([0, 0, 0], dtype="float64")
        for nbr in self.particle.neighbour_list:
            correction_term += (
                (self.particle.velocity - nbr.velocity) *
                self.cubic_spline_kernel_pos(self.particle.initial_pos - nbr.initial_pos)
            )
        correction_term *= self.OTHER_PARAMS["c_const"]
        self.particle.velocity += correction_term
        
    def update_forces(self, particle):
        
        self.all_forces = np.array([0, 0, 0], dtype="float64")
        self.all_forces += (
            particle.buoyancy +
            particle.viscosity +
            particle.surface_tension +
            particle.gravity +
            particle.vorticity_force
        )
        particle.acceleration = self.all_forces / particle.mass
        particle.velocity = particle.predicted_velocity + self.delta_time*particle.acceleration
        particle.initial_pos = particle.predicted_initial_pos + self.delta_time*particle.velocity
        
    def update(self):

        while(self.iter < self.OTHER_PARAMS["solver_iterations"]):

            self.update_constraint_function(self.particle, 4)
            self.update_constraint(self.particle, 4)
            self.update_del_position(self.particle, 4)
            self.update_position(self.particle, 4)
            
            self.iter += 1

        self.update_velocity(self.particle, 4)
        self.update_vorticity_forces(self.particle, 4)
        
        self.update_forces(self.particle)
        self.XSPH_vel_correction()
        self.XSPH_viscosity_correction()
        