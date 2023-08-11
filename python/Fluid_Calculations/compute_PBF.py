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

    def update_constraint_function(self):
        self.particle.constraint_function = (
            self.particle.mass_density / self.PARAMETERS["mass_density"] - 1
        )

    def update_constraint_grad(self, id):

        self.constraint_grad = np.array([0, 0, 0], dtype="float64")
        if self.all_particles[id] == self.particle:
            for nbr_particle in self.neighbours_list:
                self.constraint_grad += (
                    self.cubic_spline_kernel_gradient(self.particle.initial_pos - 
                                                        nbr_particle.initial_pos)
                )
        else:
            for nbr_particle in self.all_particles[id].neighbour_list:
                self.constraint_grad += (
                    -self.cubic_spline_kernel_gradient(self.all_particles[id].initial_pos - 
                                                    nbr_particle.initial_pos)
                )
        self.constraint_grad *= 1/self.PARAMETERS["mass_density"]

    def update_constraint(self):
        constraint_sum = 0
        for id, p in enumerate(self.all_particles):
            constraint_sum += (
                m.pow(np.linalg.norm(self.update_constraint_grad(id)), 2)
            )
        self.particle.constraint = self.particle.constraint_function / (constraint_sum + self.OTHER_PARAMS["relaxation_factor"])

    def update_s_corr(self, id):
        delta_q_const = np.array([0.1*self.PARAMETERS["cell_size"],
                                  0.1*self.PARAMETERS["cell_size"],
                                  0.1*self.PARAMETERS["cell_size"]], dtype="float64")
        s_corr = m.pow((
            -self.OTHER_PARAMS["k_const"]*
            (
                self.cubic_spline_kernel_pos(self.particle.initial_pos - self.neighbours_list[id].initial_pos) / 
                self.cubic_spline_kernel_pos(delta_q_const)
            )
        ), self.OTHER_PARAMS["n_const"])
        return s_corr
    
    def update_del_position(self):
        constraint_term = np.array([0, 0, 0], dtype="float64")
        for id, nbr_particle in enumerate(self.neighbours_list):
            constraint_term += (
                self.particle.constraint + nbr_particle.constraint + self.update_s_corr(id)*
                self.cubic_spline_kernel_gradient(self.particle.initial_pos - nbr_particle.initial_pos)
            )
        self.particle.del_position = 1/self.PARAMETERS["mass_density"]*constraint_term

    def update_position(self):
        self.predicted_position += self.particle.del_position

    def update_vorticity(self):
        for nbr_particle in self.neighbour_list:
            self.particle.vorticity += (
                np.cross((nbr_particle.velocity - self.particle.velocity),
                self.update_constraint_grad(self.all_particles.index(nbr_particle)))
            )

    def n_const(self):
        return self.cubic_spline_kernel_gradient(np.linalg.norm(self.particle.vorticity))
    
    def update_vorticity_force(self):
        self.particle.vorticity_force = (
            self.OTHER_PARAMS["relaxation_factor"]*
            np.cross(self.n_const()/ np.linalg.norm(self.n_const()), self.particle.vorticity)
        )

    def update_all_forces(self):
        return super().update_all_forces()
    
    def update_velocity(self):
        self.particle.velocity = 1 / self.delta_time * (self.predicted_position - self.particle.initial_pos)
    
    def update(self):
        
        self.update_all_forces()
        self.particle.acceleration = self.all_forces / self.particle.mass

        self.predicted_position = self.prediction_update(self.time_stepping, self.particle)[0]

        while(self.iter < self.OTHER_PARAMS["solver_iterations"]):

            self.update_constraint_function()
            self.update_constraint()
            self.update_del_position()
            self.update_position()

            self.choose_collision_types("Cuboid", "Normal")

        self.update_velocity()

        self.update_vorticity()
        self.update_vorticity_force()

        self.all_forces += self.particle.vorticity_force
        self.particle.acceleration = self.all_forces / self.particle.mass

        self.choose_time_stepping(self.time_stepping)
        self.XSPH_vel_correction()
        