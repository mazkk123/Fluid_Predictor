import math as m
import numpy as np
import random as rd
import re

from search_methods import SpatialHashing, CompactHashing, ZSorting, NearestNeighbour
from Collisions.box_collisions import BoxCollisions, AABB, OrientedBBox
from Collisions.sphere_collisions import SphereCollisions, CapsuleCollisions
from Particles.particles import Particle

class SPH(Particle):

    PARAMETERS = {
        "grid_separation":0.1,
        "cell_size":0.15,
        "mass": 0.1,
        "viscosity": 3.5,
        "mass_density": 998.2,
        "buoyancy":0.0,
        "tension_coefficient":0.0725,
        "tension_threshold":7,
        "pressure_const":7,
        "loss_of_speed":0.5,
        "epsilon":0.1
    }

    def __init__(self,
                 particle: Particle=None,
                 search_method: str=None,
                 hash_table: dict=None,
                 hash_value: int=None,
                 params: dict=None,
                 delta_time: float=0.02):

        if params is not None:
            self.params = params
        if particle is not None:
            self.particle = particle
        if search_method is not None:
            self.search_method = search_method
        if hash_value is not None:
            self.hash_value = hash_value
        if hash_table is not None:
            self.hash_table = hash_table
        if delta_time is not None:
            self.delta_time = delta_time

        self.neighbours_list = []
        self.update_particle_neighbours(self)
        self.gravity_const = np.array([0, -9.81, 0])

    def update_particle_neighbours(self):
        """
            find all particles hashed to the same cell in the
            global HASH_MAP param
        """
        if (len(self.hash_table[self.hash_value]!=0)):
            for items in self.hash_table[self.hash_value]:
                self.neighbours_list.append(items)
        
    def kernel_gradient(self, position: np.array=None, kernel_type:int = 0):
        """
        """
        if kernel_type==0:
            kernel_value = np.power(position * (m.pow(self.PARAMETERS["cell_size"], 2) - np.power(np.linalg.norm(position),2)), 2)
            kernel_const = -945/(32 * np.pi * m.pow(self.PARAMETERS["cell_size"], 9))
            return kernel_value*kernel_const
        if kernel_type==1:
            kernel_value = position/np.linalg.norm(position)*np.power((
                self.PARAMETERS["cell_size"] - np.linalg.norm(position)
                ), 2)
            kernel_const = -45/(np.pi*m.pow(self.PARAMETERS["cell_size"], 6))
            return kernel_value*kernel_const
        if kernel_type==2:
            pass

    def kernel_laplacian(self, position: np.array=None, kernel_type:int = 0):
        """
        """
        if kernel_type==0:
            kernel_val =( (m.pow(self.PARAMETERS["cell_size"], 2) - np.linalg.norm(position)) *
                          (3 * m.pow(self.PARAMETERS["cell_size"], 2) - 7*np.power(np.linalg.norm(position), 2))
                         )
            kernel_const = -945/(32 * np.pi * m.pow(self.PARAMETERS["cell_size"], 9))
            return kernel_val*kernel_const
        if kernel_type==1:
            kernel_val = (self.PARAMETERS["cell_size"] - np.linalg.norm(position)) * (
                self.PARAMETERS["cell_size"] - 2*np.linalg.norm(position)
            )
            kernel_const = -90/(np.pi*m.pow(self.PARAMETERS["cell_size"], 6))
            return kernel_val*kernel_const
        if kernel_type==2:
            pass

    def kernel_linear(self, position: np.array=None, kernel_type:int = 0):
        """
        """
        if kernel_type==0:
            kernel_val = np.power(m.pow(self.PARAMETERS["cell_size"], 2) - np.power(np.linalg.norm(position), 2), 3)
            kernel_const = 315/(64*np.pi*m.pow(self.PARAMETERS["cell_size"], 9))
            if np.linalg.norm(position) >= 0 and np.linalg.norm<=self.PARAMETERS["cell_size"]:
                return kernel_val*kernel_const
            else:
                return 0
        if kernel_type==1:
            kernel_val = np.power((self.PARAMETERS["cell_size"] - np.linalg.norm(position)), 3)
            kernel_const = 15/(np.pi*m.pow(self.PARAMETERS["cell_size"], 6))
            if np.linalg.norm(position) >= 0 and np.linalg.norm(position) <= self.PARAMETERS["cell_size"]:
                return kernel_val*kernel_const
            else:
                return 0
        if kernel_type==2:
            pass

    def update_mass_density(self):
        """
        """
        density = 0 
        for id, nbr_particle in enumerate(self.neighbours_list):
            kernel_value = self.kernel_linear(self.particle.initial_pos - nbr_particle.initial_pos, 0)
            density += kernel_value*self.particle.mass

        self.particle.mass_density = self.PARAMETERS["mass_density"] + density


    def update_pressure(self):
        """
        """
        pressure = self.PARAMETERS["pressure_const"]*(
            self.particle.mass_density - self.PARAMETERS["mass_density"] 
        )
        self.particle.pressure = pressure

    def update_pressure_force(self):
        """
        """
        pressure_force = np.array([0, 0, 0])
        for id, nbr_particle in enumerate(self.neighbours_list):
            pressure_avg = (self.particle.pressure + nbr_particle.pressure) / 2
            density_avg = self.particle.mass / nbr_particle.mass_density
            kernel_grad = self.kernel_gradient(
                self.particle.initial_pos - nbr_particle.initial_pos
            )
            pressure_force += kernel_grad*pressure_avg*density_avg
        self.particle.pressure_force = -1*pressure_force

    def update_viscosity(self):
        """
        """
        viscosity = np.array([0, 0, 0])
        for id, nbr_particle in enumerate(self.neighbours_list):
            vel_dif = nbr_particle.velocity - self.particle.velocity
            kernel_laplacian = (
                self.particle.initial_pos - nbr_particle.initial_pos
            )
            mass_pressure = self.particle.mass/nbr_particle.mass_density
            viscosity += vel_dif*mass_pressure*kernel_laplacian
        
        self.viscosity = viscosity*self.PARAMETERS["viscosity"]

    def update_gravity(self):
        """
        """
        self.particle.gravity = self.particle.mass_density * self.gravity_const

    def update_color_gradient(self):
        """
        """
        colour_field = np.array([0, 0, 0])
        for id, nbr_particle in enumerate(self.neighbours_list):
            pos_difference = self.particle.initial_pos - nbr_particle.initial_pos
            colour_field += (
                self.particle.mass * 1/self.particle.mass_density * self.kernel_linear(pos_difference, 0)
            )
        return colour_field

    def update_normal_field(self):
        """
        """
        normal_field = np.array([0, 0, 0])
        for id, nbr_particle in enumerate(self.neighbours_list):
            pos_difference = self.particle.initial_pos - nbr_particle.initial_pos
            normal_field += (
                self.particle.mass * 1/self.particle.mass_density * self.kernel_gradient(pos_difference, 0)
            )
        return normal_field

    def update_surface_curvature(self):
        """        
        """
        surface_curvature = np.array([0, 0, 0])
        for id, nbr_particle in enumerate(self.neighbours_list):
            pos_difference = self.particle.initial_pos - nbr_particle.initial_pos
            surface_curvature += (
                self.particle.mass * 1/self.particle.mass_density * self.kernel_laplacian(pos_difference, 0)
            )
        return surface_curvature

    def update_surface_tension(self):
        """
        """
        normal_field = self.update_normal_field()
        surface_curvature = self.update_surface_curvature()
        normal_field_magnitude = np.linalg.norm(normal_field)

        self.particle.surface_tension = (
            self.PARAMETERS["tension_coefficient"] * surface_curvature * normal_field/normal_field_magnitude 
        )

    def update_buoyancy(self):
        """
        """
        buoyancy = self.PARAMETERS["buoyancy"] * (self.particle.mass_density - self.PARAMETERS["mass_density"])
        buoyancy *= self.PARAMETERS["gravity"]
        self.buoyancy = buoyancy

    def XSPH_vel_correction(self):
        """
        """
        new_vel = np.array([0, 0, 0])
        for id, nbr_particle in enumerate(self.neighbours_list):
            average_density = self.particle.mass_density + nbr_particle.density
            new_vel += (
                2*nbr_particle.mass/(average_density) * self.kernel_linear(self.particle.initial_pos - nbr_particle.initial_pos)
            )
        self.particle.velocity += self.PARAMETERS["epsilon"] * new_vel

    def update_all_forces(self):
        """
        """
        self.update_mass_density()
        self.update_pressure()
        self.update_pressure_force()
        self.update_gravity()
        self.update_buoyancy()
        self.update_surface_tension()

        self.all_forces = self.particle.mass_density + \
                        self.particle.pressure + self.particle.pressure_force + \
                        self.gravity + self.buoyancy + \
                        self.particle.surface_tension
        
    def update(self):

        self.update_all_forces()
        self.XSPH_vel_correction()

        self.particle.acceleration = self.all_forces / self.PARAMETERS["mass_density"]
        self.particle.velocity += self.delta_time*self.particle.acceleration
        self.particle.initial_pos += self.particle.velocity*self.delta_time

