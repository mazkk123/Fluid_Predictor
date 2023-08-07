import math as m
import numpy as np
import random as rd
import re

from Fluid_Utilities.search_methods import NearestNeighbour
from Fluid_Utilities.time_stepping import ForwardEuler, EulerCromer, LeapFrog, Verlet, IndependentTime, \
                                    RegionalShortTime, RungeKutta
from Collisions.box_collisions import BoxCollisions, AABB, OrientedBBox
from Collisions.sphere_collisions import SphereCollisions, CapsuleCollisions, CylinderCollisions
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
        "epsilon":0.1,
        "neighbour_num":30
    }

    TIME_SCHEMES = {
        "Forward Euler": 0,
        "Euler Cromer": 1,
        "Leap Frog": 2,
        "Verlet": 3,
        "Independent Time": 4,
        "Regional Short Time": 5,
        "Runge Kutta": 6
    }

    COLLISION_TYPES = {
        "Cuboid":{"OrientedBBox":0, "AABB":1, "Normal":2},
        "Cylinder":1,
        "Sphere":2,
        "Capsule":3,
        "Abstract":4
    }

    def __init__(self,
                 particle: Particle=None,
                 search_method: str="Spatial Hashing",
                 hash_table: dict=None,
                 hash_value: int=None,
                 delta_time: float=0.02,
                 time_stepping: str="Euler Cromer",
                 tank_attrs: dict=None,
                 collision_type: str="box"):
        
        if tank_attrs is not None:
            self.tank_attrs = tank_attrs
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
        if time_stepping is not None:
            self.time_stepping = time_stepping
        if collision_type is not None:
            self.collision_type = collision_type
        
        self.neighbours_list = []
        self.update_particle_neighbours(self)
        self.gravity_const = np.array([0, -9.81, 0])

    def update_particle_neighbours(self):
        """
            find all particles hashed to the same cell in the
            global HASH_MAP param
        """
        if self.search_method != "Neighbour":
            if (len(self.hash_table[self.hash_value]!=0)):
                for items in self.hash_table[self.hash_value]:
                    self.neighbours_list.append(items)
        else:
            for items in NearestNeighbour(search_radius=self.PARAMETERS["cell_size"], particle=self.particle,
                                          neighbour_size=self.PARAMETERS["neighbour_num"]).find_neighbours():
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
            density += kernel_value*nbr_particle.mass

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
    
    def choose_time_stepping(self, time_step_type:str = "Euler Cromer"):
        """
        
        """
        if self.TIME_SCHEMES[time_step_type] == 0:
            ForwardEuler(
                self.particle,
                self.delta_time
            )
        if self.TIME_SCHEMES[time_step_type] == 1:
            EulerCromer(
                self.particle,
                self.delta_time
            ).exec_time_scheme(self.delta_time)
        if self.TIME_SCHEMES[time_step_type] == 2:
            LeapFrog(
                self.particle,
                self.delta_time
            ).exec_time_scheme(self.delta_time)
        if self.TIME_SCHEMES[time_step_type] == 3:
            Verlet(
                self.particle,
                self.delta_time
            ).exec_time_scheme(self.delta_time)
        if self.TIME_SCHEMES[time_step_type] == 4:
            IndependentTime().exec_time_scheme()
        if self.TIME_SCHEMES[time_step_type] == 5:
            RegionalShortTime().exec_time_scheme()
        if self.TIME_SCHEMES[time_step_type] == 6:
            RungeKutta()
        else:
            EulerCromer(
                self.particle,
                self.delta_time
            ).exec_time_scheme(self.delta_time)

    def density_prediction(self):

        density = 0 
        for id, nbr_particle in enumerate(self.neighbours_list):
            kernel_value = self.kernel_linear(self.particle.initial_pos - nbr_particle.initial_pos, 0)
            density += kernel_value*self.particle.mass

        return self.PARAMETERS["mass_density"] + density
    
    def prediction_update(self, time_step_type:str = "Euler Cromer",
                             particle: Particle=None):
        if particle is not None:
            if self.TIME_SCHEMES[time_step_type] == 0:
                ForwardEuler(
                    particle,
                    self.delta_time
                )
            if self.TIME_SCHEMES[time_step_type] == 1:
                return EulerCromer(
                    particle,
                    self.delta_time
                ).get_time_scheme_values()
            if self.TIME_SCHEMES[time_step_type] == 2:
                return LeapFrog(
                    particle,
                    self.delta_time
                ).get_time_scheme_values()
            if self.TIME_SCHEMES[time_step_type] == 3:
                return Verlet(
                    particle,
                    self.delta_time
                ).get_time_scheme_values()
            if self.TIME_SCHEMES[time_step_type] == 4:
                IndependentTime()
            if self.TIME_SCHEMES[time_step_type] == 5:
                RegionalShortTime()
            if self.TIME_SCHEMES[time_step_type] == 6:
                RungeKutta()
            else:
                return EulerCromer(
                    particle,
                    self.delta_time
                ).get_time_scheme_values()

    def choose_collision_types(self, collision_type:str = "Cuboid",
                               secondary_type:str = "Normal",
                               particle:Particle = None):
        if particle is not None:
            if isinstance(self.COLLISION_TYPES[collision_type], dict):
                if self.COLLISION_TYPES[collision_type][secondary_type] == 0:
                    OrientedBBox()
                if self.COLLISION_TYPES[collision_type][secondary_type] == 1:
                    AABB()
                if self.COLLISION_TYPES[collision_type][secondary_type] == 2:
                    BoxCollisions(
                        particle=self.particle,
                        tank_size=self.tank_attrs["dimensions"]["size"],
                        tank_location=self.tank_attrs["diemensions"]["location"],
                        speed_loss=self.PARAMETERS["loss_of_speed"]
                    ).collision_resolution()
            else:
                if self.TIME_SCHEMES[collision_type] == 1:
                    CylinderCollisions()
                if self.TIME_SCHEMES[collision_type] == 2:
                    SphereCollisions()
                if self.TIME_SCHEMES[collision_type] == 3:
                    CapsuleCollisions()
                if self.TIME_SCHEMES[collision_type] == 4:
                    pass
                else:
                    BoxCollisions()

    def update(self):

        self.update_all_forces()
        self.XSPH_vel_correction()

        self.particle.acceleration = self.all_forces / self.PARAMETERS["mass_density"]
        self.particle.next_acceleration = self.particle.acceleration

        self.choose_time_stepping(self.time_stepping, self.particle)
        self.choose_collision_types("Cuboid", "Normal")

