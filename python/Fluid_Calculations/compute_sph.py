import math as m
import numpy as np
import random as rd
import re
import sys

sys.path.append("C:\\Users\\Student\\OneDrive - Bournemouth University\\Desktop\\Personal\\Python\\Fluid_Predictor\\python\\Fluid_Utilities\\")

from Fluid_Utilities.search_methods import NearestNeighbour
from Fluid_Utilities.time_stepping import ForwardEuler, EulerCromer, LeapFrog, Verlet, IndependentTime, \
                                    RegionalShortTime, RungeKutta
from Collisions.box_collisions import BoxCollisions, AABB, OrientedBBox
from Collisions.sphere_collisions import SphereCollisions, CapsuleCollisions, CylinderCollisions
from Particles.particles import Particle

class SPH(Particle):

    PARAMETERS = {
        "grid_separation":0.1,
        "cell_size":0.3,
        "mass": 0.1,
        "viscosity": 3.5,
        "mass_density": 998.2,
        "buoyancy":0.0,
        "tension_coefficient":0.0728,
        "tension_threshold":6,
        "pressure_const":7,
        "loss_of_speed":0.5,
        "epsilon":0.1,
        "neighbour_num":30,
        "thermal_exp_coeff":0.5,
        "kinematic_visc":2.5,
        "abs_temperature":20,
        "sound_speed":300
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
                 temperature:bool=False,
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
        
        self.temperature = temperature

        self.neighbours_list = []
        self.update_particle_neighbours()
        self.gravity_const = np.array([0, -9.81, 0], dtype="float64")
        self.num_test_steps = 50

    # ------------------------------------------------------------------ TRADITIONAL SPH ------------------------------------------------------------------------

    def update_particle_neighbours(self):
        """
            find all particles hashed to the same cell in the
            global HASH_MAP param
        """
        if self.search_method != "Neighbour":
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
            if np.linalg.norm(position) <= self.PARAMETERS["cell_size"] and \
                np.linalg.norm(position) >= 0:
                if np.linalg.norm(position) == 0:
                    kernel_value = np.array([0, 0, 0], dtype="float64")
                else:
                    kernel_value = np.power(position * (m.pow(self.PARAMETERS["cell_size"], 2) - np.power(np.linalg.norm(position),2)), 2)
                kernel_const = -945.0/(32 * np.pi * m.pow(self.PARAMETERS["cell_size"], 9))
                return kernel_value*kernel_const
            else:
                return np.array([0, 0, 0], dtype="float64")
        if kernel_type==1:
            if np.linalg.norm(position) <= self.PARAMETERS["cell_size"] and \
                np.linalg.norm(position) >= 0:
                if np.linalg.norm(position) == 0:
                    kernel_value =  np.array([0, 0, 0], dtype="float64")
                else:
                    kernel_value = position/np.linalg.norm(position)*(np.power(
                                    self.PARAMETERS["cell_size"] - np.linalg.norm(position), 2))           
                kernel_const = -45.0/(np.pi*m.pow(self.PARAMETERS["cell_size"], 6))
                return kernel_value*kernel_const
            else:
                return np.array([0, 0, 0], dtype="float64")
        if kernel_type==2:
            if np.linalg.norm(position) <= self.PARAMETERS["cell_size"] and \
                np.linalg.norm(position) >= 0:
                if np.linalg.norm(position) == 0:
                    kernel_value =  np.array([0, 0, 0], dtype="float64")
                else:
                    kernel_value = (
                        position * (
                        (-3 * np.linalg.norm(position) / 2 * m.pow(self.PARAMETERS["cell_size"], 3)) + 
                        (2 / m.pow(self.PARAMETERS["cell_size"], 2)) - 
                        (self.PARAMETERS["cell_size"] / 2 * m.pow(np.linalg.norm(position), 3))
                        )
                    )           
                kernel_const = 15.0/(2*np.pi*m.pow(self.PARAMETERS["cell_size"], 3))
                return kernel_value*kernel_const
            else:
                return np.array([0, 0, 0], dtype="float64")
            
    def kernel_laplacian(self, position: np.array=None, kernel_type:int = 0):
        """
        """
        if kernel_type==0:
            if np.linalg.norm(position) <= self.PARAMETERS["cell_size"] and \
                np.linalg.norm(position) >= 0:
                if np.linalg.norm(position) == 0:
                    kernel_val = 0
                else:
                    kernel_val =( (m.pow(self.PARAMETERS["cell_size"], 2) - np.linalg.norm(position)) *
                                (3 * m.pow(self.PARAMETERS["cell_size"], 2) - 7*np.power(np.linalg.norm(position), 2))
                                )
                kernel_const = -945.0/(32 * np.pi * m.pow(self.PARAMETERS["cell_size"], 9))
                return kernel_val*kernel_const
            else:
                return 0
        if kernel_type==1:
            if np.linalg.norm(position) <= self.PARAMETERS["cell_size"] and \
                np.linalg.norm(position) >= 0:
                if np.linalg.norm(position) == 0:
                    kernel_val = 0
                else:
                    kernel_val = (self.PARAMETERS["cell_size"] - np.linalg.norm(position)) * \
                        (self.PARAMETERS["cell_size"] - 2*np.linalg.norm(position))
                kernel_const = -90 / (np.pi * m.pow(self.PARAMETERS["cell_size"], 6))
                return kernel_val*kernel_const
            else:
                return 0
        if kernel_type==2:
            if np.linalg.norm(position) <= self.PARAMETERS["cell_size"] and \
                np.linalg.norm(position) >= 0:
                if np.linalg.norm(position) == 0:
                    kernel_val = 0
                else:
                    kernel_val = self.PARAMETERS["cell_size"] - np.linalg.norm(position)

                kernel_const = 45.0/(np.pi*m.pow(self.PARAMETERS["cell_size"], 6))

                return kernel_val*kernel_const
            else:
                return 0

    def kernel_linear(self, position: np.array=None, kernel_type:int = 0):
        """
        """
        if kernel_type==0:
            kernel_val = np.power(m.pow(self.PARAMETERS["cell_size"], 2) - np.power(np.linalg.norm(position), 2), 3)
            kernel_const = 315/(64*np.pi*m.pow(self.PARAMETERS["cell_size"], 9))
            if np.linalg.norm(position) >= 0 and np.linalg.norm(position)<=self.PARAMETERS["cell_size"]:
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
            kernel_val = (m.pow(self.PARAMETERS["cell_size"], 2) - m.pow(np.linalg.norm(position), 2)) * \
                        (3 * m.pow(self.PARAMETERS["cell_size"], 2) - 7 * m.pow(np.linalg.norm(position), 2))
            kernel_const = -945/(32 * np.pi*m.pow(self.PARAMETERS["cell_size"], 9))
            if np.linalg.norm(position) >= 0 and np.linalg.norm(position) <= self.PARAMETERS["cell_size"]:
                return kernel_val*kernel_const
            else:
                return 0

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
        pressure_force = np.array([0, 0, 0], dtype="float64")
        for id, nbr_particle in enumerate(self.neighbours_list):
            pressure_avg = (self.particle.pressure + nbr_particle.pressure) / 2
            try:
                density_avg = self.particle.mass / nbr_particle.mass_density
            except ZeroDivisionError:
                density_avg = 0
            kernel_grad = self.kernel_gradient(
                self.particle.initial_pos - nbr_particle.initial_pos, 1
            )
            pressure_force += kernel_grad*pressure_avg*density_avg
        self.particle.pressure_force = -1*pressure_force

    def update_viscosity(self):
        """
        """
        viscosity = np.array([0, 0, 0], dtype="float64")
        for id, nbr_particle in enumerate(self.neighbours_list):
            vel_dif = nbr_particle.velocity - self.particle.velocity
            kernel_laplacian = self.kernel_laplacian(self.particle.initial_pos - nbr_particle.initial_pos, 2)
            try:
                mass_pressure = self.particle.mass/nbr_particle.mass_density
            except ZeroDivisionError:
                mass_pressure = 0
            viscosity += vel_dif*mass_pressure*kernel_laplacian

        self.particle.viscosity = viscosity*self.PARAMETERS["viscosity"]

    def update_gravity(self):
        """
        """
        self.particle.gravity = self.gravity_const*self.particle.mass

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
        normal_field = np.array([0, 0, 0],  dtype="float64")
        for id, nbr_particle in enumerate(self.neighbours_list):
            pos_difference = self.particle.initial_pos - nbr_particle.initial_pos
            normal_field += (
                self.particle.mass * 1/self.particle.mass_density * self.kernel_gradient(pos_difference, 0)
            )
        return normal_field

    def update_surface_curvature(self):
        """        
        """
        surface_curvature = 0
        for id, nbr_particle in enumerate(self.neighbours_list):
            surface_curvature += (
                self.particle.mass * 1/self.particle.mass_density * self.kernel_laplacian(
                self.particle.initial_pos - nbr_particle.initial_pos, 0)
            )
        return surface_curvature

    def update_surface_tension(self):
        """
        """
        normal_field = self.update_normal_field()
        surface_curvature = self.update_surface_curvature()
        normal_field_magnitude = np.linalg.norm(normal_field)
        
        if normal_field_magnitude >= self.PARAMETERS["tension_threshold"]:
            self.particle.surface_tension = (
                self.PARAMETERS["tension_coefficient"] * surface_curvature * normal_field/normal_field_magnitude 
            )

    def update_buoyancy(self):
        """
        """
        buoyancy = self.PARAMETERS["buoyancy"] * (self.particle.mass_density - self.PARAMETERS["mass_density"])
        buoyancy *= self.gravity_const
        self.buoyancy = buoyancy

    def XSPH_vel_correction(self):
        """
        """
        new_vel = np.array([0, 0, 0], dtype="float64")
        for id, nbr_particle in enumerate(self.neighbours_list):
            average_density = self.particle.mass_density + nbr_particle.mass_density
            try:
                new_vel += (
                    2*nbr_particle.mass/(average_density) * self.kernel_linear(self.particle.initial_pos - nbr_particle.initial_pos, 0)
                )
            except ZeroDivisionError:
                new_vel = np.array([0, 0, 0], dtype="float64")

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
        self.update_viscosity() 

        if self.temperature is True:
            self.update_EOS_pressure()
            self.update_body_force()
            self.update_laminar_viscosity()
            self.update_heat_conduction()

            self.all_forces = self.particle.pressure_force + \
                              self.gravity_const + self.particle.buoyancy + \
                              self.particle.surface_tension + \
                              self.particle.body_force + self.particle.laminar_viscosity + \
                              self.particle.thermal_diffusion
        
        #self.debugging_forces()

        self.all_forces = self.particle.pressure_force + \
                          self.gravity_const + self.buoyancy + \
                          self.particle.surface_tension + \
                          self.particle.viscosity
    
    def debugging_forces(self):

        print("Mass: ", self.particle.mass)
        print("Mass Density:", self.particle.mass_density)
        print("Pressure:", self.particle.pressure)
        print("Pressure Force:", self.particle.pressure_force)
        print("Buoyancy:", self.particle.buoyancy)
        print("Gravity:", self.particle.gravity)
        print("surface tension:", self.particle.surface_tension)
        print("Body Force:", self.particle.body_force)
        print("Thermal:", self.particle.thermal_diffusion)
        print("viscosity:", self.particle.viscosity)
        print("\n\n")

    def choose_time_stepping(self, time_step_type:str = "Euler Cromer"):
        """
        
        """
        if self.TIME_SCHEMES[time_step_type] == 0:
            ForwardEuler(
                self.particle,
                self.delta_time
            ).exec_time_scheme(self.delta_time)
        if self.TIME_SCHEMES[time_step_type] == 1:
            EulerCromer(
                self.particle,
                self.delta_time
            ).exec_time_scheme(self.delta_time)
            print("doing time step")
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
                ).get_time_scheme_values(self.delta_time)
            if self.TIME_SCHEMES[time_step_type] == 2:
                return LeapFrog(
                    particle,
                    self.delta_time
                ).get_time_scheme_values(self.delta_time)
            if self.TIME_SCHEMES[time_step_type] == 3:
                return Verlet(
                    particle,
                    self.delta_time
                ).get_time_scheme_values(self.delta_time)
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
                               secondary_type:str = "Normal"):
            if isinstance(self.COLLISION_TYPES[collision_type], dict):
                if self.COLLISION_TYPES[collision_type][secondary_type] == 0:
                    OrientedBBox()
                if self.COLLISION_TYPES[collision_type][secondary_type] == 1:
                    AABB()
                if self.COLLISION_TYPES[collision_type][secondary_type] == 2:
                    BoxCollisions(
                        particle=self.particle,
                        tank_size=self.tank_attrs["dimensions"]["size"],
                        tank_location=self.tank_attrs["dimensions"]["location"],
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
                    BoxCollisions(
                        particle=self.particle,
                        tank_size=self.tank_attrs["dimensions"]["size"],
                        tank_location=self.tank_attrs["dimensions"]["location"],
                        speed_loss=self.PARAMETERS["loss_of_speed"]
                    ).collision_resolution()

    def update(self):

        self.update_all_forces()
        self.XSPH_vel_correction()

        self.particle.acceleration = self.all_forces / self.PARAMETERS["mass_density"]
        self.particle.next_acceleration = self.particle.acceleration

        self.choose_collision_types("Cuboid", "Normal")
        self.choose_time_stepping(self.time_stepping)

    # -------------------------------------------------------------------- THERMAL INTEGRATION ---------------------------------------------------------------------

    def update_q(self, position:np.array,
                 nbr_position:np.array):
        return np.linalg.norm(position - nbr_position) / self.PARAMETERS["cell_size"]

    def update_alpha_factor(self):
        return 1/(np.power(np.pi, 2/3) * m.pow(self.PARAMETERS["cell_size"], 3))

    def kernel_gradient_correction(self):
        pass

    def update_heat_conduction(self):

        const_term = 1 / (self.particle.pressure * self.particle.specific_heat)

        for nbr_particle in self.neighbours_list:
            thermal_conduc_term = 4 * nbr_particle.mass * self.particle.thermal_conduc * nbr_particle.thermal_conduc
            density_term = self.particle.mass_density * nbr_particle.mass_density * (
                self.particle.thermal_conduc + nbr_particle.thermal_conduc
            )
            temperature_term = self.particle.temperature - nbr_particle.temperature
            pos_dif = self.particle.initial_pos - nbr_particle.initial_pos
            kernel_grad_term = pos_dif * self.cubic_spline_kernel(1, nbr_particle.initial_pos) / (np.power(pos_dif, 2) + 
                                                                                                  m.pow(self.update_nubla(), 2))
            self.particle.thermal_diffusion += (
                (thermal_conduc_term / density_term) * temperature_term * kernel_grad_term
            )
        self.particle.thermal_diffusion *= const_term

    def update_nubla(self):
        return 0.1 * self.PARAMETERS["cell_size"]
    
    def update_laminar_viscosity(self):
        for nbr_particle in self.neighbours_list:
            viscosity_term = 4 * nbr_particle.mass * self.PARAMETERS["kinematic_visc"] * (self.particle.initial_pos - 
                                                                                          nbr_particle.initial_pos)
            velocity_term = self.particle.velocity - nbr_particle.velocity
            kernel_grad_term = self.cubic_spline_kernel(1, nbr_particle.initial_pos)
            density_term = self.particle.mass_density + nbr_particle.mass_density
            position_dif = self.particle.initial_pos - nbr_particle.initial_pos
            positional_term = (np.power(position_dif, 2) + np.power(self.update_nubla(), 2))
            self.particle.laminar_viscosity += (
                    (viscosity_term * kernel_grad_term)  / (density_term + positional_term) * velocity_term
            )

    def artificial_viscosity(self):
        pass

    def cubic_spline_kernel(self, kernel_type:int = 0,
                            nbr_position:np.array=np.array([0, 0, 0])):
        
        if kernel_type==0:
            if self.update_q(self.particle.initial_pos, nbr_position) >=0 and \
                self.update_q(self.particle.initial_pos, nbr_position) <= 1:
                    return self.update_alpha_factor() * (1 - 3/2*np.power(self.update_q(self.particle.initial_pos, nbr_position), 2) +  \
                        3/4*np.power(self.update_q(self.particle.initial_pos, nbr_position), 3) )
            if self.update_q(self.particle.initial_pos, nbr_position) >= 1 and \
                self.update_q(self.particle.initial_pos, nbr_position) <= 2:
                return 1/4 * np.power(2 - self.update_q(self.particle.initial_pos, nbr_position), 3)
            if self.update_q(self.particle.initial_pos, nbr_position) >= 2:
                return 0
        if kernel_type==1:
            return self.update_alpha_factor() * (-3 * (1- np.power(self.update_q(self.particle.initial_pos, nbr_position), 2)) *
                    ((2*self.update_q(self.particle.initial_pos, nbr_position) + 1 )/ self.PARAMETERS["cell_size"])) * \
                    self.particle.initial_pos - nbr_position
        
    def cubic_spline_kernel_grad_c(self, position, index_comp:int = 0):
        
        q = np.linalg.norm(position) / self.PARAMETERS["cell_size"]
        if q>=0 and q<=1:
            return -3*q + 9/4*m.pow(q, 2) * position[index_comp]
        if q>=1 and q<= 2:
            return -3/4 * m.pow(2 - q, 2) * position[index_comp]
        if q>2:
            return 0
        
    def cubic_spline_kernel_pos(self, position):
        q = np.linalg.norm(position) / self.PARAMETERS["cell_size"]
        if q >= 0 and q <= 0.5:
            return 1 - 6*m.pow(q, 2) + 6*m.pow(q, 3)
        if q > 0.5 and  q<= 1:
            return 2*m.pow(1-q, 3)
        if q>1:
            return 0
        
    def cubic_spline_kernel_gradient(self, position):
        q = np.linalg.norm(position) / self.PARAMETERS["cell_size"]
        if q>=0 and q<=1:
            return -3*q + 9/4*m.pow(q, 2)
        if q>=1 and q<= 2:
            return -3/4 * m.pow(2 - q, 2)
        if q>2:
            return 0

    def cubic_spline_kernel_gradient_mag(self, position):
        q = np.linalg.norm(position) / self.PARAMETERS["cell_size"]
        if q>=0 and q<=1:
            return np.abs(-3*q + 9/4*m.pow(q, 2))
        if q>=1 and q<= 2:
            return np.abs(-3/4 * m.pow(2 - q, 2))
        if q>2:
            return 0
        
    def update_body_force(self):
        self.particle.body_force = self.PARAMETERS["thermal_exp_coeff"]*self.particle.gravity* \
                (self.particle.temperature - self.PARAMETERS["abs_temperature"])
        
    def update_EOS_pressure(self):
        self.particle.pressure = self.PARAMETERS["sound_speed"] * (self.particle.mass_density - self.PARAMETERS["mass_density"])