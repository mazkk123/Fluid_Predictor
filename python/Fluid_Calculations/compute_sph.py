import math as m
import numpy as np
import sys
import time

sys.path.append("C:\\Users\\Student\\OneDrive - Bournemouth University\\Desktop\\Personal\\Python\\Fluid_Predictor\\python\\Fluid_Utilities\\")

from Fluid_Utilities.search_methods import NearestNeighbour, SpatialHashing, CompactHashing
from Fluid_Utilities.time_stepping import ForwardEuler, EulerCromer, LeapFrog, Verlet, IndependentTime, \
                                    RegionalShortTime, RungeKutta
from Collisions.box_collisions import BoxCollisions, AABB, OrientedBBox
from Collisions.sphere_collisions import SphereCollisions, CapsuleCollisions, CylinderCollisions
from Particles.particles import Particle

from utility_calculations import UtilityCalculations

class SPH(UtilityCalculations):

    PARAMETERS = {
        "grid_separation":0.1,
        "cell_size":0.4,
        "mass": 0.1,
        "viscosity": 3.5,
        "mass_density": 998.2,
        "buoyancy":0.3,
        "tension_coefficient":0.0728,
        "tension_threshold":6,
        "pressure_const":7,
        "loss_of_speed":0.5,
        "epsilon":0.1,
        "neighbour_num":30,
        "beta_const":0.25,
        "stiffness_constant":1000,
        "alpha":0.4,
        "v_cutoff":0.02,
        "N_cutoff":0.01,
        "thermal_exp_coeff":4.988,
        "kinematic_visc":0.000006,
        "lambda_const":0.005,
        "stiffness_n":1.5,
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
                 all_particles:list = None,
                 collision_type: str="box",
                 params: dict=None, 
                 collision_types: dict=None, 
                 time_schemes: dict=None):
        
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
            
        if params is not None:
            self.PARAMETERS = params
        if collision_types is not None:
            self.COLLISION_TYPES = collision_types
        if time_schemes is not None:
            self.TIME_SCHEMES = time_schemes


        super().__init__(hash_table=self.hash_table,
                         parameters=self.PARAMETERS)
        
        self.temperature = temperature

        self.neighbours_list = []
        self.dynamic_list = []

        self.active = []
        self.semi_active = []
        self.passive = []
        
        if all_particles is not None:
            self.all_particles = all_particles

        self.incremental_step = self.PARAMETERS["cell_size"] / 4

        self.normal_field()
        
        self.update_particle_neighbours()
        self.gravity_const = np.array([0, -9.81, 0], dtype="float64")
        self.all_forces = np.array([0, 0, 0], dtype="float64")

    # ------------------------------------------------------------------ PARTICLE SEARCHES ------------------------------------------------------------------------

    def update_particle_neighbours(self):
        """
            find all particles hashed to the same cell in the
            global HASH_MAP param
        """
        if self.search_method != "Neighbour":
            self.particle.neighbour_list = []
            self.mark_active_neighbours()
            for items in self.hash_table[self.hash_value]:
                if not items==self.particle:
                    if items in self.active:
                        self.neighbours_list.append(items)
                        self.particle.neighbour_list.append(items)
                    elif items in self.semi_active:
                        self.neighbours_list.append(items)
                        self.particle.neighbour_list.append(items)
                    elif items in self.other_active:
                        self.neighbours_list.append(items)
                        self.particle.neighbour_list.append(items)
            #self.particle_query()
        else:
            for items in NearestNeighbour(search_radius=self.PARAMETERS["cell_size"], particle=self.particle,
                                          neighbour_size=self.PARAMETERS["neighbour_num"]).find_neighbours():
                self.neighbours_list.append(items)
    
    def particle_query(self):
            
        bbox_max = self.particle.initial_pos + np.array([self.PARAMETERS["cell_size"],
                                                         self.PARAMETERS["cell_size"],
                                                         self.PARAMETERS["cell_size"]], dtype="float64")
                                                            
        bbox_min = self.particle.initial_pos - np.array([self.PARAMETERS["cell_size"],
                                                         self.PARAMETERS["cell_size"],
                                                         self.PARAMETERS["cell_size"]], dtype="float64")
        
        position = np.array([0, 0, 0], dtype="float64")
        for i in np.arange(bbox_min[0], bbox_max[0], self.incremental_step):
            for j in np.arange(bbox_min[1], bbox_max[1], self.incremental_step):
                for k in np.arange(bbox_min[2], bbox_max[2], self.incremental_step):
                    
                    position[0], position[1], position[2] = i, j, k
                    hash = SpatialHashing(self.PARAMETERS["cell_size"])
                    hash_value = hash.find_hash_value_pos(position)
                    
                    if hash_value != self.hash_value and \
                        hash_value in self.hash_table.keys():
                        
                        for elem in self.hash_table[hash_value]:
                            self.dynamic_list.append(elem)

        for nbrs in self.dynamic_list:
            if np.linalg.norm(self.particle.initial_pos - 
                nbrs.initial_pos) < self.PARAMETERS["cell_size"]:
                self.neighbours_list.append(nbrs)
                self.particle.neighbour_list.append(nbrs)

    def normal_field(self):
        self.particle.normal_field = np.array([0, 0, 0], dtype="float64")
        for nbr_particle in self.neighbours_list:
            try:
                mass_d = nbr_particle.mass / nbr_particle.mass_density
            except ZeroDivisionError:
                mass_d = 0
            self.particle.normal_field += (
                mass_d * self.kernel_gradient(self.particle.initial_pos - nbr_particle.initial_pos, 0)
            )

    def mark_active_members(self, particle, depth:int=3):

        if depth==0:
            return
        
        self.find_neighbour_list(particle)
        for nbr in particle.neighbour_list:
            if any(nbr.velocity) >= self.PARAMETERS["v_cutoff"] or  \
                any(nbr.normal_field) >= self.PARAMETERS["N_cutoff"]:
                self.active.append(nbr)
                return self.mark_active_members(nbr, depth-1)
            else:
                self.passive.append(nbr)
                return self.mark_active_members(nbr, depth-1)

    def mark_active_neighbours(self):
        
        self.mark_active_members(self.particle, 2)

        self.other_active = []
        for active_p in self.active:
            self.find_neighbour_list(active_p)
            for nbr in active_p.neighbour_list:
                if any(nbr.velocity) >= self.PARAMETERS["v_cutoff"]:
                    self.other_active.append(nbr)
                else:
                    self.semi_active.append(nbr)
        
                
    # ------------------------------------------------------------------- KERNEL STEPS ----------------------------------------------------------------------------

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

    # ------------------------------------------------------------------ FORCE CALCULATIONS -----------------------------------------------------------------------

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
        self.particle.gravity = self.gravity_const

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
            self.update_dynamic_viscosity()
            self.update_momentum()

        self.all_forces += self.particle.pressure_force + \
                          self.gravity_const + self.buoyancy + \
                          self.particle.surface_tension + \
                          self.particle.viscosity
    
    def debugging_forces(self, secs):

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
        time.sleep(secs)

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
                        tank_location=self.tank_attrs["dimensions"]["location"]
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
                        tank_location=self.tank_attrs["dimensions"]["location"]
                    ).collision_resolution()

    def update(self):

        self.update_all_forces()
        self.normal_field()
        self.XSPH_vel_correction()

        self.particle.acceleration = self.all_forces / self.PARAMETERS["mass_density"]
        self.particle.next_acceleration = self.particle.acceleration

        self.choose_collision_types("Cuboid", "Normal")
        self.choose_time_stepping(self.time_stepping)

        self.adapt_to_CFL()

# ------------------------------------------------------------- ADAPTIVE TIME STEPPING ----------------------------------------------------------------------------
    def update_vel_max(self):
        try:
            self.velocity_max = np.array([0, 0, 0], dtype="float64")
            self.velocity_max = [nbr_particle.velocity for nbr_particle in self.neighbours_list]
            return np.maximum.reduce(self.velocity_max)
        except ValueError:
            self.velocity_max = self.particle.velocity
    
    def update_Vel_max(self):
        self.update_vel_max()
        self.update_force_max()
        
        self.Vel_max = self.velocity_max + \
                       np.sqrt(self.PARAMETERS["cell_size"]*
                       self.force_max)
            
    def update_force_max(self):
        self.force_max = np.array([0, 0, 0], dtype="float64")
        self.max_force_arr = []
        for nbr_particle in self.neighbours_list:
            self.force_max = (
               nbr_particle.gravity +
               nbr_particle.buoyancy +
               nbr_particle.viscosity +
               nbr_particle.pressure_force +
               nbr_particle.surface_tension
            )
            self.max_force_arr.append(self.force_max)
        try:
            self.force_max = np.maximum.reduce(self.max_force_arr)
        except ValueError:
            self.force_max = self.all_forces
        
    def CFL_condition(self):
        return self.PARAMETERS["alpha"] * \
                (self.PARAMETERS["cell_size"] / np.sqrt(self.PARAMETERS["stiffness_constant"]))
    
    def CFL_force_condition(self):
        self.update_force_max()
        return self.PARAMETERS["beta_const"] * \
                np.sqrt(self.PARAMETERS["cell_size"]/ np.linalg.norm(self.force_max))
                
    def update_velocity_divergence(self):
        self.velocity_divergence = np.array([0, 0, 0], dtype="float64")
        for nbr_particle in self.neighbours_list:
            self.velocity_divergence += (nbr_particle.mass*
               self.cubic_spline_kernel_gradient(
               self.particle.initial_pos - nbr_particle.initial_pos
               )
            )
        
    def CFL_viscosity_condition(self):
        self.update_velocity_divergence()
        return self.PARAMETERS["lambda_const"] / np.linalg.norm(self.velocity_divergence)
    
    def adapt_to_CFL(self):
        self.update_Vel_max()
        if len(self.Vel_max)==1:
            if (self.Vel_max[0] < self.PARAMETERS["alpha"]*self.PARAMETERS["sound_speed"]).any():
                self.delta_time = self.PARAMETERS["stiffness_n"]*min([self.CFL_condition(), 
                                                                    self.CFL_force_condition(),
                                                                    self.CFL_viscosity_condition()])
            else:
                self.delta_time = min([self.CFL_condition(), 
                                    self.CFL_force_condition(),
                                    self.CFL_viscosity_condition()])
        elif len(self.Vel_max)==3:
            if (self.Vel_max < self.PARAMETERS["alpha"]*self.PARAMETERS["sound_speed"]).any():
                self.delta_time = self.PARAMETERS["stiffness_n"]*min([self.CFL_condition(), 
                                                                    self.CFL_force_condition(),
                                                                    self.CFL_viscosity_condition()])
            else:
                self.delta_time = min([self.CFL_condition(), 
                                    self.CFL_force_condition(),
                                    self.CFL_viscosity_condition()])
        
# ------------------------------------------------------------------ THERMAL INTEGRATION ---------------------------------------------------------------------------
    def nubla(self):
        return self.PARAMETERS["cell_size"]*0.1
    
    def update_dynamic_viscosity(self):
        self.particle.viscosity = (
           self.particle.mass_density*self.PARAMETERS["kinematic_visc"] /
           self.PARAMETERS["mass_density"]
        )
        
    def update_momentum(self):

        for nbr_particle in self.neighbours_list:
            
            pressure_term = (
                (nbr_particle.pressure / m.pow(nbr_particle.pressure, 2)) +
                (self.particle.pressure / m.pow(self.particle.pressure, 2))
            )
            
            thermal_component = (
               self.PARAMETERS["thermal_exp_coeff"] / 
               (self.particle.mass_density*nbr_particle.mass_density)
            )
            
            visc_component = (
                4*self.particle.viscosity*nbr_particle.viscosity /
                (self.particle.viscosity + nbr_particle.viscosity)
            )
            
            velocity_component = (
               (self.particle.velocity-nbr_particle.velocity)*
               (self.particle.initial_pos-nbr_particle.initial_pos) /
               (np.power((self.particle.initial_pos-nbr_particle.initial_pos), 2) +
               m.pow(self.nubla(), 2))
            )
            
            kernel_term = (
               self.cubic_spline_kernel_gradient(
                 self.particle.initial_pos - nbr_particle.initial_pos
               )
            )
            
            self.all_forces = nbr_particle.mass* \
                              (pressure_term + \
                              thermal_component* \
                              visc_component* \
                              velocity_component)* \
                              kernel_term
                            
        self.all_forces *= -1 
        self.all_forces += self.particle.gravity + \
                           self.particle.buoyancy + \
                           self.particle.surface_tension
        
    def update_body_force(self):
        self.particle.body_force = self.PARAMETERS["thermal_exp_coeff"]*self.particle.gravity* \
                (self.particle.temperature - self.PARAMETERS["abs_temperature"])
        