import math as m
import numpy as np
import random as rd
import re
import sys
import time

sys.path.append("C:\\Users\\Student\\OneDrive - Bournemouth University\\Desktop\\Personal\\Python\\Fluid_Predictor\\python\\Fluid_Calculations\\")

from Fluid_Calculations.compute_sph import SPH
from Fluid_Calculations.compute_multiSPH import MultiSPH
from Fluid_Calculations.compute_PCSPH import PCSPH
from Fluid_Calculations.compute_WCSPH import WCSPH
from Fluid_Calculations.compute_LPSPH import LPSPH
from Fluid_Calculations.compute_IISPH import IISPH
from Fluid_Calculations.compute_DFSPH import DFSPH
from Fluid_Calculations.compute_FSISPH import FSISPH
from Fluid_Calculations.compute_VCSPH import VCSPH
from Fluid_Calculations.compute_PBF import PBF


sys.path.append("C:\\Users\\Student\\OneDrive - Bournemouth University\\Desktop\\Personal\\Python\\Fluid_Predictor\\python\\Fluid_Utilities\\")

from distributions import Random, Uniform
from search_methods import CompactHashing, SpatialHashing, ZSorting

from Particles.particles import Particle
from Particles.error_handling import *

class FluidSystem:
    
    SIMULATION_TYPES = ["SPH", "MultiSPH", "IISPH", "WCSPH", "PCSPH", 
                       "LPSPH", "DFSPH", "VCSPH", "FSISPH", "PBF"]
    NEIGHBOUR_SEARCHES = ["Neighbour", "Spatial Hashing", "Compact Hashing", "Z-Sorting"]

    ORIENTATION_TYPE = {
        "Uniform":{
            "id":0,
            "data":["Box", "Cylinder", "Sphere"]
        },
        "Random":{
            "id":1,
            "data":["Box", "Cylinder", "Sphere"]
            }
        }
    
    TANK_ATTRS = {
        "dimensions":{
            "location":np.array([0, 0, 0], dtype="float64"), 
            "size":np.array([5, 5, 5], dtype="float64")
            },
        "type":"Cuboid",
    }

    USER_PARAMETERS = {
        "initial_velocity":np.array([0, 0, 0], dtype="float64"),
        "initial_acceleration":np.array([0, 0, 0], dtype="float64"),
        "radius":2.5,
        "height":0.15,
        "grid_separation":0.005,
        "cell_size":0.4,
        "mass": 0.1,
        "viscosity": 3.5,
        "mass_density": 998.2,
        "buoyancy":0.2,
        "tension_coefficient":0.0728,
        "tension_threshold":6,
        "pressure_const":7,
        "loss_of_speed":0.5,
        "epsilon":0.1,
        "neighbour_num":30,
        "beta_const":0,
        "stiffness_constant":1000,
        "alpha":0.4,
        "v_cutoff":0,
        "N_cutoff":0,
        "thermal_exp_coeff":4.988,
        "kinematic_visc":0.000006,
        "lambda_const":0.005,
        "stiffness_n":1.5,
        "sound_speed":300
    }

    USER_TIME_SCHEMES = {
        "Forward Euler": 0,
        "Euler Cromer": 1,
        "Leap Frog": 2,
        "Verlet": 3,
        "Independent Time": 4,
        "Regional Short Time": 5,
        "Runge Kutta": 6
    }

    USER_COLLISION_TYPES = {
        "Cuboid":{"OrientedBBox":0, "AABB":1, "Normal":2},
        "Cylinder":1,
        "Sphere":2,
        "Capsule":3,
        "Abstract":4
    }

    PHASE_INFORMATION = {
        "phase_number":7,
        "mass_density":[998.2, 1010, 700, 1000, 2000, 556, 455],
        "viscosity":[3.5, 10, 23, 15, 2.4, 3.5, 7.2],
        "mass":[0.1, 0.1, 0.1]
    }

    HASH_MAP = {}

    def __init__ (self,
                  type:str = "SPH",
                  search_method:str = "Spatial Hashing",
                  num_particles:int = 10000,
                  orientation_type:str = "Uniform",
                  shape_type:str = "Box",
                  time_stepping:str = "Euler Cromer"):

        self.simulation_type = type
        self.num_particles = num_particles
        self.shape_type = shape_type
        self.orientation_type = orientation_type
        self.time_stepping = time_stepping

        self.particle_list = []
        self.search_method = search_method
  
        self.num_frames = 24

        self.init_particle_attrs()

    def update_hash(self, particle):
        """
            updating and initializing particle hash values
            using search method specified
        """
        neighbor_search = self.NEIGHBOUR_SEARCHES[self.choose_neighbour_search()]

        if neighbor_search != "Neighbour":
            
            if neighbor_search=="Spatial Hashing":
                init_hash_value = SpatialHashing(self.USER_PARAMETERS["cell_size"], self.num_particles).find_hash_value(particle)
            elif neighbor_search=="Compact Hashing":
                init_hash_value = CompactHashing(self.USER_PARAMETERS["cell_size"], self.num_particles).find_hash_value(particle)

            try:
                self.HASH_MAP[init_hash_value].append(particle)
                particle.hash_value = init_hash_value
            except KeyError:
                self.HASH_MAP[init_hash_value] = [particle]
                particle.hash_value = init_hash_value
        else:
            return

    def init_particle_attrs(self):
        """
            initalize particles with physical attributes
            when the simulation starts
        """
        for p in range(self.num_particles):
            particle = Particle(size=1, shape="sphere")
            self.particle_list.append(particle)
    
        if self.choose_orientation() is not None:
            orientation = self.choose_orientation()
            sub_id = int(self.choose_shape(orientation[0]))
            if orientation[1]==0:
                if sub_id == 0:
                    for id, positions in enumerate(Uniform(num_particles = self.num_particles,
                                                           spacing=self.USER_PARAMETERS["grid_separation"]).uniform_box_distribution()):
                        self.particle_list[id].initial_pos = positions
                        self.update_hash(self.particle_list[id])
                        self.set_init_attrs(id)
                elif sub_id == 1:
                    for id, positions in enumerate(Uniform(num_particles = self.num_particles,
                                                           radius=self.USER_PARAMETERS["radius"], height=0.2).uniform_cylinder_distribution()):
                        self.particle_list[id].initial_pos = positions
                        self.update_hash(self.particle_list[id])
                        self.set_init_attrs(id)
                elif sub_id == 2:
                    for id, positions in enumerate(Uniform(num_particles = self.num_particles,
                                                           radius=self.USER_PARAMETERS["radius"]).uniform_sphere_distribution()):
                        self.particle_list[id].initial_pos = positions
                        self.update_hash(self.particle_list[id])
                        self.set_init_attrs(id)

            if orientation[1]==1:
                if sub_id == 0:
                    for id, positions in enumerate(Random(num_particles = self.num_particles).random_box_distribution()):
                        self.particle_list[id].initial_pos = positions
                        self.update_hash(self.particle_list[id])
                        self.set_init_attrs(id)
                elif sub_id == 1:
                    for id, positions in enumerate(Random(num_particles = self.num_particles,
                                                          radius=self.USER_PARAMETERS["radius"],
                                                          height=self.USER_PARAMETERS["height"]).random_cylinder_distribution()):
                        self.particle_list[id].initial_pos = positions
                        self.update_hash(self.particle_list[id])
                        self.set_init_attrs(id)
                elif sub_id == 2:
                    for id, positions in enumerate(Random(num_particles = self.num_particles,
                                                          radius=self.USER_PARAMETERS["radius"]).random_sphere_distribution()):
                        self.particle_list[id].initial_pos = positions
                        self.update_hash(self.particle_list[id])
                        self.set_init_attrs(id)

    def set_init_attrs(self, id):
        self.particle_list[id].velocity = self.USER_PARAMETERS["initial_velocity"]
        self.particle_list[id].acceleration = self.USER_PARAMETERS["initial_acceleration"]

    def update(self):
        """
            update calls per particle basis
        """

        for p in self.particle_list:
            if self.choose_simulation_type() is not None:
                
                self.update_hash(p)

                id = self.choose_simulation_type()
                if id==0:
                    SPH(
                        particle = p,
                        search_method = self.NEIGHBOUR_SEARCHES[self.choose_neighbour_search()],
                        hash_table = self.HASH_MAP,
                        hash_value = p.hash_value,
                        time_stepping=self.time_stepping,
                        tank_attrs = self.TANK_ATTRS,
                        all_particles=self.particle_list,
                        collision_types=self.USER_COLLISION_TYPES,
                        params=self.USER_PARAMETERS,
                        time_schemes=self.USER_TIME_SCHEMES,
                        temperature= False,
                        delta_time=0.02
                    ).update()
                if id==1:
                    MultiSPH(
                        particle = p,
                        search_method = self.NEIGHBOUR_SEARCHES[self.choose_neighbour_search()],
                        hash_table = self.HASH_MAP,
                        hash_value = p.hash_value,
                        time_stepping = self.time_stepping,
                        all_particles = self.particle_list,
                        collision_types=self.USER_COLLISION_TYPES,
                        params=self.USER_PARAMETERS,
                        time_schemes=self.USER_TIME_SCHEMES,
                        tank_attrs = self.TANK_ATTRS,
                        delta_time=0.02,
                        phase_info=self.PHASE_INFORMATION
                    ).update()
                if id==2:
                    IISPH(
                        particle = p,
                        search_method = self.NEIGHBOUR_SEARCHES[self.choose_neighbour_search()],
                        hash_table = self.HASH_MAP,
                        all_particles=self.particle_list,
                        collision_types=self.USER_COLLISION_TYPES,
                        params=self.USER_PARAMETERS,
                        time_schemes=self.USER_TIME_SCHEMES,
                        hash_value = p.hash_value,
                        time_stepping=self.time_stepping,
                        tank_attrs = self.TANK_ATTRS,
                        delta_time=0.02 
                    ).update()
                if id==3:
                    WCSPH(
                        particle = p,
                        search_method = self.NEIGHBOUR_SEARCHES[self.choose_neighbour_search()],
                        hash_table = self.HASH_MAP,
                        hash_value = p.hash_value,
                        time_stepping=self.time_stepping,
                        tank_attrs = self.TANK_ATTRS,
                        all_particles=self.particle_list,
                        collision_types=self.USER_COLLISION_TYPES,
                        params=self.USER_PARAMETERS,
                        time_schemes=self.USER_TIME_SCHEMES,
                        delta_time=0.02 
                    ).update()
                if id==4:
                    PCSPH(
                        particle = p,
                        search_method = self.NEIGHBOUR_SEARCHES[self.choose_neighbour_search()],
                        hash_table = self.HASH_MAP,
                        hash_value = p.hash_value,
                        time_stepping=self.time_stepping,
                        tank_attrs = self.TANK_ATTRS,
                        all_particles=self.particle_list,
                        collision_types=self.USER_COLLISION_TYPES,
                        params=self.USER_PARAMETERS,
                        time_schemes=self.USER_TIME_SCHEMES,
                        temperature=False,
                        delta_time=0.02 
                    ).update()
                if id == 5:
                    LPSPH(
                        particle = p,
                        search_method = self.NEIGHBOUR_SEARCHES[self.choose_neighbour_search()],
                        hash_table = self.HASH_MAP,
                        hash_value = p.hash_value,
                        time_stepping=self.time_stepping,
                        tank_attrs = self.TANK_ATTRS,
                        all_particles = self.particle_list,
                        collision_types=self.USER_COLLISION_TYPES,
                        params=self.USER_PARAMETERS,
                        time_schemes=self.USER_TIME_SCHEMES,
                        delta_time=0.02   
                    ).update()
                if id == 6:
                    DFSPH(
                        particle = p,
                        search_method = self.NEIGHBOUR_SEARCHES[self.choose_neighbour_search()],
                        hash_table = self.HASH_MAP,
                        hash_value = p.hash_value,
                        all_particles=self.particle_list,
                        collision_types=self.USER_COLLISION_TYPES,
                        params=self.USER_PARAMETERS,
                        time_schemes=self.USER_TIME_SCHEMES,
                        time_stepping=self.time_stepping,
                        tank_attrs = self.TANK_ATTRS,
                        num_particles=self.num_particles,
                        delta_time=0.02   
                    ).update()
                if id == 7:
                    VCSPH(
                        particle = p,
                        search_method = self.NEIGHBOUR_SEARCHES[self.choose_neighbour_search()],
                        hash_table = self.HASH_MAP,
                        hash_value = p.hash_value,
                        num_particles=self.num_particles,
                        all_particles=self.particle_list,
                        collision_types=self.USER_COLLISION_TYPES,
                        params=self.USER_PARAMETERS,
                        time_schemes=self.USER_TIME_SCHEMES,
                        time_stepping=self.time_stepping,
                        tank_attrs = self.TANK_ATTRS,
                        delta_time=0.02   
                    ).update()
                if id == 8:
                    FSISPH(
                        particle = p,
                        search_method = self.NEIGHBOUR_SEARCHES[self.choose_neighbour_search()],
                        hash_table = self.HASH_MAP,
                        hash_value = p.hash_value,
                        time_stepping=self.time_stepping,
                        tank_attrs = self.TANK_ATTRS,
                        collision_types=self.USER_COLLISION_TYPES,
                        params=self.USER_PARAMETERS,
                        time_schemes=self.USER_TIME_SCHEMES,
                        delta_time=0.02   
                    ).update()
                if id == 9:
                    PBF(
                        particle = p,
                        search_method = self.NEIGHBOUR_SEARCHES[self.choose_neighbour_search()],
                        hash_table = self.HASH_MAP,
                        hash_value = p.hash_value,
                        time_stepping=self.time_stepping,
                        tank_attrs = self.TANK_ATTRS,
                        collision_types=self.USER_COLLISION_TYPES,
                        params=self.USER_PARAMETERS,
                        time_schemes=self.USER_TIME_SCHEMES,
                        delta_time=0.02   
                    ).update()

    def choose_simulation_type(self):
        """
            chooses the relevant simulation type to
            do further computations on.
        """
        for id, sim_type in enumerate(self.SIMULATION_TYPES):
            if self.simulation_type == sim_type:
                return id           

    def choose_neighbour_search(self):
        """
            if and only if the solver method support neighbour
            searching, then calculate the relevant neighbour
            search algorithm
        """
        for id, neighbr_search in enumerate(self.NEIGHBOUR_SEARCHES):
            if self.search_method == neighbr_search:
                return id

    def choose_orientation(self):
        """
            if and only if the solver method supports neighbour
            searching, then calculate the relevant neighbour
            search algorithm
        """
        for distr_type, shapes in self.ORIENTATION_TYPE.items():
            if self.orientation_type == distr_type:
                return (distr_type, self.ORIENTATION_TYPE[distr_type]["id"])
            
    def choose_shape(self, distr_type):
        """
            if and only if the solver method supports neighbour
            searching, then calculate the relevant neighbour
            search algorithm
        """ 
        for id, shapes in self.ORIENTATION_TYPE[distr_type].items():
            if isinstance(shapes, list):
                for shape in shapes:
                    if self.shape_type == shape:
                        return shapes.index(shape)
            else:
                continue

    
