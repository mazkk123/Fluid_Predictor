import math as m
import numpy as np
import random as rd
import re
import sys
import time

from PySide6.QtWidgets import QProgressDialog, QWidget, QMainWindow
from PySide6.QtCore import QTimer

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
from search_methods import CompactHashing, SpatialHashing
from exports import ExportUtils

from Particles.particles import Particle
from Particles.error_handling import *

class FrameProgress(QProgressDialog):
    # Operation constructor
    def __init__(self, steps):
        
        super().__init__("Render progress is ...", "Cancel", 0, 100)

        self.steps = steps

        self.canceled.connect(self.cancel)

        self.t = QTimer(self)
        self.t.timeout.connect(self.perform)
        self.t.start(1000)

        self.t.start(0)

    def perform(self):
        self.setValue(self.steps)

        if self.steps > self.maximum():
            self.t.stop()
    
    def cancel(self):
        self.t.stop()

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
        "grid_separation":0.01,
        "cell_size":0.325,
        "mass": 0.1,
        "viscosity": 3.5,
        "mass_density": 998.2,
        "buoyancy":0,
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

    USER_PHASE_INFORMATION = {
        "phase_number":7,
        "mass_density":[998.2, 1010, 700, 1000, 2000, 556, 455],
        "viscosity":[3.5, 10, 23, 15, 2.4, 3.5, 7.2],
        "mass":[0.1, 0.1, 0.1]
    }

    USER_ADDITIONAL_ATTRIBUTES = {
        "density_error":0.01*USER_PARAMETERS["mass_density"],
        "max_iterations":2,
        "max_iterations_div":6,
        "divergence_error":0.6,
        "alpha_vorticity":0.5,
        "relaxation_factor":0.5
    }

    HASH_MAP = {}

    def __init__ (self, parent:QMainWindow=None,
                  type:str = "SPH",
                  search_method:str = "Spatial Hashing",
                  num_particles:int = 25000,
                  orientation_type:str = "Uniform",
                  shape_type:str = "Box",
                  time_stepping:str = "Euler Cromer"):

        self.simulation_type = type
        self.num_particles = num_particles
        self.shape_type = shape_type
        self.orientation_type = orientation_type
        self.time_stepping = time_stepping
        self.stored_positions = {}

        self.start_playforward= False
        self.start_playback = False
        self.stop = False
        self.start_play = False
        self.finished_caching = False
        self.parent = parent

        self.parent_dir = "C:\\Users\\Student\\OneDrive - Bournemouth University\\Desktop\\Personal\\Python\\Fluid_Predictor\\"
        self.export_utility = ExportUtils(self.parent_dir, "Particle", "particle", 4, "geo", 0, self)

        self.particle_list = []
        self.search_method = search_method
  
        self.num_frames = 20
        self.frame_counter = 0

        self.init_particle_attrs()

    def update_stored_positions(self):

        if not self.finished_caching:

            for i in range(self.num_frames):
                self.stored_positions[i] = []
                self.update()

                print(f"frame {i+1} complete ...")
                self.export_utility.frame_number = i

                self.export_utility.export_data()

            self.finished_caching = True
            self.start_playforward = True

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
            particle = Particle(size=0.05, shape="sphere", colour=np.array([255, 255, 255], dtype="int32"))
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
                        phase_info=self.USER_PHASE_INFORMATION
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
                        additional_params = self.USER_ADDITIONAL_ATTRIBUTES,
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
                        additional_params = self.USER_ADDITIONAL_ATTRIBUTES,
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
                        additional_params = self.USER_ADDITIONAL_ATTRIBUTES,
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
                        additional_params = self.USER_ADDITIONAL_ATTRIBUTES,
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
                        additional_params = self.USER_ADDITIONAL_ATTRIBUTES,
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
                        additional_params = self.USER_ADDITIONAL_ATTRIBUTES,
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

                self.stored_positions.setdefault(self.frame_counter, []).append(p.initial_pos.copy())
            
        self.frame_counter += 1
                

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

    
