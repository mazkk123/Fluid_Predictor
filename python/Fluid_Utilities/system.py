import math as m
import numpy as np
import random as rd
import re
import sys

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
                       "LPSPH", "DFSPH", "FSISPH", "VCSPH", "PBD"]
    NEIGHBOUR_SEARCHES = ["Neighbour", "Spatial Hashing", "Compact Hashing", "Z-Sorting"]
    ORIENTATION_TYPE = ["Uniform", "Random"]
    
    TANK_ATTRS = {
        "dimensions":{
            "location":np.array([0, 0, 0], dtype="float64"), 
            "size":np.array([5, 5, 5], dtype="float64")
            },
        "type":"Cuboid",
    }

    ATTRS = {
        "cell_size":0.5
    }

    PHASE_INFORMATION = {
        "phase_number":3,
        "mass_density":[998.2, 1010, 700],
        "viscosity":[3.5, 10, 23],
        "mass":[0.1, 5, 6]
    }

    HASH_MAP = {}

    def __init__ (self,
                  type:str = "SPH",
                  search_method:str = "Spatial Hashing",
                  num_particles:int = 1000,
                  orientation_type:str = "Uniform"):

        self.simulation_type = type
        self.num_particles = num_particles
        self.orientation_type = orientation_type

        self.particle_list = []
        self.search_method = search_method
  
        self.num_frames = 100

        self.init_particle_attrs()

    def update_hash(self, particle):
        """
            updating and initializing particle hash values
            using search method specified
        """
        neighbor_search = self.NEIGHBOUR_SEARCHES[self.choose_neighbour_search()]

        if neighbor_search != "Neighbour":
            
            if neighbor_search=="Spatial Hashing":
                init_hash_value = SpatialHashing(self.ATTRS["cell_size"], self.num_particles).find_hash_value(particle)
            elif neighbor_search=="Compact Hashing":
                init_hash_value = CompactHashing(self.ATTRS["cell_size"], self.num_particles).find_hash_value(particle)

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
            id = self.choose_orientation()
            if id==0:
                for id, positions in enumerate(Uniform(num_particles = self.num_particles,
                                                       spacing=0.1).uniform_box_distribution()):
                    self.particle_list[id].initial_pos = positions
                    self.update_hash(self.particle_list[id])
            if id==1:
                for id, positions in enumerate(Random(num_particles = self.num_particles).random_box_distribution()):
                    self.particle_list[id].initial_pos = positions
                    self.update_hash(self.particle_list[id])
        
    def update(self):
        """
            update calls per particle basis
        """

        for p in self.particle_list:
            if self.choose_simulation_type() is not None:
                
                self.update_hash(p)

                """ for key, item in self.HASH_MAP.items():
                    print(key, ":", item) """

                id = self.choose_simulation_type()
                if id==0:
                    """ print("do SPH") """
                    SPH(
                        particle = p,
                        search_method = self.NEIGHBOUR_SEARCHES[self.choose_neighbour_search()],
                        hash_table = self.HASH_MAP,
                        hash_value = p.hash_value,
                        time_stepping="Euler Cromer",
                        tank_attrs = self.TANK_ATTRS,
                        all_particles=self.particle_list,
                        temperature= False,
                        delta_time=0.02
                    ).update()
                if id==1:
                    """ print("do Multi") """
                    MultiSPH(
                        particle = p,
                        search_method = self.NEIGHBOUR_SEARCHES[self.choose_neighbour_search()],
                        hash_table = self.HASH_MAP,
                        hash_value = p.hash_value,
                        time_stepping = "Euler Cromer",
                        tank_attrs = self.TANK_ATTRS,
                        delta_time=0.02,
                        phase_info=self.PHASE_INFORMATION
                    ).update()
                if id==2:
                    print("do IISPH")
                    IISPH(
                        particle = p,
                        search_method = self.NEIGHBOUR_SEARCHES[self.choose_neighbour_search()],
                        hash_table = self.HASH_MAP,
                        hash_value = p.hash_value,
                        time_stepping="Euler Cromer",
                        tank_attrs = self.TANK_ATTRS,
                        delta_time=0.02 
                    ).update()
                if id==3:
                    """ print("do WCSPH") """
                    WCSPH(
                        particle = p,
                        search_method = self.NEIGHBOUR_SEARCHES[self.choose_neighbour_search()],
                        hash_table = self.HASH_MAP,
                        hash_value = p.hash_value,
                        time_stepping="Euler Cromer",
                        tank_attrs = self.TANK_ATTRS,
                        all_particles=self.particle_list,
                        delta_time=0.02 
                    ).update()
                if id==4:
                    """ print("do PCSPH") """
                    PCSPH(
                        particle = p,
                        search_method = self.NEIGHBOUR_SEARCHES[self.choose_neighbour_search()],
                        hash_table = self.HASH_MAP,
                        hash_value = p.hash_value,
                        time_stepping="Euler Cromer",
                        all_particles=self.particle_list,
                        tank_attrs = self.TANK_ATTRS,
                        delta_time=0.02 
                    ).update()
                if id == 5:
                    """ print("do LPSPH") """
                    LPSPH(
                        particle = p,
                        search_method = self.NEIGHBOUR_SEARCHES[self.choose_neighbour_search()],
                        hash_table = self.HASH_MAP,
                        hash_value = p.hash_value,
                        time_stepping="Euler Cromer",
                        tank_attrs = self.TANK_ATTRS,
                        all_particles = self.particle_list,
                        delta_time=0.02   
                    ).update()
                if id == 6:
                    print("do DFSPH")
                    DFSPH(
                        particle = p,
                        search_method = self.NEIGHBOUR_SEARCHES[self.choose_neighbour_search()],
                        hash_table = self.HASH_MAP,
                        hash_value = p.hash_value,
                        time_stepping="Euler Cromer",
                        tank_attrs = self.TANK_ATTRS,
                        delta_time=0.02   
                    ).update()
                if id == 7:
                    print("do FSISPH")
                    FSISPH(
                        particle = p,
                        search_method = self.NEIGHBOUR_SEARCHES[self.choose_neighbour_search()],
                        hash_table = self.HASH_MAP,
                        hash_value = p.hash_value,
                        time_stepping="Euler Cromer",
                        tank_attrs = self.TANK_ATTRS,
                        delta_time=0.02   
                    ).update()
                if id == 8:
                    print("do VCSPH")
                    VCSPH(
                        particle = p,
                        search_method = self.NEIGHBOUR_SEARCHES[self.choose_neighbour_search()],
                        hash_table = self.HASH_MAP,
                        hash_value = p.hash_value,
                        time_stepping="Euler Cromer",
                        tank_attrs = self.TANK_ATTRS,
                        delta_time=0.02   
                    ).update()
                if id == 9:
                    print("do PBF")
                    PBF(
                        particle = p,
                        search_method = self.NEIGHBOUR_SEARCHES[self.choose_neighbour_search()],
                        hash_table = self.HASH_MAP,
                        hash_value = p.hash_value,
                        time_stepping="Euler Cromer",
                        tank_attrs = self.TANK_ATTRS,
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
            if and only if the solver method support neighbour
            searching, then calculate the relevant neighbour
            search algorithm
        """
        for id, orientation in enumerate(self.ORIENTATION_TYPE):
            if self.orientation_type == orientation:
                return id
