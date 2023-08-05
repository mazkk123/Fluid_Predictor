import math as m
import numpy as np
import random as rd
import re
import sys

from Fluid_Calculations.compute_Euler import Euler
from Fluid_Calculations.compute_SPH import SPH
from Fluid_Calculations.compute_multiSPH import MultiSPH
from Fluid_Calculations.compute_PCSPH import PCSPH
from Fluid_Calculations.compute_WCSPH import WCSPH
from Fluid_Calculations.compute_FLIP import FLIP
from Fluid_Calculations.compute_LPSPH import LPSPH
from Fluid_Calculations.compute_IISPH import IISPH

from distributions import Random, Uniform
from search_methods import *

from Particles.particles import Particle
from Particles.error_handling import *

class FluidSystem:
    
    SIMULATION_TYPES = ["SPH", "Euler", "MultiSPH", "IISPH", 
                        "FLIP", "WCSPH", "PCSPH", "LPSPH"]
    NEIGHBOUR_SEARCHES = ["Neighbour", "Spatial Hashing", "Compact Hashing",
                          "Z-Sorting"]
    ORIENTATION_TYPE = ["Uniform", "Random"]
    
    TANK_ATTRS = {
        "dimensions":{
            "location":np.array([0, 0, 0]), "size":np.array([5, 5, 5])
            },
        "type":"Cuboid"
    }

    HASH_MAP = {}

    def __init__ (self,
                  type:str = "SPH",
                  search_method:str = "Spatial Hashing",
                  num_particles:int = 1000,
                  solver_type:str = "p",
                  orientation_type:str = "Uniform"):
        
        super(FluidSystem, self).__init__()

        self.simulation_type = type
        self.solver_type = solver_type
        self.num_particles = num_particles
        self.orientation_type = orientation_type

        self.particle_list = []
        self.search_method = search_method
  
        self.do_neighbouring = False

        self.init_particle_attrs()

    def update_hash(self, particle):
        """
            updating and initializing particle hash values
            using search method specified
        """
        init_hash_value = SpatialHashing(self.PARAMETERS["cell_size", self.num_particles]).find_hash_value(particle)
        if len(self.HASH_MAP[init_hash_value])!=0:
            self.HASH_MAP[init_hash_value].append(particle)
            particle.hash_value = init_hash_value
        else:
            self.HASH_MAP[init_hash_value] = [particle]
            particle.hash_value = init_hash_value

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
                id = self.choose_simulation_type()
                if id==0:
                    print("do SPH")
                    SPH(
                        particle = p,
                        search_method = self.NEIGHBOUR_SEARCHES[self.choose_neighbour_search()],
                        hash_table = self.HASH_MAP,
                        hash_value = p.hash_value,
                        time_stepping="Euler Cromer",
                        tank_attrs = self.TANK_ATTRS,
                        delta_time=0.02
                    ).update()
                if id==1:
                    print("do Euler")
                    Euler()
                if id==2:
                    print("do Multi")
                    MultiSPH()
                if id==3:
                    print("do IISPH")
                    IISPH()
                if id==4:
                    print("do FLIP")
                    FLIP()
                if id==5:
                    print("do WCSPH")
                    WCSPH(
                        particle = p,
                        search_method = self.NEIGHBOUR_SEARCHES[self.choose_neighbour_search()],
                        hash_table = self.HASH_MAP,
                        hash_value = p.hash_value,
                        time_stepping="Euler Cromer",
                        tank_attrs = self.TANK_ATTRS,
                        delta_time=0.02 
                    ).update()
                if id==6:
                    print("do PCSPH")
                    PCSPH(
                        particle = p,
                        search_method = self.NEIGHBOUR_SEARCHES[self.choose_neighbour_search()],
                        hash_table = self.HASH_MAP,
                        hash_value = p.hash_value,
                        time_stepping="Euler Cromer",
                        tank_attrs = self.TANK_ATTRS,
                        delta_time=0.02 
                    )
                if id == 7:
                    print("do LPSPH")
                    LPSPH()

    def choose_simulation_type(self):
        """
            chooses the relevant simulation type to
            do further computations on.
        """
        for id, sim_type in enumerate(self.SIMULATION_TYPES):
            if self.simulation_type == sim_type:
                if self.solver_type == "p":
                    self.do_neighbouring = True 
                else:
                    self.do_neighbouring = False
                return id           
            else:
                self.sim_active[id] = False
                return None

    def choose_neighbour_search(self):
        """
            if and only if the solver method support neighbour
            searching, then calculate the relevant neighbour
            search algorithm
        """
        if self.do_neighbouring:
            for id, neighbr_search in enumerate(self.NEIGHBOUR_SEARCHES):
                if self.search_method == neighbr_search:
                    return id
                else:
                    return None

    def choose_orientation(self):
        """
            if and only if the solver method support neighbour
            searching, then calculate the relevant neighbour
            search algorithm
        """
        for id, orientation in enumerate(self.ORIENTATION_TYPE):
            if self.orientation_type == orientation:
                return id
            else:
                return None

