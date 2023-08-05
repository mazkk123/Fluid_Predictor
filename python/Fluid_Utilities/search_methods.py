import math as m
import numpy as np
import random as rd
import re
import sys

from Particles.particles import Particle

class SpatialHashing:
    
    PRIME_NUMBERS = [73856093, 19349663, 83492791]

    def __init__(self,
                 cell_size: float=None,
                 num_particles: int=None):
        
        self.cell_size = cell_size
        self.num_particles = num_particles

    def is_prime(self, num: int=None):
        return num%2==1

    def find_next_prime(self, num):
        
        counter = 1
        num_to_check = num

        while self.is_prime(num_to_check) is False:
            num_to_check += counter
        
        return num

    def hash_function(self, x, y, z):
        return np.array([m.floor(x/self.cell_size),
                         m.floor(y/self.cell_size),
                         m.floor(z/self.cell_size)])

    def find_hash_value(self, particle: Particle):
        
        hash_value = (
            (self.hash_function[0] * self.PRIME_NUMBERS[0]) ^ 
            (self.hash_function[1] * self.PRIME_NUMBERS[1]) ^
            (self.hash_function[2] * self.PRIME_NUMBERS[2])
            ) % self.find_next_prime(2*self.num_particles)
        

class NearestNeighbour:
    
    def __init__(self,
                 search_radius: float=None,
                 neighbour_size: int=None):
        
        self.search_radius = search_radius
        self.neighbour_size = neighbour_size
        self.neighbours = []

    def find_neighbours(self, particle: Particle=None,
                        particle_list: list=None):
        """
            find neighbouring particles given search radius
            and number of particles allowed in search radius
        """
        self.neighbours = []

        if Particle is not None:
            if len(particle_list)!=0:
                for p in particle_list:
                    if p!=particle and (p.initial_pos <= particle.initial_pos + self.search_radius
                                        and p.initial_pos > particle.initial_pos - self.search_radius):
                        if (len(self.neighbours) < self.neighbour_size):
                            self.neighbours.append(p)
        return self.neighbours


class CompactHashing(SpatialHashing):
    pass

class ZSorting:
    pass