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
                 num_particles: int=5000):
        
        if cell_size is not None:
            self.cell_size = cell_size
        if num_particles is not None:
            self.num_particles = num_particles

    def is_prime(self, n):
        if n <= 1:
            return False
        if n <= 3:
            return True
        if n % 2 == 0 or n % 3 == 0:
            return False
        i = 5
        while i * i <= n:
            if n % i == 0 or n % (i + 2) == 0:
                return False
            i += 6
        return True

    def find_next_prime(self, num):
        num += 1
        while True:
            if self.is_prime(num):
                return num
            num += 1

    def hash_function(self, x, y, z):
        return np.array([m.floor(x/self.cell_size),
                         m.floor(y/self.cell_size),
                         m.floor(z/self.cell_size)])

    def find_hash_value(self, particle: Particle):
        
        hash_value = (
            (self.hash_function(particle.initial_pos[0],
                                particle.initial_pos[1],
                                particle.initial_pos[2])[0] * self.PRIME_NUMBERS[0]) ^ 
            (self.hash_function(particle.initial_pos[0],
                                particle.initial_pos[1],
                                particle.initial_pos[2])[1] * self.PRIME_NUMBERS[1]) ^
            (self.hash_function(particle.initial_pos[0],
                                particle.initial_pos[1],
                                particle.initial_pos[2])[2] * self.PRIME_NUMBERS[2])
            ) % self.find_next_prime(2*self.num_particles)
        return hash_value
        
    def find_hash_value_pos(self, position:np.array):
        
        hash_value = (
            (self.hash_function(position[0],
                                position[1],
                                position[2])[0] * self.PRIME_NUMBERS[0]) ^ 
            (self.hash_function(position[0],
                                position[1],
                                position[2])[1] * self.PRIME_NUMBERS[1]) ^
            (self.hash_function(position[0],
                                position[1],
                                position[2])[2] * self.PRIME_NUMBERS[2])
            ) % self.find_next_prime(2*self.num_particles)
        return hash_value
        

class NearestNeighbour:
    
    def __init__(self,
                 search_radius: float=None,
                 neighbour_size: int=None):
        
        self.search_radius = search_radius
        self.neighbour_size = neighbour_size

    def find_neighbours(self, particle: Particle=None,
                        particle_list: list=None):
        """
            find neighbouring particles given search radius
            and number of particles allowed in search radius
        """
        neighbours = []

        if Particle is not None:
            if len(particle_list)!=0:
                for p in particle_list:
                    if p!=particle and (p.initial_pos <= particle.initial_pos + self.search_radius
                                        and p.initial_pos > particle.initial_pos - self.search_radius):
                        if (len(self.neighbours) < self.neighbour_size):
                            neighbours.append(p)
        return neighbours


class CompactHashing(SpatialHashing):
    pass

class ZSorting:
    pass
