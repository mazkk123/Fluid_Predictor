import math as m
import numpy as np
import random as rd
import re
import sys

class Uniform:

    def __init__(self, num_particles:int = None,
                 spacing: float=None):

        self.num_particles = num_particles
        self.spacing = spacing

    def uniform_box_distribution(self):
        """
            perform a uniform distribution of particles and returns
            a list of np arrays holding the particle positions in a list
        """        

        ppr = m.pow(self.num_particles, 1/3)
        positions = []

        for i in range(ppr):
            for j in range(ppr):
                for k in range(ppr):
                    position = np.array([i + self.spacing,
                                         j + self.spacing,
                                         k + self.spacing])
                    positions.append(position)

        return positions

    def uniform_cylinder_distribution(self):
        """
            perform a uniform distribution of particles and returns
            a list of np arrays holding the particle positions in a list
        """        
        pass

    def uniform_sphere_distribution(self):
        """
            perform a uniform distribution of particles and returns
            a list of np arrays holding the particle positions in a list
        """        
        pass

class Random:

    def __init__(self, num_particles:int = None):

        self.num_particles = num_particles

    def random_box_distribution(self):
        """
            perform a uniform distribution of particles and returns
            a list of np arrays holding the particle positions in a list
        """        

        positions = []

        for i in range(self.num_particles):
            
            x = rd.randint(self.num_particles) + rd.random()
            y = rd.randint(self.num_particles) + rd.random()
            z = rd.randint(self.num_particles) + rd.random()

            positions.append(np.array([x, y, z]))

        return positions
