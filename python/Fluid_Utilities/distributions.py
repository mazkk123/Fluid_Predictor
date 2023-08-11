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
        self.num_layers = 10
        self.radius = 5

    def uniform_box_distribution(self):
        """
            perform a uniform distribution of particles and returns
            a list of np arrays holding the particle positions in a list
        """        

        ppr = int(m.pow(self.num_particles, 1/3))
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
        num_particles_per_layer = self.num_particles / self.num_layers
        positions = []
        
        for i in range(self.num_layers):
            rho = np.linspace(0, 2*np.pi, self.num_particles_per_layer, 
                              endpoint=False)
            position = np.array([self.radius*np.cos(rho),
                                 self.radius*np.sin(rho),
                                 i])
            positions.append(position)

    def uniform_sphere_distribution(self):
        """
            perform a uniform distribution of particles and returns
            a list of np arrays holding the particle positions in a list
        """        
        
        positions = []
        
        rho = np.linspace(0, 2*np.pi, self.num_particles / 2, 
                              endpoint=False)
        phi = np.linspace(0, np.pi, self.num_particles / 2, 
                              endpoint=False)
        position = np.array([self.radius*np.cos(rho)*np.sin(phi)
                             self.radius*np.sin(rho)*np.sim(phi)
                             self.radius*np.cos(phi)])
        positions.append(position)

class Random:

    def __init__(self, num_particles:int = None):

        self.num_particles = num_particles

    def random_box_distribution(self):
        """
            perform a uniform distribution of particles and returns
            a list of np arrays holding the particle positions in a list
        """        

        positions = []
        self.radius = 5
        self.num_layers = 10 

        for i in range(self.num_particles):
            
            x = rd.randint(self.num_particles) + rd.random()
            y = rd.randint(self.num_particles) + rd.random()
            z = rd.randint(self.num_particles) + rd.random()

            positions.append(np.array([x, y, z]))

        return positions
        
    def random_cylinder_distribution(self):
        num_particles_per_layer = self.num_particles / self.num_layers
        positions = []
        
        for i in range(self.num_layers):
            for j in range(num_particles_per_layer):
                rho = j*rd.random()*2*np.pi
                position = np.array([self.radius*np.cos(rho),
                                     self.radius*np.sin(rho),
                                     i])
            positions.append(position)
        
    def random_sphere_distribution(self):
        positions = []
        
        for i in range(self.num_particles / 2):
            for j in range(self.num_particles / 2):
                rho = j*rd.random()*2*np.pi
                phi = i*rd.random()*np.pi
                position = np.array([self.radius*np.cos(rho)*np.sin(phi),
                                     self.radius*np.sin(rho)*np.sin(phi),
                                     self.radius*np.cos(phi)])
            positions.append(position)
