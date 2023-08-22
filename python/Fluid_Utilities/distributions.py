import math as m
import numpy as np
import random as rd
import re
import sys

class Uniform:

    def __init__(self, num_particles:int = None,
                 spacing: float=None, radius:float = None,
                 height:int = None):

        self.num_particles = num_particles
        if spacing is not None:
            self.spacing = spacing
        if height is not None:
            self.height = height
        if radius is not None:
            self.radius = radius

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
        num_circles = int(m.pow(self.num_particles, 1/4))
        num_points_per_circle = int(m.pow(self.num_particles, 1/4))

        radii = np.linspace(0, self.radius, num_circles)
    
        # Generate angles for points on each circle
        angles = np.linspace(0, 2 * np.pi, num_points_per_circle)
        
        # Generate heights
        heights = np.linspace(0, self.height, num_circles)
        
        # Initialize empty lists to store points
        points = []
        
        # Generate points on each circle
        for r in radii:
            for h in heights:
                for theta in angles:
                    x = r * np.cos(theta)
                    y = r * np.sin(theta)
                    z = h
                    points.append(np.array([x, y, z]))
        
        return points
    
    def uniform_sphere_distribution(self):
        """
            perform a uniform distribution of particles and returns
            a list of np arrays holding the particle positions in a list
        """        
        
        phi = np.linspace(0, np.pi, int(m.pow(self.num_particles, 1/2)))
        theta = np.linspace(0, 2 * np.pi, int(m.pow(self.num_particles, 1/2)))
        phi, theta = np.meshgrid(phi, theta)
        
        x = np.sin(phi) * np.cos(theta)
        y = np.sin(phi) * np.sin(theta)
        z = np.cos(phi)
        
        points = np.column_stack((x.ravel(), y.ravel(), z.ravel()))
        return points

class Random:

    def __init__(self, num_particles:int = None, 
                 height:int = None,
                 radius:float = None):

        self.num_particles = num_particles
        
        if height is not None:  
            self.height = height
        if radius is not None:
            self.radius = radius

    def random_box_distribution(self):
        """
            perform a uniform distribution of particles and returns
            a list of np arrays holding the particle positions in a list
        """        

        positions = []

        for i in range(self.num_particles):
            
            x = rd.randint(0, self.num_particles) + rd.random()
            y = rd.randint(0, self.num_particles) + rd.random()
            z = rd.randint(0, self.num_particles) + rd.random()

            positions.append(np.array([x, y, z] ,dtype="float64"))

        return positions
        
    def random_cylinder_distribution(self):

        radii = np.random.uniform(0, self.radius, self.num_particles)
        angles = np.random.uniform(0, 2 * np.pi, self.num_particles)
        heights = np.random.uniform(0, self.height, self.num_particles)
        
        # Convert polar coordinates to Cartesian coordinates
        x = radii * np.cos(angles)
        y = radii * np.sin(angles)
        z = heights
        
        return np.column_stack((x.ravel(), y.ravel(), z.ravel()))
        
    def random_sphere_distribution(self):

        positions = []

        # Generate random spherical coordinates
        inclinations = np.arccos(2 * np.random.rand(self.num_particles) - 1)
        azimuths = 2 * np.pi * np.random.rand(self.num_particles)

        # Convert spherical coordinates to Cartesian coordinates
        for id, inclinations in enumerate(inclinations):
            x = self.radius * np.sin(inclinations) * np.cos(azimuths[id])
            y = self.radius * np.sin(inclinations) * np.sin(azimuths[id])
            z = self.radius * np.cos(inclinations)
            position = np.array([x, y, z])
            positions.append(position)

        return positions