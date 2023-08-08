import math as m
import numpy as np
import random as rd
import re
import sys

sys.path.append("C:\\Users\\Student\\OneDrive - Bournemouth University\\Desktop\\Personal\\Python\\Fluid_Predictor\\python\\Collisions\\")

from Particles.particles import Particle
from collisions import Collisions

class BoxCollisions(Collisions):
    
    def __init__(self,
                 particle: Particle = None,
                 tank_size: np.array=None,
                 tank_location: np.array=None,
                 speed_loss:float = None):
        
        self.particle = particle
        self.tank_size = tank_size
        self.tank_location = tank_location
        self.speed_loss = speed_loss
        self.collision_detected= [False, False, False, False, False, False]

    def collision_detection(self, x, y, z):
        """
            collision detection with normal tank 
            box
        """
        if x>self.tank_location[0] + self.tank_size[0]: 
            self.collision_detected[0] = True
        if x<-1*(self.tank_location[0] + self.tank_size[0]):
            self.collision_detected[1] = True
        if y>self.tank_location[1] + self.tank_size[1]:
            self.collision_detected[2] = True
        if y<-1*(self.tank_location[1] + self.tank_size[1]):
            self.collision_detected[3] = True
        if z>self.tank_location[2] + self.tank_size[2]: 
            self.collision_detected[4] = True
        if z<-1*(self.tank_location[2] + self.tank_size[2]):
            self.collision_detected[5] = True
        else:
            self.collision_detected[0] = False
            self.collision_detected[1] = False
            self.collision_detected[2] = False

    def collision_resolution_pos(self, id):
        """
            collision resolution for position of particle near boundary thresholds.
        """
        if (id+1)%2==0:
            self.particle.initial_pos[id] = (self.tank_location[id] + self.tank_size[id]) + \
                                            (-1*self.particle.initial_pos - (self.tank_location[id] 
                                                                       + self.tank_size[id]))
        elif (id+1)%2==1:
            self.particle.initial_pos[id] = (self.tank_location[id] + self.tank_size[id]) - \
                                            (self.particle.initial_pos - (self.tank_location[id] 
                                                                       + self.tank_size[id]))


    def collision_resolution(self):
        """
            tank collision resolution if collision detected in 
            one of the x, y, z axes
        """
        self.collision_detection(self.particle.initial_pos[0], self.particle.initial_pos[1], self.particle.initial_pos[2])
        for id, elem in enumerate(self.collision_detected):
            if elem is True:
                self.particle.velocity[id] = self.particle.velocity[id]*-1*self.speed_loss
                self.collision_resolution_pos(id)

class AABB(Collisions):
    pass

class OrientedBBox(Collisions):
    pass