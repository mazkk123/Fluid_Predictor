import math as m
import numpy as np
import random as rd
import re

from collisions import Collisions

class BoxCollisions(Collisions):
    
    def __init__(self,
                 tank_size: np.array=None,
                 tank_x: float=None, tank_y: float=None,
                 tank_z: float=None, tank_width: float=None,
                 tank_height: float=None, tank_depth: float=None):
        
        self.tank_size = tank_size
        self.tank_x = tank_x
        self.tank_y = tank_y
        self.tank_z = tank_z
        self.tank_width = tank_width
        self.tank_height = tank_height
        self.tank_depth = tank_depth
        self.collision_detected= [False, False, False]

    def collision_detection(self, x, y, z):
        """
            collision detection with normal tank 
            box
        """
        if x>self.tank_x + self.tank_width or x<self.tank_x:
            self.collision_detected[0] = True
        elif y>self.tank_y + self.tank_height or y<self.tank_y:
            self.collision_detected[1] = True
        elif z>self.tank_z + self.tank_depth or z<self.tank_z:
            self.collision_detected[2] = True
        else:
            self.collision_detected[0] = False
            self.collision_detected[1] = False
            self.collision_detected[2] = False

    def collision_resolution(self, particle_vel: np.array=None):
        """
            tank collision resolution if collision detected in 
            one of the x, y, z axes
        """
        new_particle_pos, new_particle_vel = np.array([0, 0, 0]), np.array([0, 0, 0])
        for id, elem in enumerate(self.collision_detected):
            if elem is True:
                new_particle_vel[id] = particle_vel[id]*-1
        return new_particle_vel

class AABB(Collisions):
    pass

class OrientedBBox(Collisions):
    pass