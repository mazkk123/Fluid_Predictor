import math as m
import numpy as np
import random as rd
import re

from Particles.particles import Particle
from collisions import Collisions

class CapsuleCollisions(Collisions):
    def __init__(self):
        
        super().__init__()
        
        self.speed_loss = self.ATTRIBS["speed_loss"]
        self.collision_detected= [False, False, False, False, False, False]

class SphereCollisions(Collisions):
    def __init__(self):
        
        super().__init__()
        
        self.speed_loss = self.ATTRIBS["speed_loss"]
        self.collision_detected= [False, False, False, False, False, False]
    
    def detection(self):
        return self.radius > np.linalg.norm(
            self.particle.initial_pos - self.tank_location
        ) 

class CylinderCollisions(Collisions):
    def __init__(self):
        
        super().__init__()
        
        self.speed_loss = self.ATTRIBS["speed_loss"]
        self.collision_detected= [False, False, False, False, False, False]


