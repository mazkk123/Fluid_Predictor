import math as m
import numpy as np
import random as rd
import re

from Particles.particles import Particle

class Collisions:

    ATTRIBS = {
            "speed_loss":0.5
    }

    def __init__(self,
                 particle:Particle,
                 tank_size:np.array,
                 tank_location:np.array):

        self.particle = particle
        self.tank_size = tank_size
        self.tank_location = tank_location