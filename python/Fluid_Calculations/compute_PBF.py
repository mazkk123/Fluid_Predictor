import math as m
import numpy as np
import random as rd
import re

from Fluid_Calculations.compute_sph import SPH
from Particles.particles import Particle

class PBF(SPH):
    
    def __init__(self,
                particle: Particle=None,
                search_method:str = None,
                hash_table:dict = None,
                hash_value:int = None,
                time_stepping:str = "Euler Cromer",
                tank_attrs:dict = None,
                delta_time:int = 0.02):
        
        super().__init__(particle=particle,
                         search_method=search_method,
                         hash_table=hash_table,
                         hash_value=hash_value,
                         time_stepping=time_stepping,
                         tank_attrs=tank_attrs,
                         delta_time=delta_time)