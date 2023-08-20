import math as m
import numpy as np
import sys
import time

sys.path.append("C:\\Users\\Student\\OneDrive - Bournemouth University\\Desktop\\Personal\\Python\\Fluid_Predictor\\python\\")

from Fluid_Utilities.system import FluidSystem

system_obj = FluidSystem(type="MultiSPH", search_method="Spatial Hashing", num_particles=1000)

for i in range(system_obj.num_frames):
    
    print("Position particle 1 is: ", system_obj.particle_list[0].initial_pos)
    system_obj.update() 


