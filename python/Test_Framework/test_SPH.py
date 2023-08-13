import math as m
import numpy as np
import sys

sys.path.append("C:\\Users\\Student\\OneDrive - Bournemouth University\\Desktop\\Personal\\Python\\Fluid_Predictor\\python\\")

from Fluid_Utilities.system import FluidSystem

system_obj = FluidSystem(type="WCSPH",search_method="Spatial Hashing", num_particles=1000)

for i in range(system_obj.num_frames):

    system_obj.update() 
    print("Position is: ", system_obj.particle_list[0].initial_pos)


