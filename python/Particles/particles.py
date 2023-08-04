import numpy as np
import math as m
import random as rd

from typing import Union, Optional

class Particle:

    def __init__(self,
                 init_pos : Union(float, float, float) | np.array=None,
                 colour : Union(int, int, int, int) | np.array=None,
                 velocity : Union(float, float, float) | np.array=None,
                 acceleration: Union(float, float, float)| np.array=None,
                 mass: float=None, shape : str=None, size : float=None) -> None:
        """
            initialize particle properties 
        """
        self._initial_pos = init_pos
        self._colour = colour
        self._velocity = velocity
        self._shape = shape
        self._size = size
        self._acceleration = acceleration
        self._mass = mass
        self._hash_value = 0

        # particle force attributes

        self._mass_density = 0
        self._pressure = 0
        self._viscosity = np.array([0, 0, 0])
        self._pressure_force = np.array([0, 0, 0])
        self._gravity = np.array([0, 0, 0])
        self._buoyancy = np.array([0, 0, 0])
        self._surface_tension = np.array([0, 0 ,0])

# ------------------------------------ PROPERTY GETTERS --------------------------------------
    @property
    def initial_pos(self):
        return self._initial_pos
    
    @property
    def colour(self):
        return self._colour
    
    @property
    def velocity(self):
        return self._velocity
    
    @property
    def acceleration(self):
        return self._acceleration
    
    @property
    def shape(self):
        return self._shape
    
    @property
    def neighbours(self):
        return self._neighbours
    
    @property
    def size(self):
        return self._size
    
    @property
    def mass_density(self):
        return self._mass_density
    
    @property
    def pressure(self):
        return self._pressure
    
    @property
    def viscosity(self):
        return self._viscosity
    
    @property
    def pressure_force(self):
        return self._pressure_force
    
    @property
    def gravity(self):
        return self._gravity
    
    @property
    def buoyancy(self):
        return self._buoyancy
    
    @property
    def mass(self):
        return self._mass
    
    @property
    def surface_tension(self):
        return self._surface_tension
    
    @property
    def hash_value(self):
        return self._hash_value
    
# -------------------------------------- SETTER METHODS ------------------------------------------

    @initial_pos.setter
    def set_initial_pos(self, value):
        self._initial_pos = value

    @colour.setter
    def set_colour(self, value):
        self._colour = value

    @velocity.setter
    def set_velocity(self, value):
        self._velocity = value

    @acceleration.setter
    def set_acceleration(self, value):
        self._acceleration = value

    @shape.setter
    def set_shape(self, value):
        self._shape = value

    @neighbours.setter
    def set_neighbours(self, value):
        self._neighbours = value

    @size.setter
    def set_size(self, value):
        self._size = value

    @mass_density.setter
    def mass_density(self, value):
        self._mass_density = value

    @pressure.setter
    def pressure(self, value):
        self._pressure = value

    @pressure_force.setter
    def pressure_force(self, value):
        self._pressure_force = value

    @gravity.setter
    def gravity(self, value):
        self._gravity = value

    @buoyancy.setter
    def buoyancy(self, value):
        self._buoyancy = value

    @mass.setter
    def mass(self, value):
        self._mass = value

    @surface_tension.setter
    def surface_tension(self, value):
        self._surface_tension = value

    @hash_value.setter
    def hash_value(self, value):
        self._hash_value = value