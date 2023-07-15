import numpy as np
import math as m
import random as rd

from typing import Union, Optional

class Particles:

    def __init__(self,
                 init_pos : Union(float, float, float) | np.array=None,
                 colour : Union(int, int, int, int)=None,
                 velocity : Union(float, float, float)=None,
                 shape : str=None, size : float=None, neighbours : list=None) -> None:
        """
            initialize particle properties 
        """
        self._initial_pos = init_pos
        self._colour = colour
        self._velocity = velocity
        self._shape = shape
        self._size = size
        self._neighbours = neighbours

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
    def shape(self):
        return self._shape
    
    @property
    def neighbours(self):
        return self._neighbours
    
    @property
    def size(self):
        return self._size
    
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

    @shape.setter
    def set_shape(self, value):
        self._shape = value

    @neighbours.setter
    def set_neighbours(self, value):
        self._neighbours = value

    @size.setter
    def set_size(self, value):
        self._size = value