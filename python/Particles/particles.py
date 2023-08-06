import numpy as np
import math as m
import random as rd

from typing import Union, Optional

class Particle:

    def __init__(self,
                 init_pos : Union(float, float, float) | np.array=np.array([0, 0, 0]),
                 colour : Union(int, int, int, int) | np.array=np.array([0, 0, 0]),
                 velocity : Union(float, float, float) | np.array=np.array([0, 0, 0]),
                 acceleration: Union(float, float, float)| np.array=np.array([0, 0, 0]),
                 mass: float=0.1, shape : str="circle", size : int=2,
                 phase_number: int=3) -> None:
        """
            initialize particle properties 
        """
        self.initial_pos = init_pos
        self.predicted_initial_pos = init_pos

        self.colour = colour

        self.velocity = velocity
        self.predicted_velocity = velocity
        self.advected_velocity = velocity

        self.shape = shape
        self.size = size
        self.acceleration = acceleration
        self.next_acceleration = None
        self.mass = mass
        self.hash_value = 0

        # multi SPH extra attributes

        if phase_number is not None and phase_number>1:
            self.phase_number = phase_number
            self.mass_fractions = [0.1 for i in range(phase_number)]
            self.phase_pressures = [0 for i in range(phase_number)]
            self.phase_velocities = [np.array([0, 0, 0]) for i in range(phase_number)]
            self.phase_stress_tensor = [0 for i in range(phase_number)]
            self.phase_volume_fraction = [0 for i in range(phase_number)]
            self.phase_viscosities= [0 for i in range(phase_number)]
            self.drift_velocities = [np.array([0, 0, 0]) for i in range(phase_number)]
            self.interp_density = 0

        # particle force attributes

        self.mass_density = 0

        self.pressure = 0
        self.prev_pressure = 0
        self.iter_pressure = 0

        self.viscosity = np.array([0, 0, 0])
        self.pressure_force = np.array([0, 0, 0])
        self.gravity = np.array([0, 0, 0])
        self.buoyancy = np.array([0, 0, 0])
        self.surface_tension = np.array([0, 0 ,0])
        self.displacement = np.array([0, 0, 0])
        self.displacement_iter = np.array([0, 0, 0])
        self.acceleration_adv = np.array([0, 0, 0])

# ------------------------------------ PROPERTY GETTERS --------------------------------------
    
    
    def __getattr__(self, name):
        if name.startswith('list_'):
            return []

        raise AttributeError(f"'{self.__class__.__name__}' object has no attribute '{name}'")

    def __setattr__(self, name, value):
        if name.startswith('list_'):
            if not isinstance(value, list):
                raise ValueError("The attribute must be a list.")
            super().__setattr__(name, value)
        else:
            super().__setattr__(name, value)
