import numpy as np
import math as m
import random as rd

from typing import Union, Optional

class Particle:

    def __init__(self,
                 init_pos : np.array=np.array([0, 0, 0], dtype="float64"),
                 colour : np.array=np.array([0, 0, 0], dtype="float64"),
                 velocity : np.array=np.array([0.2, 0.1, 0.5], dtype="float64"),
                 acceleration: np.array=np.array([0, 0, 0], dtype="float64"),
                 mass: float=0.1, shape : str="circle", size : int=2,
                 phase_number: int=3) -> None:
        """
            initialize particle properties 
        """
        self.initial_pos = init_pos
        self.predicted_initial_pos = init_pos
        self.neighbour_list = []

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
        self.normal_field = np.array([0, 0, 0], dtype="float64")

        # multi SPH extra attributes

        if phase_number is not None and phase_number>1:
            self.phase_number = phase_number
            self.mass_fractions = [0.1 for i in range(phase_number)]
            self.phase_pressures = [0 for i in range(phase_number)]
            self.phase_stress_tensor = [0 for i in range(phase_number)]
            self.phase_volume_fraction = [1/phase_number for i in range(phase_number)]
            self.phase_viscosities= [0 for i in range(phase_number)]
            self.drift_velocities = [np.array([0, 0, 0], dtype="float64") for i in range(phase_number)]
            self.interp_density = 0
            self.aggregate_phase_viscosity = 0

        # particle force attributes

        self.mass_density = 0
        self.predicted_density = 0
        self.density_change = 0
        self.delta_x = 0
        self.pressure_correction = 0

        self.near_particles = []
        self.far_particles = []

        self.pressure = 0
        self.prev_pressure = 0
        self.displacement = 0
        self.displacement_iter = 0
        self.acceleration_adv = 0

        self.viscosity = np.array([0, 0, 0], dtype="float64")
        self.laminar_viscosity = np.array([0, 0, 0], dtype="float64")
        self.artificial_viscosity = np.array([0, 0, 0], dtype="float64")
        self.pressure_force = np.array([0, 0, 0], dtype="float64")
        self.buoyancy = np.array([0, 0, 0], dtype="float64")
        self.surface_tension = np.array([0, 0 ,0], dtype="float64")
        self.gravity = np.array([0, 0, 0], dtype="float64")
        self.divergence_factor = 0
        self.divergence = np.array([0, 0, 0], dtype="float64")
        self.stiffness_k = np.array([0, 0, 0], dtype="float64")
        self.stiffness_k_v = np.array([0, 0, 0], dtype="float64")

        self.body_force = np.array([0, 0 ,0], dtype="float64")
        self.temperature = 0
        self.thermal_conduc = 0
        self.specific_heat = 0
        self.thermal_diffusion = np.array([0, 0, 0], dtype="float64")

        # update for FSISPH model

        self.cauchy_stress_tensor = np.array([[0, 0, 0], 
                                              [0, 0, 0],
                                              [0, 0, 0]], dtype="float64")
        
        self.deviatoric_stress_tensor = np.array([[0, 0, 0], 
                                                  [0, 0, 0],
                                                  [0, 0, 0]], dtype="float64")
        self.half_velocity = velocity
        self.perp_velocity = velocity
        self.parallel_velocity = velocity
        self.interface_velocity = velocity
        self.interface_normal = np.array([0, 0, 0], dtype="float64")
        self.volume = 0
        self.sound_speed = 0
        self.bulk_modulus = 0
        self.force = np.array([0, 0, 0], dtype="float64")
        
        self.shear_modulus = np.array([[0, 0, 0], 
                                      [0, 0, 0],
                                      [0, 0, 0]], dtype="float64")
        
        self.shear_stress = np.array([[0, 0, 0], 
                                      [0, 0, 0],
                                      [0, 0, 0]], dtype="float64")
        
        self.shear_strain = np.array([[0, 0, 0], 
                                      [0, 0, 0],
                                      [0, 0, 0]], dtype="float64")
        
        # attribs for PBF model

        self.constraint_function = 0
        self.constraint = 0
        self.del_position = np.array([0, 0, 0], dtype="float64")
        self.vorticity = np.array([0, 0, 0], dtype="float64")
        self.vorticity_new = np.array([0, 0, 0], dtype="float64")
        self.vorticity_force = np.array([0, 0, 0], dtype="float64")
        self.vorticity_del = np.array([0, 0, 0], dtype="float64")
        self.stream_function = np.array([0, 0, 0], dtype="float64")


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
