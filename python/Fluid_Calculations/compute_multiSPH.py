import math as m
import numpy as np
import random as rd
import re

from compute_SPH import SPH
from Particles.particles import Particle

class MultiSPH(SPH):

    PHASE_CONSTANTS = {
        "k":0.1,
        "tau":0.5,
        "delta":0.7,
        "lambda_stiffness":7,
        "k_stiffness":7
    }

    def __init__(self, 
                 p: Particle=None,
                 search_method: str=None,
                 hash_table:dict=None,
                 hash_value:int=None,
                 params:dict = None,
                 delta_time:float = None,
                 phase_info:dict = None):
        
        super().__init__(particle=p,
                        search_method=search_method,
                        hash_table=hash_table,
                        hash_value=hash_value,
                        params=params,
                        delta_time=delta_time)

        self.phase_info = phase_info
        self.volume_fractions = [0, 0, 0]
        self.acceleration = np.array([0, 0, 0])

        self.pressure_gradient = [np.array([0, 0, 0]) for i in range(self.particle.phase_number)]
        self.volume_frac_gradient = [np.array([0, 0, 0]) for i in range(self.particle.phase_number)]
        
    def update_interp_density(self):
        for nbr_particle in self.neighbours_list:
            self.particle.interp_density += (
                self.particle.mass * self.kernel_gradient(self.particle.initial_pos - nbr_particle.initial_pos)
            )

    def update_drift_velocity(self):

        self.pressure_grad()
        self.volume_fraction_grad()
        self.update_mass_fraction()

        for i in range(self.particle.phase_number):
            pass

    def mass_frac_density(self):
        for i in range(self.particle.phase_number):
            self.frac_mass_density += (
                self.particle.mass_fractions[i] * self.phase_info["mass_density"][i]
            )

    def mass_frac_density_grad(self):
        for i in range(self.particle.phase_number):
            pass

    def mass_frac_volume_grad(self):
        for i in range(self.particle.phase_number):
            pass

    def update_volume_fraction(self):
        return self.particle.phase_volume_fraction

    def update_convective_acc(self):
        convective_acc = np.array([0, 0, 0])
        for nbr_particle in self.neighbours_list:
            convective_acc += (
                (nbr_particle.mass / nbr_particle.interp_density) *
                self.particle.velocity * self.kernel_gradient(self.particle.initial_pos - nbr_particle.initial_pos)
            )
        self.convective_acc = self.particle.velocity * convective_acc

    def update_gravity(self):
        return super().update_gravity()
    
    def update_mixture_pressure(self):

        self.update_interp_density()
        self.update_mixture_density()

        pressure_const = self.PHASE_CONSTANTS["k_stiffness"]*self.particle.mass_density/self.PHASE_CONSTANTS["lambda_stiffness"]
        density_avg = m.pow(self.particle.interp_density / self.particle.mass_density, self.PHASE_CONSTANTS["lambda_stiffness"])
        self.particle.pressure = pressure_const*(density_avg - 1)

    def update_phase_pressure(self, phase_num:int = None):
        if phase_num is not None:
            self.particle.phase_pressures[phase_num] = self.particle.pressure* \
                                        self.particle.phase_volume_fraction[phase_num]

    def update_mixture_density(self):
        for i in range(self.particle.phase_number):
            self.particle.mass_density += ( 
                self.particle.phase_volume_fraction[i]*
                self.phase_info["mass_density"][i]
            )
    
    def pressure_grad(self, phase_num:int = None):
        pressure_gradient = np.array([0, 0, 0])
        if phase_num is not None:
            for nbr_particle in self.neighbours_list:
                pressure_gradient += (
                    (nbr_particle.mass / nbr_particle.interp_density) *
                    (nbr_particle.phase_pressures[phase_num] - self.particle.phase_pressures[phase_num]) *
                    self.kernel_gradient(self.particle.initial_pos - nbr_particle.initial_pos, 1)
                )
            self.pressure_gradient[phase_num] = pressure_gradient

    def volume_fraction_grad(self, phase_num:int = None):
        volume_fraction_gradient = np.array([0, 0, 0])
        if phase_num is not None:
            for nbr_particle in self.neighbours_list:
                volume_fraction_gradient += (
                    (nbr_particle.mass / nbr_particle.interp_density) *
                    self.particle.phase_volume_fraction[phase_num] - nbr_particle.phase_volume_fraction[phase_num] *
                    self.kernel_gradient(self.particle.initial_pos - nbr_particle.initial_pos, 1)
                )
            self.volume_frac_gradient[phase_num] = volume_fraction_gradient

    def update_mixture_velocity(self):
        intermediate_vel = np.array([0, 0, 0])
        for i in range(self.particle.phase_number):
            intermediate_vel += (
                self.particle.phase_volume_fraction[i]*self.phase_info["mass_density"][i] *
                self.particle.phase_velocities[i]
            )
        self.particle.velocity = (1/self.particle.mass_density) * intermediate_vel

    def update_mass_fraction(self, phase_num:int = None):
        if phase_num is not None:
            self.particle.mass_fractions[phase_num] = self.particle.phase_volume_fraction[phase_num] * \
                            self.phase_info["mass_density"][phase_num] / self.particle.mass_density
            
    def update_acceleration(self):
        
        self.update_gravity()
        self.update_convective_acc()

        self.acceleration = self.particle.gravity - self.convective_acc - self.particle.acceleration

    def update_phase_mass_density(self, phase_num:int = None):
        if phase_num is not None:
            self.particle.phase_mass_densities[phase_num] = self.phase_info["mass_density"][phase_num]* \
                                                        self.particle.phase_volume_fraction[phase_num]

    def update_phase_pressure(self):
        pass


