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
                 delta_time:float = None,
                 phase_info:dict = None):
        
        super().__init__(particle=p,
                         search_method=search_method,
                         hash_table=hash_table,
                         hash_value=hash_value,
                         delta_time=delta_time)

        self.phase_info = phase_info
        self.volume_fractions = [0, 0, 0]
        self.acceleration = np.array([0, 0, 0])

        self.pressure_gradient = [np.array([0, 0, 0]) for i in range(self.particle.phase_number)]
        self.volume_frac_gradient = [np.array([0, 0, 0]) for i in range(self.particle.phase_number)]
        self.frac_mass_density = 0
        self.frac_mass_pressure_grad = np.array([0, 0, 0])
        self.frac_mass_volume_grad = np.array([0, 0, 0])
        self.volume_frac_mix_vel = [np.array([0, 0, 0]) for i in range(self.particle.phase_number)]
        self.mix_vel_volume_frac = [np.array([0, 0, 0]) for i in range(self.particle.phase_number)]
        self.adjusted_pressure = 0

        self.convective_tensor = np.array([0, 0, 0])
        self.viscous_tensor = np.array([0, 0, 0])
        self.pressure_laplacian = np.array([0, 0, 0])

    def update_interp_density(self):
        for nbr_particle in self.neighbours_list:
            self.particle.interp_density += (
                self.particle.mass * self.kernel_gradient(self.particle.initial_pos - nbr_particle.initial_pos)
            )

    def update_drift_velocity(self):

        self.mass_frac_density()
        self.mass_frac_pressure_grad()
        self.mass_frac_volume_grad()
        self.update_acceleration()

        tau = self.PHASE_CONSTANTS["tau"]
        delta = self.PHASE_CONSTANTS["delta"]

        for i in range(self.particle.phase_number):
            tau_density_const = tau * (self.phase_info["mass_density"][i] - self.frac_mass_density) * \
                    self.acceleration
            tau_pressure_const = tau * (self.pressure_gradient[i] - self.frac_mass_pressure_grad)
            delta_volume_frac_const = delta * ((self.volume_frac_gradient[i]/ self.particle.phase_volume_fraction[i]) \
                                               - self.frac_mass_volume_grad)

            self.particle.drift_velocities[i] = (
                tau_density_const - tau_pressure_const - delta_volume_frac_const
            )

    def update_phase_volume_frac(self):
        for i in range(self.particle.phase_number):
            self.volume_frac_grad_mixture_vel(i)
            self.mixture_vel_grad_volume_frac(i)
            
            self.particle.phase_volume_fraction[i] = (
                -1 * self.volume_frac_mix_vel[i] -
                self.mix_vel_volume_frac[i]
            )
            
            if self.volume_frac_detection(self.particle.phase_volume_fraction[i]):
                self.perform_volume_correction(i)
                self.adjust_mixture_pressure()
                   
    def update_convective_m_transfer(self):
        for nbr_particle in self.neighbours_list:
            density_volume_frac_term = 0
            for i in range(self.particle.phase_number):
                volume_frac_neighbour_term = (
                    nbr_particle.phase_volume_fraction[i] * nbr_particle.drift_velocities[i] *
                    (nbr_particle.drift_velocities[i] * self.kernel_gradient(self.particle.initial_pos - 
                                                                             nbr_particle.initial_pos, 1))
                )
                volume_frac_term = (
                    self.particle.phase_volume_fraction[i] * self.particle.drift_velocities[i] *
                    (self.particle.drift_velocities[i] * self.kernel_gradient(self.particle.initial_pos - 
                                                                              nbr_particle.initial_pos, 1))
                )
                density_volume_frac_term += self.phase_info["mass_density"][i] * (volume_frac_neighbour_term +
                                                                                 volume_frac_term)
            self.convective_tensor += (
                nbr_particle.mass / nbr_particle.interp_density *
                density_volume_frac_term
            )
        self.convective_tensor *= -1

    def update_viscous_stress_t(self):

        self.update_aggregate_viscosity()
        for nbr_particle in self.neighbours_list:
            viscous_term = ( nbr_particle.aggregate_phase_viscosity + self.particle.aggregate_phase_viscosity)
            velocity_dif_term = nbr_particle.velocity - self.particle.velocity
            positional_term = (
                (nbr_particle.initial_pos - self.particle.initial_pos) * self.kernel_gradient(self.particle.initial_pos - 
                                                                                              nbr_particle.initial_pos, 1)
            )
            denom_term = np.power(nbr_particle.initial_pos - self.particle.initial_pos, 2)
            position = positional_term / denom_term
            self.viscous_tensor += (
                nbr_particle.mass / nbr_particle.interp_density *
                viscous_term * velocity_dif_term *
                position
            )
    
    def update_aggregate_viscosity(self):
        for i in range(self.particle.phase_number):
            self.particle.aggregate_phase_viscosity += (
                self.particle.phase_volume_fraction[i] * self.phase_info["viscosity"][i]
            )

    def update_pressure_laplacian(self):
        for nbr_particle in self.neighbours_list:
            pressure_term = (self.particle.pressure + nbr_particle.pressure) / 2 * nbr_particle.interp_density
            kernel_grad = self.kernel_gradient(self.particle.initial_pos - nbr_particle.initial_pos)
            self.pressure_laplacian += (
                nbr_particle.mass * pressure_term * kernel_grad
            )

    def update_acceleration_force(self):
        
        self.update_convective_m_transfer()
        self.update_viscous_stress_t()
        self.update_pressure_laplacian()

        pressure_term = self.pressure_laplacian / self.particle.mass_density
        viscous_term = self.viscous_tensor / self.particle.mass_density
        convective_term = self.convective_tensor / self.particle.mass_density

        self.acceleration = -1*pressure_term + self.particle.gravity + viscous_term + convective_term

    def volume_frac_detection(self, volume_fraction:float = None):
        if volume_fraction is not None:
            if volume_fraction < 0:
                self.particle.phase_volume_fraction = 0
                return True
            if volume_fraction > 0:
                return False

    def perform_volume_correction(self, phase_num:int):
        volume_frac_minus_one = self.particle.phase_volume_fraction[phase_num] - 1
        self.particle.phase_volume_fraction[phase_num] -= volume_frac_minus_one

    def adjust_mixture_pressure(self):
        
        for i in range(self.particle.phase_number):
            
            k_stiffness_const = -1 * (self.PHASE_CONSTANTS["k_stiffness"] * self.phase_info["mass_density"][i] /
                                      self.PHASE_CONSTANTS["lambda_stiffness"])
            pressure_const = ((self.PHASE_CONSTANTS["lmabda_stiffness"] - 1) * m.pow((self.particle.interp_density/
                                                                                self.particle.mass_density), self.PHASE_CONSTANTS["lambda_stiffness"])
                                                                                + 1)
            self.adjusted_pressure += (
                k_stiffness_const * pressure_const * self.volume_frac_gradient[i]
            )
        self.particle.pressure += self.adjusted_pressure


    def volume_frac_grad_mixture_vel(self, phase_num:int = None):
        if phase_num is not None:
            for nbr_particle in self.neighbours_list:
                self.volume_frac_mix_vel[phase_num] += (
                    (nbr_particle.mass / nbr_particle.interp_density) *
                    self.particle.phase_volume_fraction[phase_num] + nbr_particle.phase_volume_fraction[phase_num] *
                    (nbr_particle.velocity - self.particle.velocity) *
                    self.kernel_gradient(self.particle.initial_pos - nbr_particle.initial_pos, 1)
                )

    def mixture_vel_grad_volume_frac(self, phase_num:int = None):
        if phase_num is not None:
            for nbr_particle in self.neighbours_list:
                self.mix_vel_volume_frac[phase_num] += (
                    (nbr_particle.mass / nbr_particle.interp_density) *
                    (nbr_particle.phase_volume_fraction[phase_num] * nbr_particle.drift_velocities[phase_num] +
                     self.particle.phase_volume_fraction[phase_num] * self.particle.drift_velocities[phase_num]) *
                     self.kernel_gradient(self.particle.initial_pos - nbr_particle.initial_pos)
                )

    def mass_frac_density(self):
        for i in range(self.particle.phase_number):
            self.update_mass_fraction(i)
            self.frac_mass_density += (
                self.particle.mass_fractions[i] * self.phase_info["mass_density"][i]
            )

    def mass_frac_pressure_grad(self):
        for i in range(self.particle.phase_number):
            self.pressure_grad(i)
            self.frac_mass_pressure_grad += (
                self.pressure_gradient[i] * self.particle.mass_fractions[i]
            )

    def mass_frac_volume_grad(self):
        for i in range(self.particle.phase_number):
            self.volume_fraction_grad(i)
            self.frac_mass_volume_grad += (
                self.particle.mass_fractions[i] * (self.volume_frac_gradient[i] / self.particle.phase_volume_fraction[i])
            )

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
    
    def update(self):
        self.update_acceleration_force()
        self.XSPH_vel_correction()

        self.particle.next_acceleration = self.particle.acceleration

        self.choose_time_stepping(self.time_stepping, self.particle)
        self.choose_collision_types("Cuboid", "Normal")

