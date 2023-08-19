import math as m
import numpy as np
import time

from Fluid_Calculations.compute_sph import SPH
from Particles.particles import Particle

class MultiSPH(SPH):

    PHASE_CONSTANTS = {
        "k":0.1,
        "tau":0.00000008,
        "delta":0.0004,
        "lambda_stiffness":7,
        "k_stiffness":1
    }

    def __init__(self, 
                 particle: Particle=None,
                 search_method: str=None,
                 hash_table:dict=None,
                 hash_value:int=None,
                 all_particles:list = None,
                 time_stepping:str="Euler Cromer",
                 tank_attrs:dict=None,
                 delta_time:float = None,
                 phase_info:dict = None):
        
        super().__init__(particle=particle,
                         search_method=search_method,
                         all_particles=all_particles,
                         hash_table=hash_table,
                         hash_value=hash_value,
                         tank_attrs=tank_attrs,
                         time_stepping=time_stepping,
                         delta_time=delta_time)

        self.phase_info = phase_info
        self.volume_fractions = [0, 0, 0]
        self.acceleration = np.array([0, 0, 0], dtype="float64")

        self.pressure_gradient = [np.array([0, 0, 0], dtype="float64") for i in range(self.particle.phase_number)]
        self.volume_frac_gradient = [np.array([0, 0, 0], dtype="float64") for i in range(self.particle.phase_number)]
        self.frac_mass_density = 0
        self.frac_mass_pressure_grad = np.array([0, 0, 0], dtype="float64")
        self.frac_mass_volume_grad = np.array([0, 0, 0], dtype="float64")
        self.volume_frac_mix_vel = [np.array([0, 0, 0], dtype="float64") for i in range(self.particle.phase_number)]
        self.mix_vel_volume_frac = [np.array([0, 0, 0], dtype="float64") for i in range(self.particle.phase_number)]
        self.adjusted_pressure = 0

        self.convective_tensor = np.array([0, 0, 0], dtype="float64")
        self.viscous_tensor = np.array([0, 0, 0], dtype="float64")
        self.pressure_laplacian = np.array([0, 0, 0], dtype="float64")

    def find_neighbour_list(self, particle):
        particle.neighbour_list = []
        for nbr in self.hash_table[particle.hash_value]:
            if particle is not nbr:
                particle.neighbour_list.append(nbr)

    def normalize(self, vector):
        
        if np.linalg.norm(vector)==0:
            return
        
        return vector / np.linalg.norm(vector)

    # -------------------------------------------------------------- ACCELERATION CALCS --------------------------------------------------------------
                   
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
            try:
                mass_over_interp_density = nbr_particle.mass / nbr_particle.interp_density
            except ZeroDivisionError:
                mass_over_interp_density = 0

            self.convective_tensor += (
                mass_over_interp_density *
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
            try:
                position = positional_term / denom_term
            except ZeroDivisionError:
                position = np.array([0, 0, 0], dtype="float64")
            try:
                mass_over_interp_density = nbr_particle.mass / nbr_particle.interp_density
            except ZeroDivisionError:
                mass_over_interp_density = 0
            self.viscous_tensor += (
                mass_over_interp_density *
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
            try:
                pressure_term = (self.particle.pressure + nbr_particle.pressure) / 2 * nbr_particle.interp_density
            except ZeroDivisionError:
                pressure_term = 0
            kernel_grad = self.kernel_gradient(self.particle.initial_pos - nbr_particle.initial_pos, 1)
            self.pressure_laplacian += (
                nbr_particle.mass * pressure_term * kernel_grad
            )

    def update_acceleration_force(self):
        
        self.update_convective_m_transfer()
        self.update_viscous_stress_t()
        self.update_pressure_laplacian()

        try:
            pressure_term = self.pressure_laplacian / self.particle.mass_density
        except ZeroDivisionError:
            pressure_term = np.array([0, 0, 0], dtype="float64")
        try:
            viscous_term = self.viscous_tensor / self.particle.mass_density
        except ZeroDivisionError:
            viscous_term = np.array([0, 0, 0], dtype="float64")
        try:
            convective_term = self.convective_tensor / self.particle.mass_density
        except ZeroDivisionError:
            convective_term = np.array([0, 0, 0], dtype="float64")


        self.acceleration = -1*pressure_term + self.particle.gravity + viscous_term + convective_term

    # --------------------------------------------------------------- VOLUME CORRECTION ---------------------------------------------------------------
    
    def update_phase_volume_frac(self):
        for i in range(self.particle.phase_number):
            self.volume_frac_grad_mixture_vel(i)
            self.mixture_vel_grad_volume_frac(i)
            
            self.particle.phase_volume_fraction[i] = (
                -self.volume_frac_mix_vel[i] -
                self.mix_vel_volume_frac[i]
            )
            
            if self.volume_frac_detection(self.particle.phase_volume_fraction[i]):
                self.perform_volume_correction(i)
                self.adjust_mixture_pressure()

    def volume_frac_detection(self, volume_fraction:float = None):
        if volume_fraction is not None:
            if np.linalg.norm(volume_fraction) < 0:
                volume_fraction = np.array([0, 0, 0], dtype="float64")
                return False
            if np.linalg.norm(volume_fraction) > 0:
                return True

    def perform_volume_correction(self, phase_num:int):
        self.particle.phase_volume_fraction[phase_num] = self.normalize(self.particle.phase_volume_fraction[phase_num])

    def adjust_mixture_pressure(self):

        self.adjusted_pressure = np.array([0, 0, 0], dtype="float64")
        for i in range(self.particle.phase_number):
            
            k_stiffness_const = -1 * (self.PHASE_CONSTANTS["k_stiffness"] * self.phase_info["mass_density"][i] /
                                      self.PHASE_CONSTANTS["lambda_stiffness"])
            try:
                interp_over_mass_d = self.particle.interp_density/self.particle.mass_density
            except ZeroDivisionError:
                interp_over_mass_d = np.array([0, 0, 0], dtype="float64")

            pressure_const = ((self.PHASE_CONSTANTS["lambda_stiffness"] - 1) * \
                              np.power((interp_over_mass_d), self.PHASE_CONSTANTS["lambda_stiffness"]) + 1)
            self.adjusted_pressure += (
                k_stiffness_const * pressure_const * self.volume_frac_gradient[i]
            )
        self.particle.pressure += self.adjusted_pressure

    # ------------------------------------------------------------ PARTICLE VELOCITY CALCS ------------------------------------------------------------

    def volume_frac_grad_mixture_vel(self, phase_num:int = None):
        if phase_num is not None:
            for nbr_particle in self.neighbours_list:
                try:
                    if any(nbr_particle.interp_density)==0:
                        nbr_mass_over_interp = np.array([0, 0, 0], dtype="float64")
                    else:
                        nbr_mass_over_interp = nbr_particle.mass / nbr_particle.interp_density
                except ZeroDivisionError:
                    nbr_mass_over_interp = np.array([0, 0, 0], dtype="float64")
                
                if any(nbr_mass_over_interp)==0:
                    self.volume_frac_mix_vel[phase_num] += np.array([0, 0, 0], dtype="float64")
                else:
                    self.volume_frac_mix_vel[phase_num] += (
                        nbr_mass_over_interp *
                        (self.particle.phase_volume_fraction[phase_num] + nbr_particle.phase_volume_fraction[phase_num]) / 2 *
                        (nbr_particle.velocity - self.particle.velocity) *
                        self.kernel_gradient(self.particle.initial_pos - nbr_particle.initial_pos, 1)
                    )

    def mixture_vel_grad_volume_frac(self, phase_num:int = None):
        if phase_num is not None:
            for nbr_particle in self.neighbours_list:
                try:
                    if any(nbr_particle.interp_density)==0:
                        nbr_mass_over_interp = np.array([0, 0, 0], dtype="float64")
                    else:
                        nbr_mass_over_interp = nbr_particle.mass / nbr_particle.interp_density
                except ZeroDivisionError:
                    nbr_mass_over_interp = np.array([0, 0, 0], dtype="float64")
                
                if any(nbr_mass_over_interp)==0:
                    self.volume_frac_mix_vel[phase_num] += np.array([0, 0, 0], dtype="float64")
                else:
                    self.mix_vel_volume_frac[phase_num] += (
                        (nbr_mass_over_interp) *
                        (nbr_particle.phase_volume_fraction[phase_num] * nbr_particle.drift_velocities[phase_num] +
                        self.particle.phase_volume_fraction[phase_num] * self.particle.drift_velocities[phase_num]) *
                        self.kernel_gradient(self.particle.initial_pos - nbr_particle.initial_pos, 1)
                    )

    # ---------------------------------------------------------- DRIFT VELOCITY CALCULATIONS ----------------------------------------------------------

    def update_mass_density(self):
        self.particle.mass_density = 0
        for i in range(self.particle.phase_number): 
            self.particle.mass_density += self.phase_info["mass_density"][i]

    def update_gravity(self):
        return super().update_gravity()
    
    def update_mixture_pressure(self, particle):

        particle.pressure = np.array([0, 0, 0], dtype="float64")
        try:
            interp_over_mass_d = particle.interp_density / particle.mass_density
        except ZeroDivisionError:
            interp_over_mass_d = np.array([0, 0, 0], dtype="float64")

        pressure_const = self.PHASE_CONSTANTS["k_stiffness"]*particle.mass_density/self.PHASE_CONSTANTS["lambda_stiffness"]
        density_avg = np.power(interp_over_mass_d, self.PHASE_CONSTANTS["lambda_stiffness"])
        particle.pressure = pressure_const*(density_avg - 1)

    def update_phase_pressure(self, phase_num:int = None, particle:Particle=None):
        self.update_mixture_pressure(particle)
        particle.phase_pressures[phase_num] = particle.pressure*particle.phase_volume_fraction[phase_num]

    def update_interp_density(self, particle):
        particle.interp_density = np.array([0, 0, 0], dtype="float64")
        for nbr_particle in particle.neighbour_list:
            particle.interp_density += (
                particle.mass * self.kernel_gradient(particle.initial_pos - nbr_particle.initial_pos, 0)
            )

    def update_mixture_density(self, particle):
        for i in range(self.particle.phase_number):
            particle.mass_density += ( 
                particle.phase_volume_fraction[i]*
                self.phase_info["mass_density"][i]
            )

    def update_mixture_velocity(self):
        intermediate_vel = np.array([0, 0, 0], dtype="float64")
        for i in range(self.particle.phase_number):
            intermediate_vel += (
                self.particle.phase_volume_fraction[i]*self.phase_info["mass_density"][i] *
                self.particle.phase_velocities[i]
            )
        self.particle.velocity = (1/self.particle.mass_density) * intermediate_vel
    
    def update_acceleration(self):
        
        self.update_gravity()
        self.acceleration = self.particle.gravity - self.particle.acceleration

    def update_mass_fraction(self, phase_num:int = None):
        self.particle.mass_fractions[phase_num] = self.particle.phase_volume_fraction[phase_num] * \
                            self.phase_info["mass_density"][phase_num] / self.particle.mass_density
            
    def pressure_grad(self, phase_num:int = None):

        pressure_gradient = np.array([0, 0, 0], dtype="float64")
        if phase_num is not None:
            for nbr_particle in self.neighbours_list:

                self.update_phase_pressure(phase_num, nbr_particle)
                self.update_phase_pressure(phase_num, self.particle)
                
                try:
                    if any(nbr_particle.interp_density == 0):
                        nbr_mass_over_interp = np.array([0, 0, 0], dtype="float64")
                    else:
                        nbr_mass_over_interp = nbr_particle.mass / nbr_particle.interp_density
                except ZeroDivisionError:
                    nbr_mass_over_interp = np.array([0, 0, 0], dtype="float64")

                if any(nbr_mass_over_interp) == 0:
                    pressure_gradient += np.array([0, 0, 0], dtype="float64")
                else:
                    pressure_gradient += (
                    nbr_mass_over_interp *
                    (nbr_particle.phase_pressures[phase_num] - self.particle.phase_pressures[phase_num]) *
                    self.kernel_gradient(self.particle.initial_pos - nbr_particle.initial_pos, 1)
                    )

            self.pressure_gradient[phase_num] = pressure_gradient

    def volume_fraction_grad(self, phase_num:int = None):
        volume_fraction_gradient = np.array([0, 0, 0], dtype="float64")
        if phase_num is not None:
            for nbr_particle in self.neighbours_list:

                try:
                    if any(nbr_particle.interp_density == 0):
                        nbr_mass_over_interp = np.array([0, 0, 0], dtype="float64")
                    else:
                        nbr_mass_over_interp = nbr_particle.mass / nbr_particle.interp_density
                except ZeroDivisionError:
                    nbr_mass_over_interp = np.array([0, 0, 0], dtype="float64")

                if any(nbr_mass_over_interp)==0:
                    volume_fraction_gradient += np.array([0, 0, 0], dtype="float64")
                else:
                    volume_fraction_gradient += (
                        nbr_mass_over_interp * self.particle.phase_volume_fraction[phase_num] - nbr_particle.phase_volume_fraction[phase_num] *
                        self.kernel_gradient(self.particle.initial_pos - nbr_particle.initial_pos, 1)
                    )
            self.volume_frac_gradient[phase_num] = volume_fraction_gradient

    def mass_frac_pressure_grad(self):
        for i in range(self.particle.phase_number):
            self.pressure_grad(i)

            self.frac_mass_pressure_grad += (
                self.pressure_gradient[i] * self.particle.mass_fractions[i]
            )

    def mass_frac_volume_grad(self):
        for i in range(self.particle.phase_number):
            self.volume_fraction_grad(i)

            volume_grad = self.volume_frac_gradient[i] / self.particle.phase_volume_fraction[i]
            if any(volume_grad)==0:
                self.frac_mass_volume_grad += np.array([0, 0, 0], dtype="float64")
            else:
                self.frac_mass_volume_grad += (
                    self.particle.mass_fractions[i] * volume_grad
                )

    def mass_frac_density(self):
        for i in range(self.particle.phase_number):
            self.update_mass_fraction(i)
            self.frac_mass_density += (
                self.particle.mass_fractions[i] * self.phase_info["mass_density"][i]
            )

    def update_drift_velocity(self):

        self.update_mass_density()
        self.mass_frac_density()
        self.mass_frac_pressure_grad()
        self.mass_frac_volume_grad()
        self.update_acceleration()

        tau = self.PHASE_CONSTANTS["tau"]
        delta = self.PHASE_CONSTANTS["delta"]

        for i in range(self.particle.phase_number):
            tau_density_const = tau * (self.phase_info["mass_density"][i] - self.frac_mass_density) * self.acceleration
            tau_pressure_const = tau * (self.pressure_gradient[i] - self.frac_mass_pressure_grad)
            try:
                if any(self.volume_frac_gradient[i])==0:
                    volume_frac_grad_phase = np.array([0, 0, 0], dtype="float64")
                else:
                    volume_frac_grad_phase = self.volume_frac_gradient[i] / self.particle.phase_volume_fraction[i]
            except ZeroDivisionError:
                volume_frac_grad_phase = np.array([0, 0, 0], dtype="float64")

            delta_volume_frac_const = delta * (volume_frac_grad_phase- self.frac_mass_volume_grad)

            self.particle.drift_velocities[i] = (
                tau_density_const - tau_pressure_const - delta_volume_frac_const
            )

    def debugging_attribs(self, secs):

        print("Drift velocities are: ", self.particle.drift_velocities)
        print("Mass fraction density is", self.frac_mass_density)
        print("Acceleration is:", self.acceleration)
        print("Phase volume fractions are: ", self.particle.phase_volume_fraction)
        print("Mixture pressure is: ", self.particle.pressure)
        print("Phase pressures are:", self.particle.phase_pressures)
        print("Mass density is: ", self.particle.mass_density)
        print("Interpolated density is:", self.particle.interp_density)
        print("Mass fractions are:", self.particle.mass_fractions)
        print("\n\n")
        time.sleep(secs)
    # ------------------------------------------------------------------- TIMESTEP --------------------------------------------------------------------

    def perform_calculations(self):

        for nbr in self.neighbours_list:
            self.find_neighbour_list(nbr)
            self.update_interp_density(nbr)

        self.find_neighbour_list(self.particle)
        self.update_interp_density(self.particle)
        self.update_mixture_density(self.particle)

        self.update_drift_velocity()
        self.update_phase_volume_frac()
        """ self.update_acceleration_force() """

    def update(self):

        self.perform_calculations()

        self.debugging_attribs(0.1)
        self.XSPH_vel_correction()

        self.choose_time_stepping(self.time_stepping)
        self.choose_collision_types("Cuboid", "Normal")

