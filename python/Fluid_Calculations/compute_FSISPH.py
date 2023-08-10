import math as m
import numpy as np
import random as rd
import re

from compute_SPH import SPH
from Particles.particles import Particle

class FSISPH(SPH):

    OTHER_PARAMS = {
        "alpha_const":0.25,
        "beta_const":0.5,
        "lambda_const":7
    }

    def __init__(self,
                particle: Particle=None,
                search_method:str = None,
                hash_table:dict = None,
                hash_value:int = None,
                time_stepping:str = "Euler Cromer",
                tank_attrs:dict = None,
                delta_time:int = 0.02,
                slip_condition:bool = True):
        
        super().__init__(particle=particle,
                         search_method=search_method,
                         hash_table=hash_table,
                         hash_value=hash_value,
                         time_stepping=time_stepping,
                         tank_attrs=tank_attrs,
                         delta_time=delta_time)

        self.velocity_grad_tensor = np.array([[0, 0, 0],
                                             [0, 0, 0],
                                             [0, 0, 0]], dtype="float64")
        
        self.interface_velocity = np.array([5.2, 0.1, 3.2], dtype="float64")
        self.force = np.array([0, 0, 0], dtype="float64")

        self.interface_neighbrs = []
        self.slip_condition = slip_condition

        self.trapezoidal_sub_interval_amt = 100
        self.delta_x = 0.01

    def deformation(self):
        return 1/2 * (self.velocity_grad_tensor - self.velocity_grad_tensor.transpose())
    
    def rotation(self):
        return 1/2 * (self.velocity_grad_tensor + self.velocity_grad_tensor.transpose())
    
    def update_velocity_grad(self):
        self.velocity_grad = np.array([0, 0, 0], dtype="float64")
        for nbr_particle in self.neighbours_list:
            self.velocity_grad += (
                nbr_particle.mass  * (
                self.particle.velocity - nbr_particle.velocity
                )
                * self.cubic_spline_kernel(kernel_type=1,
                                           nbr_position=nbr_particle.initial_pos)
            )
        return self.velocity_grad

    def update_velocity_grad_c(self, index_c, index_kernel):
        self.velocity_grad = np.array([0, 0, 0], dtype="float64")
        for nbr_particle in self.neighbours_list:
            self.velocity_grad[index_c] += (
                nbr_particle.mass  * (
                self.particle.velocity[index_c] - nbr_particle.velocity[index_c]
                )
                * self.cubic_spline_kernel_c(self.particle.initial_pos - nbr_particle.initial_pos, index_kernel)
            )
        return self.velocity_grad[index_c]
    
    def update_velocity_grad_t(self):
        for i in range(np.shape(self.velocity_grad_tensor)[0]):
            for j in range(np.shape(self.velocity_grad_tensor)[1]):
                self.velocity_grad_tensor[i][j] = (
                    self.update_velocity_grad_c(i, j)
                )
    
    @staticmethod
    def kronecker_delta(i, j):
        if (i==j): return 1 
        else: return 0

    @staticmethod
    def reverse_kronecker_delta(i, j):
        if (i==j): return 0
        else: return 1

    def deviatoric_stress(self):
        for i in range(np.shape(self.particle.deviatoric_stress_tensor)[0]):
            for j in range(np.shape(self.particle.deviatoric_stress_tensor)[1]):
                self.particle.deviatoric_stress_tensor[i][j] = (
                    self.particle.cauchy_stress_tensor[i][j] - 
                    self.isotropic_pressure_component() * self.kronecker_delta(i, j)
                )
        return self.particle.deviatoric_stress_tensor
    
    def cauchy_stress(self):
        for i in range(np.shape(self.particle.cauchy_stress_tensor)[0]):
            for j in range(np.shape(self.particle.cauchy_stress_tensor)[1]):
                self.particle.cauchy_stress_tensor[i][j] = (
                    -self.particle.pressure * self.kronecker_delta(i, j) +
                    (self.update_velocity_grad(i, j) + self.update_velocity_grad_c(j, i)) *
                    self.PARAMETERS["viscosity"]
                )
        return self.particle.cauchy_stress_tensor

    def isotropic_pressure_component(self):
        return (self.particle.cauchy_stress_tensor.trace()) / 3
    
    def update_deviatoric_stress(self):
        self.update_velocity_grad_t()
        for i in range(np.shape(self.particle.deviatoric_stress_tensor)[0]):
            for j in range(np.shape(self.particle.deviatoric_stress_tensor)[1]):
                self.particle.deviatoric_stress_tensor[i][j] = (
                    2 * self.PARAMETERS["viscosity"] *
                    (self.deformation()[i][j] - self.deformation().trace()*np.identity(3)) +
                    self.deviatoric_stress()[i][j]*self.rotation().transpose()[i][j] +
                    self.rotation()[i][j]*self.deviatoric_stress().transpose()[i][j]
                )

    def update_cauchy_stress(self):

        self.update_velocity_grad_t()
        self.update_deviatoric_stress()

        for i in range(np.shape(self.particle.cauchy_stress_tensor)[0]):
            for j in range(np.shape(self.particle.cauchy_stress_tensor)[1]):
                self.particle.cauchy_stress_tensor[i][j] = (
                    self.particle.deviatoric_stress_tensor[i][j] - 
                    self.particle.pressure*np.identity(3)
                )

    def update_acceleration_force(self):
        pass

    def update_shear_modulus(self, id):
        self.particle.shear_modulus = self.update_shear_stress() / self.update_shear_strain(id)

    def update_bulk_modulus(self):
        self.particle.bulk_modulus = (
            self.particle.pressure / (
            m.pow((self.particle.mass_density/ self.PARAMETERS["mass_density"]), self.OTHER_PARAMS["lambda_const"]) - 1
            )
        )

    def update_shear_stress(self):
        for i in range(np.shape(self.particle.cauchy_stress_tensor)[0]):
            for j in range(np.shape(self.particle.cauchy_stress_tensor)[1]):
                self.particle.shear_stress[i][j] = (
                    self.particle.cauchy_stress_tensor[i][j] *
                    self.reverse_kronecker_delta(i, j)
                )
        return self.particle.shear_stress
    
    def update_shear_strain(self, id):
        for i in range(np.shape(self.particle.shear_strain)[0]):
            for j in range(np.shape(self.particle.shear_strain)[1]):
                self.particle.shear_strain[i][j] = (
                    (self.neighbours_list[id].initial_pos[j] - self.particle.initial_pos[j]) / 
                    (self.neighbours_list[id].initial_pos[i] - self.particle.initial_pos[i]) 
                )
        return self.particle.shear_strain

    def update_perpendicular_velocity(self, id):
        self.particle.perp_velocity = self.particle.velocity - \
            self.particle.velocity*self.normalize(self.particle.initial_pos - self.neighbours_list[id].initial_pos)
        
        self.update_shear_modulus(id)
        self.local_domain_volume(id)

        numerator = (self.particle.shear_modulus*self.neighbours_list[id].volume* \
                    self.cubic_spline_kernel_gradient_mag(self.particle.initial_pos - 
                                                          self.neighbours_list[id].initial_pos)* \
                    self.particle.perp_velocity) + \
                    (self.neighbours_list[id].shear_modulus*self.particle.volume* \
                    self.cubic_spline_kernel_pos(self.neighbours_list[id].initial_pos - 
                                                 self.particle.initial_pos)* \
                    self.neighbours_list[id].perp_velocity)
        
        denom = (self.particle.shear_modulus*self.neighbours_list[id].volume* \
                self.cubic_spline_kernel_gradient_mag(self.particle.initial_pos - 
                                                      self.neighbours_list[id].initial_pos)) + \
                (self.neighbours_list[id].shear_modulus*self.particle.volume* \
                self.cubic_spline_kernel_pos(self.neighbours_list[id].initial_pos - 
                                             self.particle.initial_pos))
        
        self.particle.perp_velocity = (
            numerator / denom
        )


    def local_domain_volume(self, id):
        trapezoidal_amt = self.trapezoidal_sub_interval_amt
        while (self.trapezoidal_sub_interval_amt > 0):
            
            if self.trapezoidal_sub_interval_amt == trapezoidal_amt:
                self.particle.volume += self.cubic_spline_kernel_pos(self.particle.initial_pos - 
                                                                     self.neighbours_list[id].initial_pos)

            self.particle.volume += 2*self.cubic_spline_kernel_pos(self.particle.initial_pos - 
                                                                     self.neighbours_list[id].initial_pos)

            self.trapezoidal_sub_interval_amt -= 1
        self.particle.volume *= 0.5*self.delta_x

    def update_parallel_velocity(self, id):
        self.particle.parallel_velocity = self.particle.velocity* \
            self.normalize(self.particle.initial_pos - self.neighbours_list[id].initial_pos)
        
        self.update_bulk_modulus()
        self.local_domain_volume(id)

        numerator = (
            (self.particle.bulk_modulus*self.neighbours_list[id].volume*self.cubic_spline_kernel_gradient_mag(
            self.particle.initial_pos - self.neighbours_list[id].initial_pos) * self.particle.parallel_velocity) +
            (self.neighbours_list[id].bulk_modulus*self.neighbours_list[id].volume*self.cubic_spline_kernel_pos(
                self.neighbours_list[id] - self.particle.initial_pos) * self.neighbours_list[id].parallel_velocity)
        )

        denominator = (
            (self.particle.bulk_modulus*self.neighbours_list[id].volume*self.cubic_spline_kernel_gradient_mag(
            self.particle.initial_pos - self.neighbours_list[id].initial_pos)) + 
            (self.neighbours_list[id].bulk_modulus*self.particle.volume*self.cubic_spline_kernel_gradient_mag(
            self.neighbours_list[id].initial_pos - self.particle.initial_pos))  
        )

        self.particle.parallel_velocity = (
            numerator / denominator
        )
        
    def update_slip_condition(self, id):
        self.interface_velocity = (
                self.particle.parallel_velocity*self.normalize(self.particle.initial_pos) + 
                self.particle.perp_velocity
            )
        self.update_stress_state()

    def update_no_slip_condition(self, id):
        self.interface_velocity = (
                self.particle.parallel_velocity*self.normalize(self.particle.initial_pos - 
                                                               self.neighbours_list[id].initial_pos) + 
                self.particle.perp_velocity
            )
        self.update_stress_state()

    def update_stress_state(self):
        if self.slip_condition is True:
            numerator = (
            (self.outward_facing_unit_norm(self.particle.initial_pos -
                                          self.neighbours_list[id].initial_pos) * self.particle.cauchy_stress_tensor * 
            self.outward_facing_unit_norm(self.particle.initial_pos -
                                          self.neighbours_list[id].initial_pos) * self.neighbours_list[id].mass_density) +
            (self.outward_facing_unit_norm(self.neighbours_list[id].initial_pos -
                                          self.particle.initial_pos) * self.neighbours_list[id].cauchy_stress_tensor *
                                          self.neighbours_list[id].mass_density)
            )

            denominator = ( self.particle.mass_density + self.neighbours_list[id].mass_density)

            self.particle.cauchy_stress_tensor = (
                np.identity(3) * (numerator / denominator)
            )
        if self.slip_condition is not True:
            numerator = (
            (self.particle.cauchy_stress_tensor * self.neighbours_list[id].mass_density) +
              (self.neighbours_list[id].cauchy_stress_tensor*self.particle.mass_density) 
            )

            denominator = ( self.particle.mass_density + self.neighbours_list[id].mass_density)

            self.particle.cauchy_stress_tensor = numerator / denominator

    def update_conservative_diffusion(self):
        pass

    def update_interface_neighbours(self):
        pass
    
    def outward_facing_unit_norm(self, position):
        return self.normalize(position)
    
    def update_interface_normals(self):
        pass

    def update_artificial_viscosity_slip(self):
        pass
    
    def update_sound_speed(self):
        pass
    
    def update_artificial_viscosity_no_slip(self, id):
        viscosity_term = (
            self.particle.velocity - self.neighbours_list[id].velocity *
            self.particle.initial_pos - self.neighbours_list[id].initial_pos *
            self.PARAMETERS["cell_size"] / 
            m.pow(np.linalg.norm(
            self.particle.initial_pos - self.neighbours_list[id].initial_pos
            ), 2)
        )
        viscosity = max(viscosity_term, 0)
        const_term = -self.force / (self.particle.mass_density - self.neighbours_list[id].mass_density)
        self.particle.artificial_viscosity = (
            const_term * 
            (self.OTHER_PARAMS["alpha_const"])
        )

    def update_force(self, id):
        self.force = (
            (self.particle.acceleration + self.neighbours_list[id].acceleration) *
            (self.particle.half_velocity - self.neighbours_list[id].half_velocity) / 
            (2*self.particle.acceleration * (
                self.particle.half_velocity - self.interface_velocity
            ) - 2*self.neighbours_list[id].acceleration*(
                self.interface_velocity - self.neighbours_list[id].half_velocity
            ))
        )
        return self.force

    def update_half_velocity(self):
        self.particle.half_velocity = self.particle.velocity + self.particle.acceleration*self.delta_time/2

    def update_mass_density(self):
        mass_density = np.array([0, 0, 0], dtype="float64")
        for nbr_particle in self.neighbours_list:
            mass_density += (
                nbr_particle.mass / nbr_particle.mass_density *
                (self.particle.half_velocity - self.interface_velocity) *
                self.cubic_spline_kernel(kernel_type=1,
                                         nbr_position=nbr_particle.initial_pos)
            )
        self.particle.mass_density *= -2*mass_density

    @staticmethod
    def normalize(position):
        return position / np.linalg.norm(position)
    
    def update_interface_velocity(self):
        
        self.update_parallel_velocity()
        self.update_perpendicular_velocity()
        
        if self.slip_condition is True:
            self.update_slip_condition()
        elif self.slip_condition is not True:
            self.update_no_slip_condition()
        else:
            self.interface_velocity = self.particle.velocity/ np.sum([nbr_particle.velocity for nbr_particle in self.neighbours_list])

    def update_velocity(self):
        velocity = np.array([0, 0, 0], dtype="float64")
        for id, nbr_particle in enumerate(self.neighbours_list):
            velocity += (
                self.update_force(id) * nbr_particle.mass *
                self.particle.acceleration *
                (self.particle.half_velocity - self.interface_velocity)
            )
        self.particle.velocity = velocity*2

