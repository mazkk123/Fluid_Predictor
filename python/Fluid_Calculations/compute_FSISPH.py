import math as m
import numpy as np
import time

from Fluid_Calculations.compute_sph import SPH
from Particles.particles import Particle

class FSISPH(SPH):

    OTHER_PARAMS = {
        "alpha_const":0.25,
        "beta_const":0.5,
        "lambda_const":7,
        "deformation":0.7,
        "diffusion_coefficient":0.05
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
        
        self.force = np.array([[0, 0, 0],
                               [0, 0, 0],
                               [0, 0, 0]], dtype="float64")

        self.interface_neighbrs = []
        self.slip_condition = slip_condition

    # ---------------------------------------------------------------- TENSOR ATTRIBUTES ----------------------------------------------------------------------

    def deformation(self):
        return 1/2 * (self.velocity_grad_tensor - self.velocity_grad_tensor.transpose())
    
    def rotation(self):
        return 1/2 * (self.velocity_grad_tensor + self.velocity_grad_tensor.transpose())
    
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

    def update_deviatoric_stress(self):
        self.update_velocity_grad_t()
        trace_matrix = self.deformation().trace()*np.identity(3)
        for i in range(np.shape(self.particle.deviatoric_stress_tensor)[0]):
            for j in range(np.shape(self.particle.deviatoric_stress_tensor)[1]):
                self.particle.deviatoric_stress_tensor[i][j] = (
                    2 * self.PARAMETERS["viscosity"] *
                    (self.deformation()[i][j] - trace_matrix[i][j]) +
                    self.deviatoric_stress()[i][j]*self.rotation().transpose()[i][j] +
                    self.rotation()[i][j]*self.deviatoric_stress().transpose()[i][j]
                )

    def update_cauchy_stress(self):

        self.update_velocity_grad_t()
        self.update_deviatoric_stress()

        trace_pressure_matrix = self.particle.pressure*np.identity(3)
        for i in range(np.shape(self.particle.cauchy_stress_tensor)[0]):
            for j in range(np.shape(self.particle.cauchy_stress_tensor)[1]):
                self.particle.cauchy_stress_tensor[i][j] = (
                    self.particle.deviatoric_stress_tensor[i][j] - 
                    trace_pressure_matrix[i][j]
                )

    def isotropic_pressure_component(self):
        return (self.particle.cauchy_stress_tensor.trace()) / 3

    # ------------------------------------------------------------- MATHEMATICAL CALCULATIONS -----------------------------------------------------------------

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

    def update_velocity_grad_c(self, index_c, kernel_index):
        self.velocity_grad = np.array([0, 0, 0], dtype="float64")
        for nbr_particle in self.neighbours_list:
            self.velocity_grad[index_c] += (
                nbr_particle.mass  * (
                self.particle.velocity[index_c] - nbr_particle.velocity[index_c]
                )
                * self.cubic_spline_kernel_grad_c(self.particle.initial_pos - nbr_particle.initial_pos, kernel_index)
            )
        return self.velocity_grad[index_c]
    
    def update_velocity_grad_t(self):
        for i in range(np.shape(self.velocity_grad_tensor)[0]):
            for j in range(np.shape(self.velocity_grad_tensor)[1]):
                self.velocity_grad_tensor[i][j] = (
                    self.update_velocity_grad_c(i, j)
                )

    def update_acceleration_force(self):
        
        if self.slip_condition is True:
            self.update_artificial_viscosity_slip()
        elif self.slip_condition is not True:
            self.update_artificial_viscosity_no_slip()
        
        for nbr_particle in self.neighbours_list:
            try:
                divisor = self.particle.cauchy_stress_tensor / (self.particle.mass_density*nbr_particle.mass_density)
            except ZeroDivisionError:
                divisor = np.array([[0, 0, 0],
                                    [0, 0, 0],
                                    [0, 0, 0]], dtype="float64")
            tensor_term = ( 
            divisor - (self.particle.artificial_viscosity / 2)
            )
            
            kernel_term = (
                    self.new_cubic_spline_kernel_gradient(
                    np.abs( self.particle.initial_pos -  nbr_particle.initial_pos)
                    )
            )
        
            self.particle.acceleration += np.dot(tensor_term, kernel_term)

    def update_shear_modulus(self):

        self.update_shear_stress()
        for id, nbr_particle in enumerate(self.neighbours_list):
            self.particle.shear_modulus += self.particle.shear_stress / self.update_shear_strain(id)

    def update_bulk_modulus(self):
        try:
            const_term = np.power((self.particle.mass_density/ self.PARAMETERS["mass_density"]), self.OTHER_PARAMS["lambda_const"]) - 1
        except ZeroDivisionError:
            const_term = 0
        self.particle.bulk_modulus = (
            self.particle.pressure / const_term
        )

    # ---------------------------------------------------------------- FORCE CALCULATIONS -------------------------------------------------------------------------

    def update_perpendicular_velocity(self):

        intermediate_perp_vel = np.array([0, 0, 0], dtype="float64")
        for nbr_particle in self.neighbours_list:
            intermediate_perp_vel = self.particle.velocity - \
                self.particle.velocity*self.normalize(self.particle.initial_pos - nbr_particle.initial_pos)
        
            self.update_shear_modulus()
            self.update_volume()

            numerator = (self.particle.shear_modulus*nbr_particle.volume* \
                        self.cubic_spline_kernel_gradient_mag(self.particle.initial_pos - 
                                                            nbr_particle.initial_pos)* \
                        intermediate_perp_vel) + \
                        (nbr_particle.shear_modulus*self.particle.volume* \
                        self.cubic_spline_kernel_pos(nbr_particle.initial_pos - 
                                                    self.particle.initial_pos)* \
                        nbr_particle.perp_velocity)
            
            denom = (self.particle.shear_modulus*nbr_particle.volume* \
                    self.cubic_spline_kernel_gradient_mag(self.particle.initial_pos - 
                                                        nbr_particle.initial_pos)) + \
                    (nbr_particle.shear_modulus*self.particle.volume* \
                    self.cubic_spline_kernel_pos(nbr_particle.initial_pos - 
                                                self.particle.initial_pos))

            if any(denom[0])==0:
                perp_velocity = np.array([[0, 0, 0],
                                          [0, 0, 0]
                                          [0, 0, 0]], dtype="float64")
            else:
                perp_velocity = numerator / denom
            
            self.particle.perp_velocity += perp_velocity

    def update_volume(self):
        try:
            self.particle.volume = self.particle.mass / self.particle.mass_density
        except ZeroDivisionError:
            self.particle.volume = 0

    def update_parallel_velocity(self):
        
        self.update_bulk_modulus()
        self.update_volume()

        for id, nbr_particle in enumerate(self.neighbours_list):
            self.particle.parallel_velocity = self.particle.velocity* \
                self.normalize(self.particle.initial_pos - nbr_particle.initial_pos)

            numerator = (
                (self.particle.bulk_modulus*nbr_particle.volume*self.cubic_spline_kernel_gradient_mag(
                self.particle.initial_pos - nbr_particle.initial_pos) * self.particle.parallel_velocity) +
                (nbr_particle.bulk_modulus*nbr_particle.volume*self.cubic_spline_kernel_pos(
                 nbr_particle.initial_pos - self.particle.initial_pos) * nbr_particle.parallel_velocity)
            )

            denominator = (
                (self.particle.bulk_modulus*nbr_particle.volume*self.cubic_spline_kernel_gradient_mag(
                self.particle.initial_pos - nbr_particle.initial_pos)) + 
                (nbr_particle.bulk_modulus*self.particle.volume*self.cubic_spline_kernel_gradient_mag(
                nbr_particle.initial_pos - self.particle.initial_pos))  
            )

            if any(denominator)==0:
                parallel_velocity = np.array([0, 0, 0], dtype="float64")
            else:
                parallel_velocity = numerator / denominator

            self.particle.parallel_velocity += parallel_velocity
        
    def update_slip_condition(self):
        
        self.update_stress_state()

        self.interface_velocity = (
                self.particle.parallel_velocity*self.normalize(self.particle.initial_pos) + 
                self.particle.perp_velocity
            )

    def update_no_slip_condition(self):

        self.update_stress_state()

        for nbr_particle in self.neighbours_list:
            self.interface_velocity = (
                    self.particle.parallel_velocity*self.normalize(self.particle.initial_pos - 
                                                                nbr_particle.initial_pos) + 
                    self.particle.perp_velocity
            )

    def update_stress_state(self):
        for nbr_particle in self.neighbours_list:
            if self.slip_condition is True:
                numerator = (
                (self.outward_facing_unit_norm(self.particle.initial_pos -
                                            nbr_particle.initial_pos) * self.particle.cauchy_stress_tensor * 
                self.outward_facing_unit_norm(self.particle.initial_pos -
                                            nbr_particle.initial_pos) * nbr_particle.mass_density) +
                (self.outward_facing_unit_norm(nbr_particle.initial_pos -
                                            self.particle.initial_pos) * nbr_particle.cauchy_stress_tensor *
                                            nbr_particle.mass_density)
                )

                denominator = self.particle.mass_density + nbr_particle.mass_density

                if all(denominator)==0:
                    cauchy_stress = np.array([0, 0, 0], dtype="float64")
                else:
                    cauchy_stress = numerator / denominator

                self.particle.cauchy_stress_tensor = (
                    np.identity(3) * cauchy_stress
                )

            if self.slip_condition is not True:
                numerator = (
                (self.particle.cauchy_stress_tensor * nbr_particle.mass_density) +
                (nbr_particle.cauchy_stress_tensor*self.particle.mass_density) 
                )

                if denominator==0:
                    cauchy_stress = np.array([0, 0, 0], dtype="float64")
                else:
                    cauchy_stress = numerator / denominator

                self.particle.cauchy_stress_tensor += cauchy_stress

    # --------------------------------------------------------------- CONSERVATIVE DIFFUSION ---------------------------------------------------------------------

    def update_conservative_diffusion_mass_density(self):
        mass_density = 0
        for id, nbr_particle in enumerate(self.neighbours_list):
            if self.particle.sound_speed == nbr_particle.sound_speed:
                numerator = (
                     (nbr_particle.mass_density - self.particle.mass_density) *
                     (self.particle.sound_speed - nbr_particle.sound_speed) *
                     (self.PARAMETERS["cell_size"]*(self.particle.initial_pos -
                     nbr_particle.initial_pos) * self.cubic_spline_kernel_gradient(
                      self.particle.initial_pos - nbr_particle.initial_pos
                    ))
                )
                denominator = (
                   nbr_particle.mass_density*(m.pow(np.linalg.norm(
                    self.particle.initial_pos - nbr_particle.initial_pos), 2) +
                    self.OTHER_PARAMS["deformation"])
                )
                if any(denominator)==0:
                    mass_density_term = np.array([0, 0, 0], dtype="float64")
                else:
                    mass_density_term = numerator / denominator
                mass_density += mass_density_term
                
        self.particle.mass_density = -self.OTHER_PARAMS["diffusion_coefficient"] * mass_density
        
    def update_conservative_diffusion_energy(self):
        velocity = 0
        for id, nbr_particle in enumerate(self.neighbours_list):
            if self.update_sound_speed(id) == nbr_particle.sound_speed:
                numerator = (
                     (nbr_particle.velocity - self.particle.velocity) *
                     (self.update_sound_speed(id) - nbr_particle.sound_speed) *
                     (self.PARAMETERS["cell_size"]*(self.particle.initial_pos -
                     nbr_particle.initial_pos) * self.cubic_spline_kernel_gradient(
                      self.particle.initial_pos - nbr_particle.initial_pos
                    ))
                )
                denominator = (
                   nbr_particle.mass_density*(m.pow(np.linalg.norm(
                    self.particle.initial_pos - nbr_particle.initial_pos), 2) +
                    self.OTHER_PARAMS["deformation"])
                )
                velocity += numerator / denominator
                
        self.particle.velocity = -self.OTHER_PARAMS["diffusion_coefficient"] * \
                       velocity
    
    # ------------------------------------------------------------------ FORCE ATTRIBUTES ------------------------------------------------------------------------

    def outward_facing_unit_norm(self, position):
        return self.normalize(position)
    
    def update_interface_normals(self):
        self.particle.interface_normal = np.array([0, 0, 0], dtype="float64")
        for id, nbr_particle in enumerate(self.neighbours_list):
            if (self.particle.sound_speed != nbr_particle.sound_speed).all():
                try:
                    mass_d = nbr_particle.mass / nbr_particle.mass_density
                except ZeroDivisionError:
                    mass_d = np.array([0, 0, 0], dtype="float64")
                self.particle.interface_normal = (
                    mass_d * self.cubic_spline_kernel_gradient(self.particle.initial_pos -
                                                               nbr_particle.initial_pos)
                )
        self.particle.interface_normal *= -1
        
    def interface_normal(self):
        self.update_interface_normals()
        return self.normalize(self.particle.interface_normal)

    def update_artificial_viscosity_slip(self):
        self.update_artificial_viscosity_no_slip()
        for nbr_particle in self.neighbours_list:
            self.particle.artificial_viscosity += (
                self.particle.artificial_viscosity
                * (self.interface_normal() * 
                self.normalize(self.particle.velocity -
                            nbr_particle.velocity)) *
                (self.normalize(nbr_particle.interface_normal) * 
                self.normalize(self.particle.velocity -
                            nbr_particle.velocity))
            )

    def update_sound_speed(self):
        for nbr_particle in self.neighbours_list:
            self.particle.sound_speed += (
                np.sqrt(self.OTHER_PARAMS["lambda_const"] * 
                        (self.particle.pressure - 
                        nbr_particle.pressure) / 
                        (self.particle.mass_density - 
                        nbr_particle.pressure))
            )
        return self.particle.sound_speed

    def update_artificial_viscosity_no_slip(self):
        
        for nbr_particle in self.neighbours_list:
            viscosity_term = (
                self.particle.velocity - nbr_particle.velocity *
                self.particle.initial_pos - nbr_particle.initial_pos *
                self.PARAMETERS["cell_size"] / 
                m.pow(np.linalg.norm(
                self.particle.initial_pos - nbr_particle.initial_pos
                ), 2)
            )
            viscosity = np.maximum(viscosity_term, np.array([0, 0, 0], dtype="float64"))
            const_term = -self.force / (self.particle.mass_density - nbr_particle.mass_density)
            self.particle.artificial_viscosity = (
                const_term * 
                (self.OTHER_PARAMS["alpha_const"])*
                self.particle.sound_speed*
                viscosity + 

                self.OTHER_PARAMS["beta_const"]*
                np.power(viscosity, 2)
            )

    def update_force(self):
        for nbr_particle in self.neighbours_list:
            self.particle.force += (
                (self.particle.acceleration + nbr_particle.acceleration) *
                (self.particle.half_velocity - nbr_particle.half_velocity) / 
                (2*self.particle.acceleration * (
                    self.particle.half_velocity - self.interface_velocity
                ) - 2*nbr_particle.acceleration*(
                    self.interface_velocity - nbr_particle.half_velocity
                ))
            )

    def update_half_velocity(self, particle):
        particle.half_velocity = particle.velocity + particle.acceleration*self.delta_time/2

    def update_mass_density(self, particle):

        self.update_half_velocity(particle)
        mass_density = np.array([0, 0, 0], dtype="float64")
        for nbr_particle in particle.neighbour_list:
            try:
                mass_d = nbr_particle.mass / nbr_particle.mass_density
            except ZeroDivisionError:
                mass_d = np.array([0, 0, 0], dtype="float64")
            mass_density += (
                mass_d *
                (particle.half_velocity - self.interface_velocity) *
                self.cubic_spline_kernel_gradient(particle.initial_pos - nbr_particle.initial_pos)
            )
        particle.mass_density *= -2*mass_density
    
    def update_interface_velocity(self):
        
        self.update_parallel_velocity()
        self.update_perpendicular_velocity()
        
        if self.slip_condition is True:
            self.update_slip_condition()
        elif self.slip_condition is not True:
            self.update_no_slip_condition()
        else:
            self.interface_velocity = self.particle.velocity/ np.sum([nbr_particle.velocity for nbr_particle in self.neighbours_list])

    def update_energy_conservation(self):
        
        velocity = np.array([[0, 0, 0],
                             [0, 0, 0],
                             [0, 0, 0]], dtype="float64")
        
        for id, nbr_particle in enumerate(self.neighbours_list):
            velocity += (
                self.particle.force * nbr_particle.mass *
                self.particle.acceleration *
                (self.particle.half_velocity - self.interface_velocity)
            )
        self.particle.velocity = velocity*2

    def update_net_forces(self):

        for nbr in self.neighbours_list:
            self.find_neighbour_list(nbr)
            self.update_mass_density(nbr)

        self.update_conservative_diffusion_mass_density()

        self.find_neighbour_list(self.particle)
        self.update_mass_density(self.particle)

        self.update_cauchy_stress()
        self.update_sound_speed()
        self.update_pressure()

        self.update_interface_velocity()

        self.update_force()

        if self.slip_condition is True:
            self.update_artificial_viscosity_slip()
        elif self.slip_condition is not True:
            self.update_artificial_viscosity_slip()
    
    # ------------------------------------------------------------------- UPDATE CALLS ---------------------------------------------------------------------------

    def debugging_forces(self, secs):

        print(f"Mass density is {self.particle.mass_density}")
        print(f"Acceleration is {self.particle.acceleration}")
        print(f"Interface velocity is {self.particle.interface_velocity}")
        print(f"Parallel velocity is {self.particle.parallel_velocity}")
        print(f"Perpendicular velocity is {self.particle.perp_velocity}")
        print(f"Cauchy stress is {self.particle.cauchy_stress_tensor}")
        print(f"Deviatoric stress is {self.particle.deviatoric_stress_tensor}")
        print("\n\n")
        time.sleep(secs)

    def update_velocity(self):

        self.acceleration = np.array([0, 0, 0], dtype="float64")
            
        for nbr_particle in self.neighbours_list:
            self.acceleration += nbr_particle.mass * (
                self.particle.acceleration - nbr_particle.acceleration
            )
        
        self.particle.velocity += self.delta_time*self.acceleration

    def update(self):

        self.update_net_forces()
        self.update_acceleration_force()
        self.update_velocity()

        self.particle.initial_pos += self.particle.velocity*self.delta_time
        self.XSPH_vel_correction()

        self.choose_collision_types("Cuboid", "Normal")
        self.debugging_forces(0.1)