import math as m
import numpy as np
import random as rd
import re

from Fluid_Calculations.compute_sph import SPH
from Particles.particles import Particle

class IISPH(SPH):

    OTHER_PARAMS = {
        "lambda_stiffness":7,
        "k_stiffness":7,
        "max_iterations":2,
        "relaxation_factor":0.5
    }

    def __init__(self, 
                 particle: Particle=None,
                 search_method: str=None,
                 hash_table:dict=None,
                 time_stepping:str = "Euler Cromer",
                 tank_attrs:dict = None,
                 hash_value:int=None,
                 params:dict = None,
                 delta_time:float = None):
        
        super().__init__(particle=particle,
                        search_method=search_method,
                        time_stepping=time_stepping,
                        tank_attrs=tank_attrs,
                        hash_table=hash_table,
                        hash_value=hash_value,
                        params=params,
                        delta_time=delta_time)
        
        self.predict_advection()

    def calculate_density_error(self):
        return self.particle.mass_density - self.PARAMETERS["mass_density"]

    def cubic_spline_kernel_gradient(self, position):
        q = np.linalg.norm(position) / self.PARAMETERS["cell_size"]
        kernel_const = 1 / (np.pi*m.pow(self.PARAMETERS["cell_size"], 4))
        if q>=0 and q<=1:
            kernel_val = 9/4*m.pow(q, 2) - 3*q
            return kernel_val*kernel_const
        if q>=1 and q<=2:
            kernel_val = -3/4*m.pow((2 - q), 2)
            return kernel_val*kernel_const
        if q>=2:
            return 0
    
    # ---------------------------------------------------------------------- PREDICT FORCES --------------------------------------------------------------------------------

    def update_mass_density(self, particle):
        """
        """
        density = 0 
        for nbr_particle in particle.neighbour_list:
            kernel_value = self.kernel_linear(particle.initial_pos - nbr_particle.initial_pos, 0)
            density += kernel_value*nbr_particle.mass

        particle.mass_density = self.PARAMETERS["mass_density"] + density
    
    def update_viscosity(self, particle):
        """
        """
        viscosity = np.array([0, 0, 0], dtype="float64")
        for id, nbr_particle in enumerate(particle.neighbour_list):
            vel_dif = nbr_particle.velocity - particle.velocity
            kernel_laplacian = self.kernel_laplacian(particle.initial_pos - nbr_particle.initial_pos, 2)
            try:
                mass_pressure = particle.mass/nbr_particle.mass_density
            except ZeroDivisionError:
                mass_pressure = 0
            viscosity += vel_dif*mass_pressure*kernel_laplacian

        particle.viscosity = viscosity*self.PARAMETERS["viscosity"]

    def update_gravity(self, particle):
        """
        """
        particle.gravity = self.gravity_const
    
    def update_normal_field(self, particle):
        """
        """
        normal_field = np.array([0, 0, 0],  dtype="float64")
        for id, nbr_particle in enumerate(particle.neighbour_list):
            pos_difference = particle.initial_pos - nbr_particle.initial_pos
            normal_field += (
                particle.mass * 1/particle.mass_density * self.kernel_gradient(pos_difference, 0)
            )
        return normal_field

    def update_surface_curvature(self, particle):
        """        
        """
        surface_curvature = 0
        for id, nbr_particle in enumerate(particle.neighbour_list):
            surface_curvature += (
                particle.mass * 1/particle.mass_density * self.kernel_laplacian(
                particle.initial_pos - nbr_particle.initial_pos, 0)
            )
        return surface_curvature
    
    def update_surface_tension(self, particle):
        """
        """
        normal_field = self.update_normal_field(particle)
        surface_curvature = self.update_surface_curvature(particle)
        normal_field_magnitude = np.linalg.norm(normal_field)
        
        if normal_field_magnitude >= self.PARAMETERS["tension_threshold"]:
            particle.surface_tension = (
                self.PARAMETERS["tension_coefficient"] * surface_curvature * normal_field/normal_field_magnitude 
            )

    def update_buoyancy(self, particle):
        """
        """
        buoyancy = self.PARAMETERS["buoyancy"] * (particle.mass_density - self.PARAMETERS["mass_density"])
        buoyancy *= self.gravity_const
        particle.buoyancy = buoyancy

    # ----------------------------------------------------------------------- PREDICT ADVECTION -----------------------------------------------------------------------------

    def predict_advection(self):

        self.mass_density_adv = 0
        self.acceleration_adv = 0

        for nbr in self.neighbours_list:
            self.predict_displacement(nbr)
            self.predict_velocity_advection(nbr)
            self.update_acceleration_advection(nbr)

        self.predict_displacement(self.particle)
        self.predict_velocity_advection(self.particle)
        self.update_acceleration_advection(self.particle)

        self.particle.pressure  = self.particle.prev_pressure*0.5
        self.update_mass_density_advection()
    
    def predict_displacement(self, particle):
        displacement = np.array([0, 0, 0])
        for id, nbr_particle in enumerate(particle.neighbour_list):
            density_div = -nbr_particle.mass/ m.pow(particle.mass_density, 2)
            kernel_gradient = self.cubic_spline_kernel_gradient(particle.initial_pos - nbr_particle.initial_pos)
            displacement += density_div*kernel_gradient

        particle.displacement = displacement*m.pow(self.delta_time, 2)

    def find_neighbour_list(self, particle):
        particle.neighbour_list = []
        for nbr in self.hash_table[particle.hash_value]:
            if particle is not nbr:
                particle.neighbour_list.append(nbr)
    
    def predict_velocity_advection(self, particle):

        self.find_neighbour_list(particle)
        self.update_mass_density(particle)
        self.update_gravity(particle)
        self.update_buoyancy(particle)
        self.update_viscosity(particle)
        self.update_surface_tension(particle)

        self.all_forces = particle.gravity + \
                            particle.buoyancy + \
                            particle.viscosity + \
                            particle.surface_tension
        
        particle.advected_velocity = particle.velocity + self.delta_time*(self.all_forces / particle.mass)     
    
    def update_mass_density_advection(self):

        advection_amt = np.array([0, 0, 0])
        for id, nbr_particle in enumerate(self.neighbours_list):
            mass_advection = nbr_particle.mass * (self.particle.advected_velocity - self.particle.velocity)
            kernel_gradient = self.cubic_spline_kernel_gradient(self.particle.initial_pos - nbr_particle.initial_pos)
            advection_amt += mass_advection*kernel_gradient
        self.mass_density_adv = advection_amt*self.delta_time + self.particle.mass_density

    def update_acceleration_advection(self, particle):
        acc_adv = np.array([0, 0 , 0])
        for id, nbr_particle in enumerate(particle.neighbour_list):
            acc_adv += (
                nbr_particle.mass * (particle.displacement - nbr_particle.displacement ) *
                self.cubic_spline_kernel_gradient(particle.initial_pos - nbr_particle.initial_pos)
            )
        particle.acceleration_adv = acc_adv
    
    # ------------------------------------------------------------------------ PRESSURE SOLVE ---------------------------------------------------------------------------------

    def pressure_solve(self):
        iteration = 0
        while self.calculate_density_error() < 0.01*self.PARAMETERS["mass_density"] \
                and iteration<self.OTHER_PARAMS["max_iterations"]:
            
            self.update_displacement_iter()
            self.update_iter_pressure()
            self.particle.pressure = self.particle.iter_pressure

            iteration +=1 
            
        self.particle.prev_pressure = self.particle.pressure
        self.update_pressure_force()

    def update_displacement_iter(self, particle):
        displacement = np.array([0, 0, 0])
        for id, nbr_particle in enumerate(self.neighbours_list):
            density_div = -1* nbr_particle.mass/ m.pow(particle.mass_density, 2)
            kernel_gradient = self.cubic_spline_kernel_gradient(particle.initial_pos - nbr_particle.initial_pos)
            displacement += density_div*particle.pressure*kernel_gradient

        particle.displacement_iter = displacement*m.pow(self.delta_time, 2)
    
    def neighbour_dis_diff(self):
        
        for nbr in self.neighbours_list:
            self.find_neighbour_list(nbr)
            for nbrs_nbr in nbr.neighbour_list:
                self.update_displacement_iter(nbrs_nbr)
                
        nbr_displacement_diff = np.array([0, 0, 0], dtype="float64")
        for nbr in self.neighbours_list:
            for nbrs_nbr in nbr.neighbour_list:
                nbr_displacement_diff += (
                    nbrs_nbr.displacement_iter - nbr.displacement_iter
                )
        return nbr_displacement_diff
        
    def displacement_diff(self):
        for nbr in self.neighbours_list:
            self.update_displacement_iter(nbr)
        
        self.update_displacement_iter(self.particle)
        displacement_diff = np.array([0, 0, 0], dtype="float64")
        for nbr in self.neighbours_list:
            displacement_diff += (
                self.displacement_iter - nbr.displacement_iter
            )
        return displacement_diff
        
    def update_iter_pressure(self, particle):
        relaxation_const = (1 - self.OTHER_PARAMS["relaxation_factor"])*particle.pressure
        acceleration_adv_const = self.OTHER_PARAMS["relaxation_factor"] * 1/self.acceleration_adv
        disp_final = np.array([0, 0, 0])
        for id, nbr_particle in enumerate(particle.neighbour_list):
            displacement_difference = self.displacement_diff()
            nbr_displacement_difference = self.neighbour_dis_diff()
            disp_final += (displacement_diff - nbr_displacement_diff) * \
                        self.cubic_spline_kernel_gradient(self.particle.initial_pos - nbr_particle.initial_pos) * \
                        nbr_particle.mass
            
        self.particle.prev_pressure = (self.mass_density_adv - self.PARAMETERS["mass_density"] - disp_final) \
                                    *acceleration_adv_const + relaxation_const
    
    # -------------------------------------------------------------------------- UPDATE STEP ----------------------------------------------------------------------------------

    def update(self):
        
        self.pressure_solve()
        
        self.particle.velocity = self.particle.advected_velocity + self.delta_time*self.particle.pressure_force/self.particle.mass
        self.particle.initial_pos += self.particle.velocity*self.delta_time
        self.XSPH_vel_correction()

        self.choose_collision_types("Cuboid", secondary_type="Normal", particle=self.particle)