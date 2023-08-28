import math as m
import numpy as np
import random as rd
import re
import time

from Fluid_Calculations.compute_sph import SPH
from Particles.particles import Particle

class IISPH(SPH):

    OTHER_PARAMS = {
        "max_iterations":2,
        "relaxation_factor":0.5
    }

    def __init__(self, 
                 particle: Particle=None,
                 search_method: str=None,
                 all_particles:list = None,
                 time_schemes:dict = None,
                 collision_types:dict = None,
                 params:dict = None,
                 hash_table:dict=None,
                 time_stepping:str = "Euler Cromer",
                 tank_attrs:dict = None,
                 additional_params:dict = None,
                 hash_value:int=None,
                 delta_time:float = None):
        
        super().__init__(particle=particle,
                        search_method=search_method,
                        all_particles=all_particles,
                        time_stepping=time_stepping,
                        time_schemes=time_schemes,
                        params=params,
                        collision_types=collision_types,
                        tank_attrs=tank_attrs,
                        hash_table=hash_table,
                        hash_value=hash_value,
                        delta_time=delta_time)
        
        self.additional_params = additional_params

        self.OTHER_PARAMS["max_iterations"] = self.additional_params["max_iterations"]
        self.OTHER_PARAMS["relaxation_factor"] = self.additional_params["relaxation_factor"]

        self.predict_advection()
        self.iteration = 0

    def calculate_density_error(self):
        return self.particle.mass_density - self.PARAMETERS["mass_density"]

    # ----------------------------------------------------------------------- PREDICT ADVECTION -----------------------------------------------------------------------------

    def predict_nbr_advection(self, particle, depth:int=3):

        if depth==0:
            return
        
        self.find_neighbour_list(particle)
        self.predict_velocity_advection(particle)
        self.predict_displacement(particle)
        self.update_acceleration_advection(particle)

        for nbr in particle.neighbour_list:
            return self.predict_nbr_advection(nbr, depth-1)

    def predict_advection(self):

        self.mass_density_adv = 0
        self.predict_nbr_advection(self.particle, 2)

        self.particle.pressure  = self.particle.prev_pressure*0.5
        self.update_mass_density_advection()
    
    def predict_displacement(self, particle):
        displacement = 0
        for id, nbr_particle in enumerate(particle.neighbour_list):
            try:
                density_div = -nbr_particle.mass/ m.pow(particle.mass_density, 2)
            except ZeroDivisionError:
                density_div = 0
            kernel_gradient = self.cubic_spline_kernel_gradient(particle.initial_pos - nbr_particle.initial_pos)
            displacement += density_div*kernel_gradient

        particle.displacement = displacement*m.pow(self.delta_time, 2)
    
    def predict_velocity_advection(self, particle):

        self.find_neighbour_list(particle)
        self.update_predicted_mass_density(particle)
        self.update_predicted_gravity(particle)
        self.update_predicted_buoyancy(particle)
        self.update_predicted_viscosity(particle)
        self.update_predicted_surface_tension(particle)

        self.all_forces = particle.gravity + \
                            particle.buoyancy + \
                            particle.viscosity + \
                            particle.surface_tension
        
        particle.advected_velocity = particle.velocity + self.delta_time*(self.all_forces / particle.mass)     
    
    def update_mass_density_advection(self):

        advection_amt = 0
        for id, nbr_particle in enumerate(self.neighbours_list):
            mass_advection = nbr_particle.mass * (self.particle.advected_velocity - nbr_particle.advected_velocity)
            kernel_gradient = self.cubic_spline_kernel_gradient(self.particle.initial_pos - nbr_particle.initial_pos)
            advection_amt += mass_advection*kernel_gradient
        self.mass_density_adv = advection_amt*self.delta_time + self.particle.mass_density

    def update_acceleration_advection(self, particle):
        acc_adv = 0
        for id, nbr_particle in enumerate(particle.neighbour_list):
            acc_adv += (
                nbr_particle.mass * (particle.displacement - nbr_particle.displacement ) *
                self.cubic_spline_kernel_gradient(particle.initial_pos - nbr_particle.initial_pos)
            )
        particle.acceleration_adv = acc_adv
    
    # ------------------------------------------------------------------------ PRESSURE SOLVE ---------------------------------------------------------------------------------

    def debugging_forces(self, secs):

        print("Mass Density:", self.particle.mass_density)
        print("Pressure:", self.particle.pressure)
        print("Pressure Force:", self.particle.pressure_force)
        print("Buoyancy:", self.particle.buoyancy)
        print("Gravity:", self.particle.gravity)
        print("viscosity:", self.particle.viscosity)
        print("Delta time is: ", self.delta_time)
        print("\n\n")
        time.sleep(secs)

    def pressure_solve(self):

        self.iteration = 0

        while self.calculate_density_error() > 0.01*self.PARAMETERS["mass_density"] \
                and self.iteration<self.OTHER_PARAMS["max_iterations"]:

            self.update_iter_pressure()

            self.iteration +=1 

        self.particle.prev_pressure = self.particle.pressure
        self.update_pressure_force()

    def update_pressure_force(self):
        return super().update_pressure_force()
    
    def update_displacement_iter(self, particle):
        displacement = 0
        for id, nbr_particle in enumerate(self.neighbours_list):
            try:
                density_div = -1* nbr_particle.mass/ m.pow(nbr_particle.mass_density, 2)
            except ZeroDivisionError:
                density_div = 0
            kernel_gradient = self.cubic_spline_kernel_gradient(particle.initial_pos - nbr_particle.initial_pos)
            displacement += density_div*nbr_particle.pressure*kernel_gradient

        particle.displacement_iter = displacement*m.pow(self.delta_time, 2)
    
    def neighbour_dis_diff(self):
        
        for nbr in self.neighbours_list:
            self.find_neighbour_list(nbr)
            for nbrs_nbr in nbr.neighbour_list:
                self.update_displacement_iter(nbrs_nbr)
                
        nbr_displacement_diff = 0
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

        displacement_diff = 0
        for nbr in self.neighbours_list:
            displacement_diff += (
                self.particle.displacement_iter - nbr.displacement_iter
            )
        return displacement_diff
        
    def update_iter_pressure(self):
        relaxation_const = (1 - self.OTHER_PARAMS["relaxation_factor"])*self.particle.pressure

        acc_adv = 0
        try:
            if self.particle.acceleration_adv == 0:
                acc_adv = 0
            else:
                1/self.particle.acceleration_adv
        except ZeroDivisionError:
            acc_adv = 0

        acceleration_adv_const = self.OTHER_PARAMS["relaxation_factor"] * acc_adv
        disp_final = 0
        for nbr_particle in self.particle.neighbour_list:
            displacement_diff = self.displacement_diff()
            nbr_displacement_diff = self.neighbour_dis_diff()
            disp_final += (displacement_diff - nbr_displacement_diff) * \
                        self.cubic_spline_kernel_gradient(self.particle.initial_pos - nbr_particle.initial_pos) * \
                        nbr_particle.mass
            
        self.particle.pressure = (self.PARAMETERS["mass_density"] - self.mass_density_adv - disp_final) \
                                    *acceleration_adv_const + relaxation_const
    
    # -------------------------------------------------------------------------- UPDATE STEP ----------------------------------------------------------------------------------

    def update(self):
        
        self.pressure_solve()
        
        self.particle.velocity = self.particle.advected_velocity + self.delta_time*self.particle.pressure_force/self.particle.mass
        self.particle.initial_pos += self.particle.velocity*self.delta_time
        
        self.XSPH_vel_correction()
        self.choose_collision_types("Cuboid", "Normal")

        #self.adapt_to_CFL()