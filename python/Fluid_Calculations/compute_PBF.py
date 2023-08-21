import math as m
import numpy as np
import time

from Fluid_Calculations.compute_sph import SPH
from Particles.particles import Particle

class PBF(SPH):
    
    OTHER_PARAMS = {
        "relaxation_factor":0.7,
        "k_const":0.1,
        "n_const":4,
        "solver_iterations":2,
        "c_const":0.01
    }

    def __init__(self,
                particle: Particle=None,
                search_method:str = None,
                hash_table:dict = None,
                hash_value:int = None,
                time_stepping:str = "Euler Cromer",
                tank_attrs:dict = None,
                all_particles:list = None,
                delta_time:int = 0.02):
        
        super().__init__(particle=particle,
                         search_method=search_method,
                         hash_table=hash_table,
                         hash_value=hash_value,
                         time_stepping=time_stepping,
                         tank_attrs=tank_attrs,
                         delta_time=delta_time)
        
        self.all_particles = all_particles
        self.constraint_grad = np.array([0, 0, 0], dtype="float64")
        self.iter = 0
        self.depth_lvl = 3
        self.main_depth_lvl = 3
        
        self.predict_positions(self.particle, self.main_depth_lvl)
        
    def predict_positions(self, particle, depth:int=4):
        
        if depth==0:
            return
            
        self.find_neighbour_list(particle)
        self.update_predicted_mass_density(particle)
        self.update_predicted_gravity(particle)
        self.update_predicted_buoyancy(particle)
        self.update_predicted_viscosity(particle)
        self.update_predicted_surface_tension(particle)
        
        self.all_forces += (
            particle.gravity +
            particle.buoyancy +
            particle.viscosity +
            particle.surface_tension
        )
        
        particle.predicted_velocity = particle.velocity + self.delta_time* \
                                        self.all_forces
        particle.predicted_initial_pos = particle.initial_pos + self.delta_time* \
                                        particle.predicted_velocity
        for nbr in particle.neighbour_list:
            return self.predict_positions(nbr, depth-1)
    
    def update_constraint_function(self, particle):
        
        particle.constraint_function = (
            particle.mass_density / self.PARAMETERS["mass_density"] - 1
        )

    def update_constraint_grad(self, particle):

        self.constraint_grad = np.array([0, 0, 0], dtype="float64")            
        for nbr_particle in particle.neighbour_list:
            self.constraint_grad += (
                self.cubic_spline_kernel_gradient(particle.predicted_initial_pos - 
                                            nbr_particle.predicted_initial_pos)
                )
        self.constraint_grad *= 1/self.PARAMETERS["mass_density"]
        return self.constraint_grad

    def update_constraint(self, particle):
        
        constraint_sum = 0
        for p in particle.neighbour_list:
            constraint_sum += (
                m.pow(np.linalg.norm(self.update_constraint_grad(p)), 2)
            )
        particle.constraint = particle.constraint_function / (constraint_sum + self.OTHER_PARAMS["relaxation_factor"])

    def update_s_corr(self, particle, id):
        delta_q_const = np.array([0.1*self.PARAMETERS["cell_size"],
                                  0.1*self.PARAMETERS["cell_size"],
                                  0.1*self.PARAMETERS["cell_size"]], dtype="float64")
        s_corr = 0
        s_corr = m.pow((
            -self.OTHER_PARAMS["k_const"]*
            (
                self.cubic_spline_kernel_pos(particle.predicted_initial_pos - 
                                             particle.neighbour_list[id].predicted_initial_pos) / 
                self.cubic_spline_kernel_pos(delta_q_const)
            )
        ), self.OTHER_PARAMS["n_const"])
        return s_corr
    
    def update_del_position(self, particle):
        
        constraint_term = np.array([0, 0, 0], dtype="float64")
        for id, nbr_particle in enumerate(particle.neighbour_list):
            constraint_term += (
                self.particle.constraint + nbr_particle.constraint + self.update_s_corr(particle, id)*
                self.cubic_spline_kernel_gradient(particle.predicted_initial_pos - nbr_particle.predicted_initial_pos)
            )
        particle.del_position = 1/self.PARAMETERS["mass_density"]*constraint_term

    def update_position(self, particle):
        particle.predicted_initial_pos += particle.del_position
            
    def update_vorticity(self, particle):
        for nbr_particle in particle.neighbour_list:
            particle.vorticity += (
                np.cross((nbr_particle.velocity - particle.velocity),
                          self.update_constraint_grad(nbr_particle))
            )

    def n_const(self, particle):
        return self.new_cubic_spline_kernel_gradient(np.linalg.norm(particle.vorticity))
    
    def update_vorticity_force(self, particle):
        try:
            curl_term = np.cross(self.n_const(particle)/ np.linalg.norm(self.n_const(particle)), particle.vorticity)
        except np.AxisError:
            curl_term = np.array([0, 0, 0], dtype="float64")
        particle.vorticity_force = (
            self.OTHER_PARAMS["relaxation_factor"]*curl_term
        )
        
    def update_vorticity_forces(self, particle):

        self.update_vorticity(particle)
        self.update_vorticity_force(particle)
    
    def update_velocity(self, particle):
        particle.predicted_velocity = 1 / self.delta_time * (particle.predicted_initial_pos - particle.initial_pos)
    
    def XSPH_viscosity_correction(self):
        correction_term = np.array([0, 0, 0], dtype="float64")
        for nbr in self.particle.neighbour_list:
            correction_term += (
                (self.particle.velocity - nbr.velocity) *
                self.cubic_spline_kernel_pos(self.particle.initial_pos - nbr.initial_pos)
            )
        correction_term *= self.OTHER_PARAMS["c_const"]
        self.particle.velocity += correction_term
        
    def update_forces(self, particle):
        
        self.all_forces = np.array([0, 0, 0], dtype="float64")
        self.all_forces += (
            particle.buoyancy +
            particle.viscosity +
            particle.surface_tension +
            particle.gravity +
            particle.vorticity_force
        )
        particle.acceleration = self.all_forces / particle.mass
        particle.velocity = particle.predicted_velocity + self.delta_time*particle.acceleration
        particle.initial_pos = particle.predicted_initial_pos + self.delta_time*particle.velocity
    
    def update_secondary_forces(self, particle, depth:int = 2):

        if depth==0:
            return
        
        self.update_velocity(particle)
        self.update_vorticity_forces(particle)

        for nbr in particle.neighbour_list:
            return self.update_secondary_forces(nbr, depth-1)
    
    def debugging_forces(self, secs):
        print("Vorticity is:", self.particle.vorticity)
        print("Vorticity force is: ", self.particle.vorticity_force)
        print("\n\n")
        time.sleep(secs)

    def update_constraints(self, particle, depth:int=4):

        if depth==0:
            return
        
        self.find_neighbour_list(particle)
        self.update_constraint_function(particle)
        self.update_constraint(particle)
        self.update_del_position(particle)
        self.update_position(particle)

        for nbr in particle.neighbour_list:
            return self.update_constraints(nbr, depth-1)

    def update(self):

        while(self.iter < self.OTHER_PARAMS["solver_iterations"]):

            self.update_constraints(self.particle, self.depth_lvl)
            
            self.iter += 1

        self.update_secondary_forces(self.particle, self.depth_lvl)
        
        self.update_forces(self.particle)
        self.debugging_forces(0.01)

        self.XSPH_vel_correction()
        self.XSPH_viscosity_correction()

        self.adapt_to_CFL()
        