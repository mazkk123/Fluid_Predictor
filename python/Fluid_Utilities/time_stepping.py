import math as m
import numpy as np
import random as rd

from Particles.particles import Particle

class ForwardEuler:
    def __init__(self,
                 particle: Particle,
                 delta_time: float):
        
        self.particle = particle
        self.delta_time = delta_time

        self.exec_time_scheme()
    
    def exec_time_scheme(self):

        self.particle.velocity += self.delta_time*self.particle.acceleration
        self.particle.initial_pos += self.particle.velocity*self.delta_time

class EulerCromer:
    
    def __init__(self,
                 particle: Particle,
                 delta_time: float):
        
        self.particle = particle

        self.velocity = np.array([0, 0 ,0])
        self.position = np.array([0, 0 ,0])

        self.exec_time_scheme()
        self.exec_pseudo_time_scheme(delta_time)

    def exec_time_scheme(self, delta_time):

        self.particle.velocity += delta_time*self.particle.acceleration
        self.particle.initial_pos += self.particle.velocity*delta_time

    def exec_pseudo_time_scheme(self, delta_time):

        self.velocity += delta_time*self.particle.acceleration
        self.position += self.particle.velocity*delta_time
        
    def get_time_scheme_values(self, delta_time):

        self.exec_pseudo_time_scheme(delta_time)
        return (self.particle.initial_pos, self.particle.velocity)

# ----------------------------------------------------------------------- ADVANCED TIME STEPPING ---------------------------------------------------------------
    
class LeapFrog(EulerCromer):
    def __init__(self,
                 particle: Particle,
                 delta_time: float):
        
        super().__init__()

        self.particle = particle
        self.delta_time = delta_time

        self.exec_time_scheme()
        self.values = self.get_time_scheme_values(self.delta_time/2)
    
    def exec_Leap_time_scheme(self):

        self.particle.velocity += self.delta_time*self.particle.acceleration
        self.particle.initial_pos += self.particle.velocity*self.delta_time

    def get_time_scheme_values(self, delta_time):
        return super().get_time_scheme_values(delta_time)
    
class Verlet(EulerCromer):

    def __init__(self,
                 particle: Particle,
                 delta_time: float):
        
        super().__init__()

        self.particle = particle
        self.delta_time = delta_time

        self.exec_time_scheme()
    
    def exec_time_scheme(self):

        self.particle.velocity += self.delta_time * (self.particle.acceleration + self.particle.next_acceleration) * \
                                    self.delta_time
        self.particle.initial_pos += self.particle.velocity*self.delta_time + \
                                    0.5 * self.particle.acceleration * m.pow(self.delta_time, 2)
        
    def get_time_scheme_values(self, delta_time):
        return super().get_time_scheme_values(delta_time)

# --------------------------------------------------------------------- ADAPTIVE TIME STEPPING --------------------------------------------------------------------
class IndependentTime(EulerCromer):
    
    CFL_NUM = 0.3

    def __init__(self,
                 particle:Particle = None,
                 cell_size:float = None):

        super().__init__()

        if cell_size is not None:
            self.delta_time = self.CFL_NUM*cell_size/np.linalg.norm(particle.velocity)
        
        self.exec_time_scheme()

    def exec_time_scheme(self):

        return super().exec_time_scheme(self.delta_time)

class RegionalShortTime(IndependentTime):
    
    def __init__(self,
                 particle: Particle=None,
                 cell_size:float = None):
        
        super().__init__(particle, cell_size)

class RungeKutta:
    pass