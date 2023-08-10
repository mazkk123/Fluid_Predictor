# Fluid_Predictor
Using a mixture of CNN (Convolutional Neural Networks)
and machine learning methodologies, this project aims to 
predict particle positions, velocities and accelerations
in multiple SPH approximations of a fluid.

## Weighting Kernels
In an SPH approximation of fluids, different weighting 
kernels are used as approximators for computing various
forces within a neighbourhood vicinity of a particle in
the system.

### KGC - Kernel Gradient Correction
In large scale simulations and for SPH stability, KGC
is required to reduce weighting kernel steps from being 
too large by rescaling them based on a kernel correction 
factor.

### Poly 6 
For mass density and surface tension calculations the poly
6 kernel is used and generates the most stable results.

### Spiky
The spiky kernel type is mostly used for pressure force and
pressure calculations based on extensive experimentation.

### Viscosity
As the name suggests, the viscosity kernel is used to 
approximate the viscosity force/ laminar viscosity that
is acting on a particle.

### Cubic Spline
Some SPH implementations choose to use a full Cubic Spline
approximation, whereby all force computations utilize a
cubic spline weighting kernel in their calculations.

## Forces
A variety of physical forces exist per particle which
contrbute to their velocity and position staps between
timesteps

## Presure Forces
These are unknown forces that exist between particles
contributing to their intermolecular state.

### Pressure

## Advective Forces
These forces contribute to the general flow of a particle
by affecting its velocity and position.

### Viscosity
This determines the clumping between particle neighbourhoods
in the system. It is the intermolecular friction that 
exists between particles which determines whether a fluid 
appears water like or globular.

#### Artificial Viscosity
An artificial viscosity term can be added in certain situations
to ensure better stability in the simulation

#### Laminar Viscosity
this is the main type of viscosity that exists for laminar 
fluid flows. These type of flows are quite stable and non
turbulent.

### Thermal Diffusion
This force determines how a fluid should react in different
temperatures and should the fluid expand in these cases.

### Buoyancy
A force that controls the upward trajectory of the fluid. 
Liquids have no buoyancy and therefore stick to the ground.
Buoyant fluids such as gasses rise and diffuse through a 
containing space due to the contribution of this force.

### Surface Tension
This is the friction that exists between the fluid and a 
solid boundary wall. Controlling this force can create splashy
or more subdued droplets when a fluid is poured or dropped
onto another fluid interface.

### Body Forces
Additional body force exists to provide the fluid with some
mass and initial density when it reacts with different interfaces.

## Fluid Models
Below are the fluid models used within this framework
to the user's disposal.

### SPH - Smooth Particle Hydrodynamics
The first implementation type is using the average
SPH model with different nearest neighbour search algorithms
to the user's disposal. The algorithm works by using 
weighting kernels to approximate the average physical forces
present between different neighbourhood clusters of 
particles approximated using a nearest neighbour search
algorithms and hash tables.

### Multi SPH - Mixing particles with attributes
A multi SPH approach introduces the presence of different 
attributes in the traditional SPH approach to fluids. Instead
of treating fluids of uniform mass, mass density and physical
forces, particles within a neighbourhood are split further
into different phases within the same neighbourhood. Further 
computations are applied when particles with different attributes
coexist in the same neighbourhood cluster

### WCSPH - Weakly Compressible SPH
This method extends the computations undertaken in the traditional
SPH by guaranteeing better density fluctuation and advanced time
step computations. By doing so, it ensures that fluids can achieve,
not a fully incompressible state, but one that is weakly 
compressible and accurately depicts time stepping and density 
fluctuations for more particles. 

### PCSPH - Predictive Corrective SPH
Pressure instability is commonplace with the traditional SPH 
approximation of fluids. To achieve further fluid incompressibility,
velocity and position values are first used to predict pressure 
and pressure force. While density error of the fluid is within a  
threshold, pressure correction values are added to the pressure until 
this condition is invalidated. The benefit of this approximation is to 
stabilize pressure between discrete time steps within the simulation.

### LPSPH - Local Poisson SPH
Another implementation that improves pressure instabilities between
time steps to further achieve incompressibility. In contrast to 
the PCSPH approximation, this implementation categorizes particles 
within a local domain into pressure neighbourhoods of a particle. 
Delineating particles beyond this threshold as far, and those within
as near. In addition to the traditional SPH scheme, if the density
error is within a threshold, computations on far and near particles within
their respective pressure neighbourhoods are undergone. The sum of
individual pressures is added to the main pressure calculation as
the compensation step to ensure pressure stability.

### IISPH - Implicit Incompressible SPH
When particles reach solid boundaries or frontiers,
the normal SPH model assumes uniform velocities and particle
separation. With this model, drift velocity and velocity advection
is taken into account whenever particles interact with solid 
boundaries and reach above a certain velocity threshold.

### FSISPH - Fluid Solid Interface SPH
Enforcing incompressibility and density variance near interfacial
regions in mixing fluids is a limitation of the traditional SPH
model. Here, intermediate velocity, slip conditions and conservative
diffusion is introduced near interfacial mixing regions to precipitate
more turbulent behaviour when different fluids mix together.

### DFSPH - Divergence Free SPH
In this model, there is better splashy behaviour near solid
boundaries. 2 additional correction steps are undertaken 
before every velocity update which ensure that divergence is
zero and predicted densities match the rest density of the fluid
in question.

### VCSPH - Vorticity Confinement SPH
An extended version of the DFSPH model, which introduces a vorticity
step using the curl of the fluid and a stream function to update 
velocity based on the vorticity of the fluid. Using the last velocity 
update step in the DFSPH model, we compute vorticity and stream 
function forces and reupdate velocity values based on these. This
structure also implements look-up tables for weighted kernels for faster
computation times.

### PBF - Particle Based Fluid
Another model which introduces vorticity confinement and a constraint 
function to the fluid. As a result, particles are more close together 
when solid objects are dropped into the fluid tank. An XSPH viscosity 
term is also introduced to ensure more coherent motion.

# Neighbourhood search algorithms

## Spatial Hashing
This nearest neighbour search algorithm is the most
memory optimized using O(1) complexity. It spatially 
partitions a grid of particles based on cell size and
particle separation amounts. Particles hashed to the same
cell are considered neighbourhood pairs and this is 
computed per frame and on initialization of the simulation

## Squared Distance 
The simplest to implement but the most performance un-
optimized method with O(n^2) complexity. Particles within
a certain squared distance threshold of one another are 
candidates to be neighbourhood pairs. The computational 
complexity of this method is because finding a neighbouring
pair requires one to search all particles in the system 
during simulation.

## Compact Hashing
Very similar to the Spatial Hashing algorithm, but particle
data is stored into linked lists instead of generic list
datatypes. The benefit of this implementation is towards
performance and stability of the SPH simulation. Using linked
lists removes any unused data spaces between hashed particles
to the same list using a traditional Spatial Hashing approach.

# Time Stepping

## Normal time steps
In traditional time stepping, velocity and positions calculations
are computed per particle basis using a constant incremental time
step amount delta t.

### Euler Cromer
This basic time step implementation suggests that particle 
velocities should be updated based on the sum of all forces divided
by the mass of a particle in question. Then using this acceleration
to update a particle's velocity and position linearly

### Leap Frog
Through this approach, time steps are calculated in 1/2 steps 
forward and backward. Velocity has to be precomputed, which is then
used further to approximate position and velocity in a frame and half
time based on a frame and a half before.

### Verlet
This scheme uses the acceleration of a particle in a full time 
step forward added with the current acceleration at this time step
to approximate the velocity and position at the next time step
for a particle.

## Adaptive Time Stepping
In adaptive time stepping, to consider the amount that a particle
should be able to move between a time step, a further CFL - Couratt
Fredewy Levy condition states that the time step delta t between 2
frames should not be greater than the maximum cell size spacing 2
particle cell neighbourhoods in the system. The CFL condition has to
be satisfied for larger scale simulations operating near boundary walls.

### Regional/Independent Time
With this time stepping technique, delta t is instead approximated
using the average mass density divided by the cell length separating 2
particle neighbourhoods. This is multiplied by a CFL factor to constrain
the allowed position of a particle to be within a reasonable scale, and 
not beyond the cell spacing amount.

### Runga Kutta

# Multithreading and parallelization
This project uses multi-threading and multiprocessing by 
using the python Threading and queue libraries to distribute 
force computations across multiple threads and processes, then 
reassemble the value computed the quickest. In the future, a 
GPU implementation with CUDA is expected to attain even faster
results and simulate greater numbers of particles.

## Importing data
It is important that data from other packages are 
imported back into this application as either txt or 
JSON formats as well. Failure to do so can cause 
numerical instabilities in the calculations implemented
into the application.








