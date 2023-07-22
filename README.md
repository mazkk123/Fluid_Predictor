# Fluid_Predictor
Using a mixture of CNN (Convolutional Neural Networks)
and machine learning methodologies, this project aims to 
predict particle positions, velocities and accelerations
in multiple SPH approximations of a fluid.

## Traditional SPH - Smooth Particle Hydrodynamics
The first implementation type is using the average
SPH model with different nearest neighbour search algorithms
to the user's disposal. The algorithm works by using 
weighting kernels to approximate the average physical forces
present between different neighbourhood clusters of 
particles approximated using a nearest neighbour search
algorithms and hash tables.

## Multi SPH - Mixing particles with attributes
A multi SPH approach introduces the presence of different 
attributes in the traditional SPH approach to fluids. Instead
of treating fluids of uniform mass, mass density and physical
forces, particles within a neighbourhood are split further
into different phases within the same neighbourhood. Further 
computations are applied when particles with different attributes
coexist in the same neighbourhood cluster

## WCSPH - Weakly Compressible SPH
This method extends the computations undertaken in the traditional
SPH by guaranteeing better density fluctuation and advanced time
step computations. By doing so, it ensures that fluids can achieve
- not a fully incompressible state - but one that is weakly 
compressible and accurately depicts time stepping and density 
fluctuations for more particles. 

## PCSPH - Predictive Corrective SPH

## IISPH - Implicit Incompressible SPH
When particles reach solid boundaries or frontiers,
the normal SPH model assumes uniform velocities and particle
separation. With this model, drift velocity and velocity advection
is taken into account whenever particles interact with solid 
boundaries and reach above a certain velocity threshold

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

# Multithreading and parallelization
This project uses multi-threading and multiprocessing by 
using the python Threading and queue libraries to distribute 
force computations across multiple threads and processes, then 
reassemble the value computed the quickest. In the future, a 
GPU implementation with CUDA is expected to attain even faster
results and simulate greater numbers of particles.

# GUI creation
The graphical interface for this project was designed
using pyside6 and qtcreator. It is inspired by the interfaces 
implemented into modern DCC applications such as Maya, Blender 
and Houdini. Focused purely on fluid simulations, much of the 
tree view and nodal graph representations are overlooked for this
implementation

## Menu Bars
In the menu bars, users are given options to save,
copy, redo, import and export. The import and export 
options are purely to enable users to easily transfer
their simulation data into other DCC packages such as 
Houdini. This application exports data into JSON and txt
formats currently, but further data formats are seen to
be added in the future.

## Importing data
It is important that data from other packages are 
imported back into this application as either txt or 
JSON formats as well. Failure to do so can cause 
numerical instabilities in the calculations implemented
into the application.

## Frame buffer system
To facilitate faster frame rates, Data is stored into
a temporary frame buffer whenever a simulation is run. 
The visualized display is essentially a cached state of 
JSON data that is being fetched from the frame buffer.

## Dockable widgets
All widgets in this project are tab dockable just for 
better UX and UI design. Custom functions exist so that
widgets can be minimized and rescaled whenever a user 
wishes to fully display their simulation options or fully 
display the simulation screen fully.






