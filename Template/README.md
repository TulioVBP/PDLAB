# How to set up an experiment
Open `experiment_template.m` to start to set up the experiment. Follow each section below.

## Parameters
Here we define the parameters of the simulation
### Material properties 
The material properties are adimensioneless (so be consistent with units)

* Young modulus 
* Density 
* Poisson's ratio  
* Energy release rate
* Horizon
### Numerical parameters
* Mesh ratio
* Influence function
* Partial area algorithm (PA-AC, PA-HHB, FA)
* Surface effect correction algorithm (None, Volume method)
### Mesh parameters
The `mesh.generateMesh()` returns a rectangular mesh with `a` length and `b` height. If you want to generate your own geometry, please edit the file `generateMesh.m` to return your specific mesh. The function should return

* Position matrix `x = [xi yi]`
* Element area vector `A = [Ai]`
### Initial cracks
Initial cracks are introduced in the matrix `damage.crackIn`, where each row is composed by a line segment (first collumn is the initial point,second column the final point).
###  Model choice
Here we define

* If damage is on (turn it off if quasi-static solver is chosen)
* Model's name
* Solver (Quasi-Static or Dynamic-Explicit)

## Simulation
The simulation conditions are established here.

+ `b_parallelComp` is true if one wants to run the code in parallel;
+ `stress_app` is a variable that specifies how the boundary conditions are applied. If equal to `'-'`, `mesh.boundaryCondition()` function will look for `pc` variable with the personalized constraints (in the future, this will be the only option).
+ `prescribedBC()` is a user defined function that returns the a structure `pc` with the variables
    1. `pc.disp = [nodeId, boolean for u_x, boolean for u_y, u_x, u_y]`, which the matrix with displacement constraints (e.g., prescribing displacement 3 in x direction for node 56: `pc.disp = [56, true, 0, 3, 0]`) 
    2. `pc.vel = [nodeId, boolean for v_x, boolean for v_y, v_x, v_y]`, which is the matrix with velocity constraints
    3. `pc.bodyForce = [nodeId, b_x, b_y] `, which is the matrix the prescribed body force
+ `mesh.boundaryCondition()` defines the body force vector, the constrained dof and the damage free nodes (`noFailZone`)
+ `neighborhood.generateFamily_v2()` generates the family matrix (each row i collects the neighbors of node i), the partial area matrix, and the surface effect correction matrix. It saves a file called `family1.mat` in the FAMILY FILES folder
+ `history = models.historyDependency()` initializes the history dependence matrix, which depends on the model used

## Solver
There are two types of solver

* Dynamic Explicit Time solver
    1. Define the time increment `dt`, final time `t_tot` and the frequency data is stored in `data_dump` (e.g., if `data_dump = 4`, the solver will store data every 4 time steps). The greater `data_dump` is, the faster the simulation will go since the solver won't evaluate quantities as energy for more time steps.
    2. Choose a filename to save the output
* Quasi-Static Implicit Time Solver
    1. Define the number of load steps `n_tot`

## Post-processing
The embedded post-processing functions will plot the outputs, such as energy, damage, displacement field and strains.

