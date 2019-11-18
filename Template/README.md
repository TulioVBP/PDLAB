# How to set up an experiment
Open `experiment_template.m` to start to set up the experiment. Follow each section below.

## Parameters
Here we define the parameters of the simulation
### Material's properties 
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