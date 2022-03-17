# Neutral and adaptive differentiation in spatial graphs

This repository contains code used for the paper 

> *Neutral and adaptive differentiation in spatial graphs*, Boussange, V. & Pellissier, L. (2022)

- `code/` contains all scripts related to the simulation runs
    - `graph_utils` contains scripts to calculate graph properties
    - `simulations` contains the script to generate the simulations related to the paper settings
- `figure/` contains all scripts to create graphical output and crunch the raw simulation results

All scripts are written in the Julia programming language. You can download the software at https://julialang.org .
The scripts can be executed out of the box by activating the environment stored in the `code/simulations/Project.toml` file.
To activate the environment, type in the Julia REPL

```julia
julia>] activate code/simulations
```
or run 
```
> julia --project=code/simulations name_of_the_script.jl
```

