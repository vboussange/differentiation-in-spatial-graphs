# Neutral and adaptive differentiation in spatial graphs

This repository contains the code used for the paper 

> *Eco-evolutionary model on spatial graphs reveals how habitat structure affects phenotypic differentiation*, Boussange, V. & Pellissier, L. (2022)

- `code/` contains all scripts related to the simulation runs
    - `graph_utils` contains scripts to calculate graph properties
    - `simulations` contains the scripts to generate the simulations related to the paper settings.
- `figure/` contains all scripts to generate the manuscript figures and crunch the raw simulation results.

All scripts are written in the Julia programming language. A short description of the purpose of each script is placed in each script preamble.
The scripts can be executed out of the box by activating the environment stored in the `code/simulations/Project.toml` file.

To activate the environment in an interactive session, type in the Julia REPL

```julia
julia>] activate .
julia>] instantiate
```
To simply run a script, type in the terminal
```
> julia --project=. name_of_the_script.jl
```

