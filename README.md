# Apophis
Chemistry Analysis Package

## Installation
```julia
using Pkg
Pkg.add(path=".../Apophis")
```
## Quick Start
### Setting the State
```julia
julia> using Apophis

julia> init(:H2, 1000.0, 1.59e6, H2 = 29, N2 = 56, O2 = 15)
# Creates an instance of the Gas struct based on the given mechanism and initial conditions

julia> gas.initial.temperature
1000.0

julia> gas.initial.mass_fractions
10-element Vector{Float64}:
 0.0
 0.29
 0.0
 0.15
 0.0
 0.0
 0.56
 0.0
 0.0
 0.0
```
### Chemical Kinetics
```julia
julia> step!(gas, gas.initial.mass_fractions, gas.initial.temperature)
# Computes state and intermediate variables based on the given mass fractions and temperature

julia> gas.intermediate.production_rate
10-element Vector{Float64}:
  2.2300800003459409e-10
 -2.2300796859512652e-10
  2.6106322069659418e-21
 -2.2300793715696425e-10
  0.0
  0.0
  0.0
  2.2300793715565893e-10
  0.0
  0.0
```
