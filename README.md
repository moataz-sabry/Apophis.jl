# Apophis
A package for the sensitivity analysis of gas-phase chemical and kinetics

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
 "H"
 "H2"
 "O"
 "O2"
 "OH"
 "H2O"
 "N2"
 "HO2"
 "H2O2"
 "AR"
 
julia> gas.initial.mass_fractions
10-element Vector{Float64}:
 0.0   # H
 0.29  # H2
 0.0   # O
 0.15  # O2
 0.0   # OH
 0.0   # H2O
 0.56  # N2
 0.0   # HO2
 0.0   # H2O2
 0.0   # AR
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
## References
- [Chemkin User Manual](https://www3.nd.edu/~powers/ame.60636/chemkin2000.pdf)
