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

julia> init(:H2, 1000.0, 1.59e6, H2 = 0.29, N2 = 0.56, O2 = 0.15)
# Creates an instance of the Gas struct based on the given mechanism and initial conditions

julia> gas.initial.temperature # T
1000.0

julia> gas.mechanism.molecular_weight # W
10-element Vector{Float64}:
  1.00784  # H
  2.01568  # H2
 15.9994   # O
 31.9988   # O2
 17.00724  # OH
 18.01508  # H2O
 28.0134   # N2
 33.00664  # HO2
 34.01448  # H2O2
 39.948    # AR
```
### Chemical Kinetics
```julia
julia> step!(gas, gas.initial.mass_fractions, gas.initial.temperature)
# Computes state and intermediate variables based on the given mass fractions and temperature

julia> gas.current.temperature # equals initial temperature
1000.0

julia> gas.intermediate.production_rate # ω̇
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
  
julia> gas.intermediate.temperature_change_rate # Ṫ
1-element Vector{Float64}:
 -0.12101709573654974
```
## References
- [Chemkin User Manual](https://www3.nd.edu/~powers/ame.60636/chemkin2000.pdf)
