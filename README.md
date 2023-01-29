# Apophis
A package for the sensitivity analysis of gas-phase chemical and kinetics

## Installation
```julia
pkg> add https://github.com/moataz-sabry/Apophis.jl
```
## Quick Start
### Reading a Mechanism
```julia
julia> using Apophis

julia> Gas(:GRI3; T = 1000, P = 1013250, Y = "CH4: 0.05, O2: 0.20, N2: 0.75")
# Creates an instance of the Mechanism struct based on the given mechanism data files
```
### Setting the State
```julia
julia> init(:H2, 1000.0, 1.59e6; H2 = 0.29, N2 = 0.56, O2 = 0.15)
# Creates an instance of the Gas struct based on the given mechanism and initial conditions

julia> gas.state.T # T [K]
1000.0

julia> molecular_weights(gas) # W [g/mole]
10-element Vector{Float64}:
  1.00784  # H
  2.01568  # H₂
 15.9994   # O
 31.9988   # O₂
 17.00724  # OH
 18.01508  # H₂O
 28.0134   # N₂
 33.00664  # HO₂
 34.01448  # H₂O₂
 39.948    # Ar
```
### Chemical Kinetics
```julia
julia> update!(gas)
# Computes state and intermediate variables based on the given mass fractions and temperature

julia> production_rates(gas) # ω̇ [mole/(cm³⋅s)]
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
- [Chemkin Transport Manual](https://www3.nd.edu/~powers/ame.60636/transport.pdf)
