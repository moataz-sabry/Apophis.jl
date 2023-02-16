# Apophis
Welcome to Apophis! This package is designed to provide sensitivity analysis of combustion processes using efficient and accurate derivatives. It is built to be fast, reliable, and easy to use, with a range of features and tools that can help to get the most out of combustion simulations.

## Installation
```julia
] add Apophis
```
## Quick Start
### Reading a Mechanism
```julia
using Apophis

gas = Gas(:GRI3; T = 1000, P = 101325, Y = "CH4: 0.05, O2: 0.20, N2: 0.75")
# Creates a Gas based on the given mechanism data files

temperature(gas) ## T [K]
1000.0

density(gas) ## ρ [kg/m³]
0.337205
```
### Changing the State
```julia
TρY!(gas, 1250.0, 0.35, mass_fractions(gas))

pressure(gas) ## P [Pa]
131461.83
```
### Chemical Kinetics
```julia
update(gas)
# Updates internal variables based on the current gas state

production_rates(gas) # ω̇ [kmol/(m³⋅s)]
53-element Vector{Float64}:
  0.0
  3.7072899145013144e-11
  1.3302904784766728e-16
 -8.845602918023088e-8
  0.0
  0.0
  8.845602904821038e-8
  0.0
  0.0
  0.0
  ⋮
  0.0
  0.0
  0.0
 -1.434894871257002e-16
  0.0
  0.0
  0.0
  0.0
  0.0
  
```
## To-Do
- Add routines for other reaction auxillary parameters
- Implement diffusion routines

## References
- [Chemkin User Manual](https://www3.nd.edu/~powers/ame.60636/chemkin2000.pdf)
- [Chemkin Transport Manual](https://www3.nd.edu/~powers/ame.60636/transport.pdf)
