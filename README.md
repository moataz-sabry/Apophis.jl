# Apophis
A package for sensitivity analysis of gas chemical and kinetics

## Installation
```julia
] add Apophis
```
## Quick Start
### Reading a Mechanism
```julia
using Apophis

gas = Gas(:GRI3; T = 1000, P = 1013250, Y = "CH4: 0.05, O2: 0.20, N2: 0.75")
# Creates an instance of the Mechanism struct based on the given mechanism data files

temperature(gas) ## T [K]
1000.0

density(gas) ## ρ [g/cm³]
0.000337205
```
### Changing the State
```julia
TρY!(gas, 1250.0, 0.0035, mass_fractions(gas))

pressure(gas) ## P [dyn/cm²]
1.314618e7
```
### Chemical Kinetics
```julia
update(gas)
# Updates internal variables based on the current gas state

production_rates(gas) # ω̇ [mol/(cm³⋅s)]
53-element Vector{Float64}:
  0.0
  3.231563056183965e-8
  1.6467013518389592e-13
 -2.8148788209138157e-6
  0.0
  0.0
  2.814878669719108e-6
  ⋮
  0.0
  
```
## To-Do
- Change default units to SI
- Add routines for other reaction auxillary parameters
- Correctly implement pressure-dependant reactions
- Implement diffusion routines

## References
- [Chemkin User Manual](https://www3.nd.edu/~powers/ame.60636/chemkin2000.pdf)
- [Chemkin Transport Manual](https://www3.nd.edu/~powers/ame.60636/transport.pdf)
