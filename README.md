# Apophis
Welcome to Apophis! This package is designed to provide sensitivity analysis of combustion processes using efficient and accurate derivatives. It is built to be fast, reliable, and easy to use; but still early in development!

## Installation
```julia
pkg> add Apophis
```
## Quick Start
### Reading a Mechanism
```julia
julia> using Apophis
julia> using Unitful: K, Pa

julia> gas = Gas(:GRI3; T = 1000K, P = 101325Pa, Y = "CH4: 0.05, O2: 0.20, N2: 0.75")
# Creates a Gas based on the given mechanism data files

              T: 1000 K  //  P: 101325 Pa  //  ρ: 0.3372 kg/m³
 --------- ------ ------ --------- ------------ -------------- ------------
  Species      Y      X         C           cₚ              h            s 
               –      –   kmol/m³   J/(kmol⋅K)         J/kmol   J/(kmol⋅K) 
 --------- ------ ------ --------- ------------ -------------- ------------
       N2   0.75   0.74    0.0090        32762    2.14699e+07       228089
       O2   0.20   0.17    0.0021        34883    2.27068e+07       243586
      CH4   0.05   0.09    0.0011      73616.7   -3.59484e+07       248279
     ⋮       ⋮      ⋮        ⋮          ⋮             ⋮             ⋮
 --------- ------ ------ --------- ------------ -------------- ------------
                                                             50 rows omitted
julia> temperature(gas) # T [K]
1000.0

julia> density(gas) # ρ [kg/m³]
0.337205
```
### Changing the State
```julia
julia> set!(gas; T = 1250.0, ρ = 0.35, Y = mass_fractions(gas)); # T && (P || ρ) && (Y || X || C)

julia> pressure(gas) # P [Pa]
131461.83
```
### Chemical Kinetics
```julia
julia> update(gas);
# Updates internal variables based on the current gas state

julia> production_rates(gas; view=false) # ω̇ [kmol/(m³⋅s)]
53-element Vector{Float64}:
  0.0
  8.034465079230524e-7
  1.6467013518389382e-12
 -2.8148788209137565e-5
  0.0
  0.0
  2.8148786697190487e-5
  0.0
  0.0
  0.0
  ⋮
  0.0
  0.0
  0.0
 -2.6060314085092194e-12
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