module Apophis

using Base.Iterators: filter, flatten, partition, reverse, take
using Base: Fix1, Fix2, OneTo, rest
using PrettyTables
using SparseArrays
using SplitApplyCombine: combinedimsview, mapview
using Unitful
using YAML: load_file
using Zygote

abstract type AbstractSpecies{N<:Number} end
abstract type AbstractReaction{N<:Number} end

const R = 8.31446261815324e3 # J K⁻¹ kmol⁻¹
const Rc = 1.987261815324 # cal K⁻¹ mol⁻¹
const kB = 1.380649e-23 # J K⁻¹
const Pa = 101325 # Pa
const Tᵣ = 300 # K
const d = 0.14

const Maybe{T} = Union{T, Nothing}
const Energies{N<:Number} = NamedTuple{(:val, :dT), Tuple{Vector{N}, Vector{N}}}
const Thermodynamics{N<:Number} = NamedTuple{(:cₚ, :h, :s), NTuple{3, Energies{N}}}
const Rates{N<:Number} = NamedTuple{(:val, :dT, :dC), NTuple{3, Vector{N}}}
const SpeciesRates{N<:Number} = NamedTuple{(:ω̇,), Tuple{Rates{N}}}
const ReactionRates{N<:Number} = NamedTuple{(:kf, :kr, :q), NTuple{3, Rates{N}}}

include("species.jl")
include("gas.jl")
include("reactions.jl")
include("reader.jl")
include("calc.jl")
include("sensitivity.jl")
include("utils.jl")

end