module Apophis

using Base.Iterators: filter, flatten, reverse, take, partition
using Base: Fix1, Fix2, OneTo, rest
using Combinatorics
using LinearAlgebra
using SparseArrays
using SplitApplyCombine: combinedimsview, mapview
using Unitful
using YAML: load_file
using Zygote
# using ThreadsX

abstract type AbstractReaction{N} end

const R = 8.31446261815324e7 # erg/(K*mol)
const Rc = 1.987261815324 # cal/(K*mol)
const kB = 1.380649e-16 # erg/K
const Pa = 1013250.0 # dyn/cm^2
const Tᵣ = 300.0 # K
const d = 0.14

const Maybe{T} = Union{T, Nothing} where {T}
const Diffusions = (:binary_diffusions, :mean_diffusions)
const Energies{N} = NamedTuple{(:val, :dT), Tuple{Vector{N}, Vector{N}}} where {N<:Number}
const Thermodynamics{N} = NamedTuple{(:cₚ, :h, :s), NTuple{3, Energies{N}}} where {N<:Number}
const Rates{N} = NamedTuple{(:val, :dT, :dC), NTuple{3, Vector{N}}} where {N<:Number}
const SpeciesRates{N} = NamedTuple{(:ω̇,), Tuple{Rates{N}}} where {N<:Number}
const ReactionRates{N} = NamedTuple{(:kf, :kr, :q), NTuple{3, Rates{N}}} where {N<:Number}

include("species.jl")
include("state.jl")
include("reactions.jl")
include("reader.jl")
include("calc.jl")
include("utils.jl")

end