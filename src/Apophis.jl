module Apophis

export gas
export init
export step!
export equilibrate
export readmechanism
export sensitivity
export subspaces
export monte

using LinearAlgebra
using SparseArrays: sparse, SparseMatrixCSC
using DifferentialEquations
using Sundials
using Distributions

const Rc = 1.987261815324 # * u"cal * (K * mol)"
const R = 8.31446261815324e7 # * u"erg / (K * mol)"
const Pₐ = 1013250.0 # * u"dyn / cm^2"
const d = 0.14

include("types.jl")
include("utils.jl")
include("complex.jl")
include("forward.jl")
include("solvers.jl")
include("subspaces.jl")

end
