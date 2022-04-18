module Apophis

export gas
export init
export step!
export jac
export equilibrate

#import BenchmarkTools: @btime ## only for tests?
using LinearAlgebra
using SparseArrays: sparse, SparseMatrixCSC
using Symbolics: build_function, jacobian, derivative, Num
using DifferentialEquations: ODEProblem, solve
using Sundials: CVODE_BDF

const Rc = 1.987261815324 #* u"cal * (K * mol)"
const R = 8.31446261815324e7 # * u"erg / (K * mol)"
const Pₐ = 1013250.0 #* u"dyn / cm^2"
const d = 0.14

include("structs.jl")
include("utils.jl")
include("forward.jl")
#include("/Users/sabry/Developer/Apophis/jac/H2.jl")

#include("../jac/H2.jl")

#precompile(jac, (Float64, Vector{Float64})) ## works?

end
