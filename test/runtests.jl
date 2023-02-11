using Apophis, PyCall, Test, Unitful
ct = pyimport("cantera")

include("tests/utils.jl")
include("tests/reader.jl")
#include("complex.jl")

@testset verbose = true "Apophis vs. Cantera" begin
    for mech in (:H2, :GRI12)#)
        test_reader(mech)
    end
end