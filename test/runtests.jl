using Apophis, PyCall, Test, Unitful
ct = pyimport("cantera")

include("tests/utils.jl")
include("tests/calc.jl")
include("tests/reader.jl")
#include("complex.jl")

@testset verbose = true "Apophis vs. Cantera" begin
    for mech in (:H2, :GRI12)#)
        @testset "Mechanism: $mech" begin
            @testset verbose = true "Reader" test_reader(mech)
            @testset verbose = true "Values" test_values(mech)
        end
    end
end