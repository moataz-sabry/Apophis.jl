using Apophis, PyCall, Test, Unitful
ct = pyimport("cantera")

include("tests/utils.jl")
include("tests/reader.jl")
include("tests/values.jl")
include("tests/derivatives.jl")

@testset verbose = true "Apophis vs. Cantera" begin
    for mech in (:H2, :GRI12)#)
        @testset verbose = true "Mechanism: $mech" begin
            @testset verbose = false "Reader" test_reader(mech)
            #@testset verbose = false "Values" test_values(mech)
            #@testset verbose = true "Derivatives" test_derivatives(mech)
        end
    end
end