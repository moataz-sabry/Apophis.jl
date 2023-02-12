using Apophis, PyCall, Test, Unitful
ct = pyimport("cantera")

include("tests/utils.jl")
include("tests/reader.jl")
include("tests/values.jl")
include("tests/derivatives.jl")

@testset verbose = true "Apophis vs. Cantera" begin
    for mech in (:H2, :GRI12, :GRI2, :GRI3, :ITV)#)
        @testset verbose = true "Mechanism: $mech" begin
            @testset "Reader" test_reader(mech)
            @testset "Values" test_values(mech)
            #@testset "Derivatives" test_derivatives(mech)
        end
    end
end