using Test

@testset verbose = true "Apophis.jl" begin

    @testset "H2" begin
        H2 = readmechanism(:H2)
        @test H2.molecular_weight[1] == 1.00784
    end

    @testset "GRI3" begin
        GRI3 = readmechanism(:GRI3)
        @test GRI3.molecular_weight[1] == 1.00784
    end
end