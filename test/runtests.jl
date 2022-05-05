using Apophis, LinearAlgebra
using Test

include("tests.jl")

@testset verbose = true "Apophis.jl" begin
    @testset verbose = true "mechanism parameters" begin end

    @testset verbose = true "intermediate variables & derivatives 'AD vs. CSD'" begin # for loooop
        tests = 1:1
        mechanism = [:H2, :SynGas, :GRI12, :GRI2, :GRI3]#, :ITV]

        for m in mechanism
            @testset verbose = true "$m" begin

                for _ in tests
                    T = rand(750.0:2500.0)
                    P = rand(0.8e6:1.59e6)
                    ns = length(readmechanism(m, Float64).species)
                    fractions = normalize(rand(ns), 1)

                    init(m, T + 0im, P + 0im; rand=fractions)
                    complexgas = gas

                    init(m, T, P; rand=fractions)
                    realgas = gas

                    @testset verbose = true "values" begin
                        checkvalues(realgas, complexgas)
                    end
                    @testset verbose = true "derivatives" begin
                        checkderivatives(realgas, complexgas)
                    end
                end
            end
        end
    end
end