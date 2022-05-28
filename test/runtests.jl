using Apophis
using Test

include("tests.jl")

@testset verbose = true "Apophis.jl" begin
    @testset verbose = true "mechanism parameters" begin end

    @testset verbose = true "intermediate variables & derivatives 'AD vs. CSD'" begin # for loooop
        tests = 1:2
        mechanism = [:H2, :SynGas, :GRI12, :GRI2, :GRI3, :ITV]

        for m in mechanism
            @testset verbose = false "$m" begin

                for t in tests
                    T = rand(750.0:1250.0)
                    P = rand(0.8e6:1.59e6)
                    fuel = rand(0.05:0.3)
                    oxygen = rand(0.15:0.5)
                    nitrogen = 1.0 - fuel - oxygen

                    init(m, T + 0im, P + 0im, H2=fuel, N2=nitrogen, O2=oxygen)
                    complexgas = gas

                    init(m, T, P, H2=fuel, N2=nitrogen, O2=oxygen)
                    realgas = gas

                    @testset verbose = false "values" begin
                        checkvalues(realgas, complexgas)
                    end
                    @testset "derivatives" begin
                        checkderivatives(realgas, complexgas)
                    end
                end
            end
        end
    end
end