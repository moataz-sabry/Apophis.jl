using Apophis
using Test

include("tests.jl")

@testset verbose = true "Apophis.jl" begin
    @testset verbose = true "mechanism parameters" begin end

    @testset verbose = true "intermediate variables & derivatives 'AD vs. CSD'" begin # for loooop
        tests = 1:5

        @testset verbose = false "H2" begin

            for t in tests
                T = rand(750.0:1250.0)
                P = rand(0.8e6:1.59e6)
                fuel = rand(0.05:0.3)
                oxygen = rand(0.15:0.5)
                nitrogen = 1.0 - fuel - oxygen

                init(:H2, T + 0im, P + 0im, H2=fuel, N2=nitrogen, O2=oxygen)
                complexgas = gas

                init(:H2, T, P, H2=fuel, N2=nitrogen, O2=oxygen)
                realgas = gas

                @testset verbose = false "values" begin
                    checkvalues(realgas, complexgas)
                end
                @testset "derivatives" begin
                    checkderivatives(realgas, complexgas)
                end
            end
        end

        @testset verbose = false "SynGas" begin

            for t in tests
                T = rand(750.0:1250.0)
                P = rand(0.8e6:1.59e6)
                fuel = rand(0.05:0.3)
                oxygen = rand(0.15:0.5)
                nitrogen = 1.0 - fuel - oxygen

                init(:SynGas, T + 0im, P + 0im, H2=fuel, N2=nitrogen, O2=oxygen)
                complexgas = gas

                init(:SynGas, T, P, H2=fuel, N2=nitrogen, O2=oxygen)
                realgas = gas

                @testset verbose = false "values" begin
                    checkvalues(realgas, complexgas)
                end
                @testset "derivatives" begin
                    checkderivatives(realgas, complexgas)
                end
            end
        end

        @testset verbose = false "GRI 1.2" begin

            for t in tests
                T = rand(750.0:1250.0)
                P = rand(0.8e6:1.59e6)
                fuel = rand(0.05:0.3)
                oxygen = rand(0.15:0.5)
                nitrogen = 1.0 - fuel - oxygen

                init(:GRI12, T + 0im, P + 0im, CH4=fuel, N2=nitrogen, O2=oxygen)
                complexgas = gas

                init(:GRI12, T, P, CH4=fuel, N2=nitrogen, O2=oxygen)
                realgas = gas

                @testset verbose = false "values" begin
                    checkvalues(realgas, complexgas)
                end
                @testset "derivatives" begin
                    checkderivatives(realgas, complexgas)
                end
            end
        end

        @testset verbose = false "GRI 2.0" begin

            for t in tests
                T = rand(750.0:1250.0)
                P = rand(0.8e6:1.59e6)
                fuel = rand(0.05:0.3)
                oxygen = rand(0.15:0.5)
                nitrogen = 1.0 - fuel - oxygen

                init(:GRI2, T + 0im, P + 0im, CH4=fuel, N2=nitrogen, O2=oxygen)
                complexgas = gas

                init(:GRI2, T, P, CH4=fuel, N2=nitrogen, O2=oxygen)
                realgas = gas

                @testset verbose = false "values" begin
                    checkvalues(realgas, complexgas)
                end
                @testset "derivatives" begin
                    checkderivatives(realgas, complexgas)
                end
            end
        end

        @testset verbose = false "GRI 3.0" begin

            for t in tests
                T = rand(750.0:1250.0)
                P = rand(0.8e6:1.59e6)
                fuel = rand(0.05:0.3)
                oxygen = rand(0.15:0.5)
                nitrogen = 1.0 - fuel - oxygen

                init(:GRI3, T + 0im, P + 0im, CH4=fuel, N2=nitrogen, O2=oxygen)
                complexgas = gas

                init(:GRI3, T, P, CH4=fuel, N2=nitrogen, O2=oxygen)
                realgas = gas

                @testset verbose = false "values" begin
                    checkvalues(realgas, complexgas)
                end
                @testset "derivatives" begin
                    checkderivatives(realgas, complexgas)
                end
            end
        end

        @testset verbose = false "ITV" begin

            for t in 1:2
                T = rand(750.0:1250.0)
                P = rand(0.8e6:1.59e6)
                fuel = rand(0.05:0.3)
                oxygen = rand(0.15:0.5)
                nitrogen = 1.0 - fuel - oxygen

                init(:ITV, T + 0im, P + 0im, CH4=fuel, N2=nitrogen, O2=oxygen)
                complexgas = gas

                init(:ITV, T, P, CH4=fuel, N2=nitrogen, O2=oxygen)
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