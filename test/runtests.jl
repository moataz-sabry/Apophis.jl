using Apophis, PyCall, Test, Unitful

include("cantera.jl")
include("complex.jl")

function test_calculations(mech::Symbol)
    @testset verbose = true "Mechanism: $mech" begin
        gas_Apophis = Gas(mech)
        gas_Cantera = run_cantera(mech)
        rnd = rand(length(gas_Apophis.mechanism.species))
        for _ in 1:1
            T = rand(300.0:3000.0)
            P = rand(0.5Apophis.Pa:2Apophis.Pa)
            Y = rnd / sum(rnd)
            
            TPY!(gas_Apophis, T, P, Y) |> update
            update(gas_Apophis, :dT)
            gas_Cantera.TPY = T, P * ustrip(u"Pa", 1u"dyn/cm^2"), Y
            #show(stdout, "text/plain", gas_Cantera.third_body_concentrations)
            #show(stdout, "text/plain", [gas_Cantera.third_body_concentrations map(r -> r isa Union{Apophis.ThreeBodyReaction, Apophis.FallOffReaction} ? Apophis.total_molar_concentration(gas_Apophis.state.C, r.enhancement_factors) : NaN, reactions(gas_Apophis))])
            #show(stdout, "text/plain", gas_Cantera.third_body_concentrations)
            test_species_calculations(gas_Apophis, gas_Cantera, T)
            #test_species_derivatives(gas_Apophis, gas_Cantera)
        end
    end
end

function test_calculations_derivatives(mech::Symbol)
    @testset verbose = true "Mechanism: $mech" begin
        real_gas = Gas(mech)
        complex_gas = Gas(mech; as=ComplexF64)
        rnd = rand(length(real_gas.mechanism.species))
        for _ in 1:1
            Tc = (rand(300.0:3000.0)) + 0im
            Pc = (rand(0.5Apophis.Pa:2Apophis.Pa)) + 0im
            Yc = (rnd / sum(rnd)) .+ 0im
            
            T, P, Y = real(Tc), real(Pc), real(Yc)
            TPY!(real_gas, T, P, Y)
            TPY!(complex_gas, Tc, Pc, Yc)
            check_derivatives(real_gas, complex_gas)
        end
    end
end

@testset verbose = true "Apophis vs. Cantera" begin
    for mech in (:H2, :GRI12, :GRI2, :GRI3, :ITV)
            @testset verbose = true "Reader" test_reader(mech)
            #@testset verbose = true "Calculations" test_calculations(mech)
            test_calculations_derivatives(mech)
    end
end