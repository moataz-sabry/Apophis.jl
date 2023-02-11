########################################################################  Species  ########################################################################

function test_enthalpies(species_Apophis, species_Cantera, T)
    species_Apophis_enthalpy = enthalpy(species_Apophis; in=u"J/(kmol*K)")
    species_Cantera_enthalpy = species_Cantera.partial_molar_enthalpies
    try
        @test species_Apophis_enthalpy ≈ species_Cantera_enthalpy
    catch e
        error("Enthalpies comparison: Apophis ($species_Apophis_enthalpy) & Cantera ($species_Cantera_enthalpy) for species $(species_Apophis.formula) at temperature $T")
        rethrow(e)
    end
end

test_entropies(species_Apophis, species_Cantera, T) = @test entropy(species_Apophis; in=u"J/kmol") ≈ species_Cantera.thermo.s(T)
test_capacities_pressure(species_Apophis, species_Cantera, T) = @test heat_capacity_pressure(species_Apophis; in=u"J/kmol") ≈ species_Cantera.thermo.cp(T)

function test_production_rates(gas_Apophis, gas_Cantera)
    @testset "production rates" begin
        for (i, v) in enumerate(gas_Cantera.net_production_rates)
            @test isapprox(Apophis.production_rates(gas_Apophis, in=u"kmol/m^3/s")[i], v, rtol = 0.005)
        end
    end
end

########################################################################  Reactions  ########################################################################

function test_equilibrium_constants(gas_Apophis, gas_Cantera)
    @testset verbose = false "equilibrium constants" begin 
        for (i, v) in enumerate(gas_Cantera.equilibrium_constants)
            reaction_Apophis = reaction(gas_Apophis, i)
            ∑v = reaction_Apophis.reaction_order
            Kc = Apophis.equilibrium_constants(reaction_Apophis, gas_Cantera.T) * ustrip(uparse("(kmol/m^3)^$∑v/s"), 1uparse("(mol/cm^3)^$∑v/s"))
            @test Kc ≈ v
        end
    end
end

function test_total_molar_concentrations(gas_Apophis, gas_Cantera)
    @testset verbose = false "total molar concentrations" begin 
        for (i, v) in enumerate(gas_Cantera.third_body_concentrations)
            if !isnan(v)
                reaction_Apophis = reaction(gas_Apophis, i)
                M = Apophis.total_molar_concentration(gas_Apophis.state.C, reaction_Apophis.enhancement_factors) * ustrip(uparse("kmol/m^3"), 1uparse("mol/cm^3"))
                @test M ≈ v rtol = 0.0001
            end
        end
    end
end

function test_forward_rate_constants(gas_Apophis, gas_Cantera)   
    @testset verbose = false "forward rate constants" begin 
        for (i, v) in enumerate(gas_Cantera.forward_rate_constants)
            reaction_Apophis = gas_Apophis.mechanism.reactions[i]
            ∑νᵣ = sum(v -> last(v) |> abs, reaction_Apophis.reactants) + (reaction_Apophis isa Apophis.ElementaryReaction ? 0.0 : 1.0)
            unitify = 10.0^3(1-∑νᵣ)
            unitifyM = 10.0^3
            M = reaction_Apophis isa Apophis.ThreeBodyReaction ? unitifyM * Apophis.total_molar_concentration(gas_Apophis.state.C, reaction_Apophis.enhancement_factors) : one(Float64)
            @test reaction_Apophis.rates.kf.val[] * unitify * M ≈ v rtol = 0.005
        end
    end
end

function test_forward_rate_of_progress(gas_Apophis, gas_Cantera)   
    @testset verbose = false "forward rates of progress" begin 
        for (i, v) in enumerate(gas_Cantera.forward_rates_of_progress)
            reaction = Apophis.reaction(gas_Apophis, i)
            Π = Apophis.step(reaction.reactants, gas_Apophis.state.C)
            M = reaction isa Apophis.ThreeBodyReaction ? Apophis.total_molar_concentration(gas_Apophis.state.C, reaction.enhancement_factors) : one(Float64)
            #println(gas_Cantera.reaction(i-1), " – ", v, " – ", M)
            @test reaction.rates.kf.val[] * M * Π * ustrip(uparse("kmol/m^3/s"), 1uparse("mol/cm^3/s")) ≈ v rtol = 0.005
        end
    end
end

function test_reverse_rate_constants(gas_Apophis, gas_Cantera)
    @testset verbose = false "reverse rate constants" begin 
        for (i, v) in enumerate(gas_Cantera.reverse_rate_constants)
            ∑νᵣ = sum(v -> last(v) |> abs, gas_Apophis.mechanism.reactions[i].reactants) + (gas_Apophis.mechanism.reactions[i] isa Apophis.ElementaryReaction ? 0.0 : 1.0)
            unitifykf = 10.0^3(1-∑νᵣ)
            n = gas_Apophis.mechanism.reactions[i].reaction_order
            unitify = 10.0^-3n
            unitifyM = 10.0^3
            M = gas_Apophis.mechanism.reactions[i] isa Apophis.ThreeBodyReaction ? unitifyM * Apophis.total_molar_concentration(gas_Apophis.state.C, gas_Apophis.mechanism.reactions[i].enhancement_factors) : one(Float64)
            @test gas_Apophis.mechanism.reactions[i].rates.kr.val[] * unitifykf * unitify * M ≈ v rtol = 0.005
        end
    end
end

function test_reverse_rate_of_progress(gas_Apophis, gas_Cantera)
    @testset verbose = false "reverse rates of progress" begin 
        for (i, v) in enumerate(gas_Cantera.reverse_rates_of_progress)
            reaction = Apophis.reaction(gas_Apophis, i)
            Π = Apophis.step(reaction.products, gas_Apophis.state.C)
            M = reaction isa Apophis.ThreeBodyReaction ? Apophis.total_molar_concentration(gas_Apophis.state.C, reaction.enhancement_factors) : one(Float64)
            @test reaction.rates.kr.val[] * M * Π * ustrip(uparse("kmol/m^3/s"), 1uparse("mol/cm^3/s")) ≈ v rtol = 0.005
        end
    end
end

function test_progress_rates(gas_Apophis, gas_Cantera)   
    @testset verbose = false "progress_rate" begin
        for (i, v) in enumerate(gas_Cantera.net_rates_of_progress)
            @test isapprox(progress_rate(reaction(gas_Apophis, i), in=u"kmol/m^3/s"), v, rtol = 0.005)
        end
    end
end

function test_production_rates_dT(gas_Apophis, gas_Cantera)
    @testset verbose = false "production rates: dT" begin
        for (i, v) in enumerate(gas_Cantera.net_test_production_rates_ddT)
            @test isapprox(production_rate(species(gas_Apophis, i), :dT) * 10^3, v, rtol = 0.005)
        end
    end
end

function test_production_rates_dC(gas_Apophis, gas_Cantera)
    @testset verbose = false "production rates: dC" begin
        dC = gas_Cantera.net_test_production_rates_ddC
        for i in axes(dC, 1), j in axes(dC, 2)
            @test isapprox(production_rate(species(gas_Apophis, i), :dC)[i, j] * 10^3, dC[i, j], rtol = 0.005)
        end
    end
end

function test_species_calculations(gas_Apophis, gas_Cantera, T)
    OneToK = eachindex(species(gas_Apophis))
    @testset verbose = true "Values" begin
        @testset "enthalpies" begin for k in OneToK test_enthalpies(species(gas_Apophis, k), gas_Cantera.species(k-1), T) end end 
        @testset "entropies" begin for k in OneToK test_entropies(species(gas_Apophis, k), gas_Cantera.species(k-1), T) end end 
        @testset "capacities_pressure" begin for k in OneToK test_capacities_pressure(species(gas_Apophis, k), gas_Cantera.species(k-1), T) end end
        test_equilibrium_constants(gas_Apophis, gas_Cantera)
        test_total_molar_concentrations(gas_Apophis, gas_Cantera)
        test_forward_rate_of_progress(gas_Apophis, gas_Cantera)
        #test_forward_rate_constants(gas_Apophis, gas_Cantera)
        #test_reverse_rate_of_progress(gas_Apophis, gas_Cantera)
        #test_reverse_rate_constants(gas_Apophis, gas_Cantera)
        #test_progress_rates(gas_Apophis, gas_Cantera)
        #test_production_rates(gas_Apophis, gas_Cantera)
    end
end

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
            test_species_calculations(gas_Apophis, gas_Cantera, T)
            #test_species_derivatives(gas_Apophis, gas_Cantera)
        end
    end
end
