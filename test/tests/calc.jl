########################################################################  Species  ########################################################################

function test_species_enthalpies(gas_Apophis, gas_Cantera)
    species_Apophis_enthalpies = enthalpies(gas_Apophis; in=u"J/kmol")
    species_Cantera_enthalpies = gas_Cantera.partial_molar_enthalpies
    for k in eachindex(species_Apophis_enthalpies)
        try 
            @test species_Apophis_enthalpies[k] ≈ species_Cantera_enthalpies[k]
        catch e
            println("Apophis: $(species_Apophis_enthalpies[k]) & Cantera: $(species_Cantera_enthalpies[k]) for species $(gas_Apophis.species[k].formula) at temperature $T")
            rethrow(e)
        end
    end
end

function test_species_entropies(gas_Apophis, gas_Cantera) ## gas_Cantera.partial_molar_entropies; Apparently, only for 300K
    species_Apophis_entropies = entropies(gas_Apophis; in=u"J/(kmol*K)")
    for k in eachindex(species_Apophis_entropies)
        species_Cantera_entropy = gas_Cantera.species(k-1).thermo.s(gas_Cantera.T)
        try 
            @test species_Apophis_entropies[k] ≈ species_Cantera_entropy
        catch e
            println("Apophis: $(species_Apophis_entropies[k]) & Cantera: $(species_Cantera_entropies[k]) for species $(gas_Apophis.species[k].formula) at temperature $T")
            rethrow(e)
        end
    end
end

function test_species_heat_capacities_pressure(gas_Apophis, gas_Cantera)
    species_Apophis_heat_capacity_pressure = heat_capacities_pressure(gas_Apophis; in=u"J/(kmol*K)")
    species_Cantera_heat_capacity_pressure = gas_Cantera.partial_molar_cp
    for k in eachindex(species_Apophis_heat_capacity_pressure)
        try
            @test species_Apophis_heat_capacity_pressure[k] ≈ species_Cantera_heat_capacity_pressure[k]
        catch e
            println("Apophis: $(species_Apophis_heat_capacity_pressure[k]) & Cantera: $(species_Cantera_heat_capacity_pressure[k]) for species $(gas_Apophis.species[k].formula) at temperature $T")
            rethrow(e)
        end
    end
end

function test_species_production_rates(gas_Apophis, gas_Cantera)
    species_production_rates_Apophis = production_rates(gas_Apophis, in=u"kmol/m^3/s")
    species_production_rates_Cantera = gas_Cantera.net_production_rates
    for k in eachindex(species_production_rates_Apophis)
        try
            @test isapprox(species_production_rates_Apophis[k], species_production_rates_Cantera[k], rtol = 0.005)
        catch e
            print("Apophis: $(species_production_rates_Apophis[k]) & Cantera: $(species_production_rates_Cantera[k]) for species $(gas_Apophis.species[k].formula) at temperature $T")
            rethrow(e)
        end
    end
end

########################################################################  Reactions  ########################################################################

function test_equilibrium_constants(gas_Apophis, gas_Cantera)
    reaction_equilibrium_constants = gas_Cantera.equilibrium_constants
    for i in eachindex(reaction_equilibrium_constants)
        reaction_equilibrium_constant_Cantera = reaction_equilibrium_constants[i]
        reaction_Apophis = reaction(gas_Apophis, i)
        ∑v = reaction_Apophis.reaction_order

        reaction_equilibrium_constant_Apophis = Apophis.equilibrium_constants(reaction_Apophis, gas_Cantera.T) * ustrip(uparse("(kmol/m^3)^$∑v/s"), 1uparse("(mol/cm^3)^$∑v/s"))
        @test reaction_equilibrium_constant_Apophis ≈ reaction_equilibrium_constant_Cantera
    end
end

function test_total_molar_concentrations(gas_Apophis, gas_Cantera)
    for (i, v) in enumerate(gas_Cantera.third_body_concentrations)
        if !isnan(v)
            reaction_Apophis = reaction(gas_Apophis, i)
            M = Apophis.total_molar_concentration(gas_Apophis.state.C, reaction_Apophis.enhancement_factors) * ustrip(uparse("kmol/m^3"), 1uparse("mol/cm^3"))
            @test M ≈ v rtol = 0.0001
        end
    end
end

function test_forward_rate_constants(gas_Apophis, gas_Cantera)   
    @testset verbose = false "forward rate constants" begin 
        for (i, v) in enumerate(gas_Cantera.forward_rate_constants)
            reaction_Apophis = gas_Apophis.mechanism.reactions[k]
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
            ∑νᵣ = sum(v -> last(v) |> abs, gas_Apophis.mechanism.reactions[k].reactants) + (gas_Apophis.mechanism.reactions[k] isa Apophis.ElementaryReaction ? 0.0 : 1.0)
            unitifykf = 10.0^3(1-∑νᵣ)
            n = gas_Apophis.mechanism.reactions[k].reaction_order
            unitify = 10.0^-3n
            unitifyM = 10.0^3
            M = gas_Apophis.mechanism.reactions[k] isa Apophis.ThreeBodyReaction ? unitifyM * Apophis.total_molar_concentration(gas_Apophis.state.C, gas_Apophis.mechanism.reactions[k].enhancement_factors) : one(Float64)
            @test gas_Apophis.mechanism.reactions[k].rates.kr.val[] * unitifykf * unitify * M ≈ v rtol = 0.005
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

function test_species_values(gas_Apophis, gas_Cantera)
    @testset "– Enthalpy" test_species_enthalpies(gas_Apophis, gas_Cantera)
    @testset "– Entropy" test_species_entropies(gas_Apophis, gas_Cantera)
    @testset "– Heat Capacity" test_species_heat_capacities_pressure(gas_Apophis, gas_Cantera)
    @testset "– Production Rate" test_species_production_rates(gas_Apophis, gas_Cantera)
end

function test_reaction_values(gas_Apophis, gas_Cantera)
    @testset "– Equilibrium Constant" test_equilibrium_constants(gas_Apophis, gas_Cantera)
    @testset "– Total Molar Concentrations" test_total_molar_concentrations(gas_Apophis, gas_Cantera)

end

function test_values(mech::Union{String, Symbol})
    gas_Apophis = Gas(mech)
    gas_Cantera = run_cantera(mech)
    for _ in 1:1
        T = rand(300.0:3000.0)
        P = rand(0.5Apophis.Pa:2Apophis.Pa)

        rnd = (rand ∘ length ∘ species)(gas_Apophis)
        Y = rnd / sum(rnd)
        
        TPY!(gas_Apophis, T, P, Y) |> update
        gas_Cantera.TPY = T, P * ustrip(u"Pa", 1u"dyn/cm^2"), Y

        @testset "Species" test_species_values(gas_Apophis, gas_Cantera)
        @testset "Reactions" test_reaction_values(gas_Apophis, gas_Cantera)
    end
end
