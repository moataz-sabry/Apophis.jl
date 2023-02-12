########################################################################  Species  ########################################################################

function test_species_enthalpies(gas_Apophis, gas_Cantera)
    species_Apophis_enthalpies = enthalpies(gas_Apophis)
    species_Cantera_enthalpies = gas_Cantera.partial_molar_enthalpies
    for k in eachindex(species_Apophis_enthalpies)
        @test species_Apophis_enthalpies[k] ≈ species_Cantera_enthalpies[k]
    end
end

function test_species_entropies(gas_Apophis, gas_Cantera) ## gas_Cantera.partial_molar_entropies; Apparently, only for 300K
    species_Apophis_entropies = entropies(gas_Apophis)
    for k in eachindex(species_Apophis_entropies)
        species_Cantera_entropy = gas_Cantera.species(k-1).thermo.s(gas_Cantera.T)
        @test species_Apophis_entropies[k] ≈ species_Cantera_entropy
    end
end

function test_species_heat_capacities_pressure(gas_Apophis, gas_Cantera)
    species_Apophis_heat_capacity_pressure = heat_capacities_pressure(gas_Apophis)
    species_Cantera_heat_capacity_pressure = gas_Cantera.partial_molar_cp
    for k in eachindex(species_Apophis_heat_capacity_pressure)
        @test species_Apophis_heat_capacity_pressure[k] ≈ species_Cantera_heat_capacity_pressure[k]
    end
end

function test_species_production_rates(gas_Apophis, gas_Cantera)
    species_production_rates_Apophis = production_rates(gas_Apophis)
    species_production_rates_Cantera = gas_Cantera.net_production_rates
    for k in eachindex(species_production_rates_Apophis)
        @test species_production_rates_Apophis[k] ≈ species_production_rates_Cantera[k] rtol = 5e-2
    end
end

########################################################################  Reactions  ########################################################################

function test_reaction_equilibrium_constants(gas_Apophis, gas_Cantera)
    reaction_equilibrium_constants = gas_Cantera.equilibrium_constants
    for i in eachindex(reaction_equilibrium_constants)
        reaction_equilibrium_constant_Cantera = reaction_equilibrium_constants[i]
        reaction_equilibrium_constant_Apophis = Apophis.equilibrium_constants(Apophis.reaction(gas_Apophis, i), gas_Cantera.T)
        @test reaction_equilibrium_constant_Apophis ≈ reaction_equilibrium_constant_Cantera rtol = 5e-2
    end
end
function total_molar_concentrations(gas_Apophis, gas_Cantera)
    reactions_third_body_concentrations_Cantera = gas_Cantera.third_body_concentrations
    for i in eachindex(gas_Cantera.third_body_concentrations)
        reaction_third_body_concentrations_Cantera = reactions_third_body_concentrations_Cantera[i]
        reaction = Apophis.reaction(gas_Apophis, i)
        if !isnan(reaction_third_body_concentrations_Cantera)
            reaction_third_body_concentrations_Apophis = Apophis.total_molar_concentration(molar_concentrations(gas_Apophis), reaction.enhancement_factors)
            @test reaction_third_body_concentrations_Apophis ≈ reaction_third_body_concentrations_Cantera rtol = 5e-2
        end
    end
end

function test_reaction_forward_rate_constants(gas_Apophis, gas_Cantera)   
    for (i, forward_rate_constant_Cantera) in enumerate(gas_Cantera.forward_rate_constants)
        reaction = Apophis.reaction(gas_Apophis, i)
        M = reaction isa Apophis.ThreeBodyReaction ? Apophis.total_molar_concentration(molar_concentrations(gas_Apophis), reaction.enhancement_factors) : one(Float64)
        @test M * Apophis.forward_rate(reaction) ≈ forward_rate_constant_Cantera rtol = 5e-2
    end 
end

function test_reaction_forward_rate_of_progress(gas_Apophis, gas_Cantera)   
    for (i, forward_rate_constant_Cantera) in enumerate(gas_Cantera.forward_rates_of_progress)
        reaction = Apophis.reaction(gas_Apophis, i)
        Π = Apophis.step(reaction.reactants, molar_concentrations(gas_Apophis))
        M = reaction isa Apophis.ThreeBodyReaction ? Apophis.total_molar_concentration(molar_concentrations(gas_Apophis), reaction.enhancement_factors) : one(Float64)
        @test M * Π * Apophis.forward_rate(reaction) ≈ forward_rate_constant_Cantera rtol = 5e-2
    end
end

function test_reaction_reverse_rate_constants(gas_Apophis, gas_Cantera)
    for (i, reverse_rate_constant_Cantera) in enumerate(gas_Cantera.reverse_rate_constants)
        reaction = Apophis.reaction(gas_Apophis, i)
        M = reaction isa Apophis.ThreeBodyReaction ? Apophis.total_molar_concentration(molar_concentrations(gas_Apophis), reaction.enhancement_factors) : one(Float64)
        @test M * Apophis.reverse_rate(reaction) ≈ reverse_rate_constant_Cantera rtol = 5e-2
    end
end

function test_reaction_reverse_rate_of_progress(gas_Apophis, gas_Cantera)
    for (i, reverse_rate_constant_Cantera) in enumerate(gas_Cantera.reverse_rates_of_progress)
        reaction = Apophis.reaction(gas_Apophis, i)
        Π = Apophis.step(reaction.products, molar_concentrations(gas_Apophis))
        M = reaction isa Apophis.ThreeBodyReaction ? Apophis.total_molar_concentration(molar_concentrations(gas_Apophis), reaction.enhancement_factors) : one(Float64)
        @test M * Π * Apophis.reverse_rate(reaction) ≈ reverse_rate_constant_Cantera rtol = 5e-2
    end
end

_test_reaction_progress_rates(progress_rate_Apophis, progress_rate_Cantera) = @test progress_rate_Apophis ≈ progress_rate_Cantera rtol = 5e-2
test_reaction_progress_rates(gas_Apophis, gas_Cantera) = foreach(_test_reaction_progress_rates, progress_rates(gas_Apophis), gas_Cantera.net_rates_of_progress)

function test_species_values(gas_Apophis, gas_Cantera)
    @testset "– Enthalpy" test_species_enthalpies(gas_Apophis, gas_Cantera)
    @testset "– Entropy" test_species_entropies(gas_Apophis, gas_Cantera)
    @testset "– Heat Capacity" test_species_heat_capacities_pressure(gas_Apophis, gas_Cantera)
    @testset "– Production Rate" test_species_production_rates(gas_Apophis, gas_Cantera)
end

function test_reaction_values(gas_Apophis, gas_Cantera)
    @testset "– Equilibrium Constant" test_reaction_equilibrium_constants(gas_Apophis, gas_Cantera)
    @testset "– Total Molar Conc." total_molar_concentrations(gas_Apophis, gas_Cantera)
    @testset "– Forward Rate" test_reaction_forward_rate_constants(gas_Apophis, gas_Cantera)
    @testset "– Forward Progress Rate" test_reaction_forward_rate_of_progress(gas_Apophis, gas_Cantera)
    @testset "– Reverse Rate" test_reaction_reverse_rate_constants(gas_Apophis, gas_Cantera)
    @testset "– Reverse Progress Rate" test_reaction_reverse_rate_of_progress(gas_Apophis, gas_Cantera)
    @testset "– Progress Rate" test_reaction_progress_rates(gas_Apophis, gas_Cantera)
end

function test_values(mech::Union{String, Symbol})
    gas_Apophis = Gas(mech)
    gas_Cantera = run_cantera(mech)
    for _ in 1:rand(1:3)
        T = rand(300.0:3000.0)
        P = rand(0.5Apophis.Pa:2Apophis.Pa)

        rnd = (rand ∘ length ∘ species)(gas_Apophis)
        Y = rnd / sum(rnd)
        
        TPY!(gas_Apophis, T, P, Y) |> update
        gas_Cantera.TPY = T, P, Y

        # Test the species values
        @testset "Species" test_species_values(gas_Apophis, gas_Cantera)
        
        # Test the reactions values
        @testset "Reactions" test_reaction_values(gas_Apophis, gas_Cantera)
    end
end
