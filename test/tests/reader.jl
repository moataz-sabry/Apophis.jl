########################################################################  Species  ########################################################################

function test_total_species(gas_Apophis, gas_Cantera)
    total_species_Apophis = species(gas_Apophis) |> length
    total_species_Cantera = gas_Cantera.n_total_species
    @test total_species_Apophis == total_species_Cantera
end

function _test_species_name(species_Apophis, species_Cantera)
    species_name_Apophis = species_Apophis.formula
    species_name_Cantera = Symbol(species_Cantera.name)
    @test species_name_Apophis == species_name_Cantera
end

test_species_name(gas_Apophis, gas_Cantera) = foreach(_test_species_name, species(gas_Apophis), gas_Cantera.species())

function test_species_molecular_weights(gas_Apophis, gas_Cantera)
    molecular_weights_Apophis = molecular_weights(gas_Apophis)
    molecular_weights_Cantera = gas_Cantera.molecular_weights
    for k in eachindex(molecular_weights_Apophis)
        try
            @test molecular_weights_Apophis[k] ≈ molecular_weights_Cantera[k] rtol = 0.001
        catch e
            println("Molecular weight: Apophis ($molecular_weight_Apophis) & Cantera ($molecular_weight_Cantera) for species $(species(gas_Apophis, k))")
            rethrow(e)
        end
    end
end

function _test_species_components(species_Apophis, species_Cantera)
    for (key, value) in species_Apophis.components
        species_component_Apophis = value
        species_component_Cantera = species_Cantera.composition[(uppercasefirst ∘ lowercase ∘ String)(key)]
        try
            @test species_component_Cantera == species_component_Apophis
        catch e
            println("Composition: Apophis ($value) & Cantera ($(species_Cantera.composition[key])) for species $(species_Apophis.formula), component $key")
            rethrow(e)
        end
    end
end

test_species_components(gas_Apophis, gas_Cantera) = foreach(_test_species_components, species(gas_Apophis), gas_Cantera.species())

function _test_species_max_temp(species_Apophis, species_Cantera)
    species_max_temp_Apophis = species_Apophis.nasa_polynomial.Tmax
    species_max_temp_Cantera = species_Cantera.thermo.max_temp
    try
        @test species_max_temp_Apophis == species_max_temp_Cantera
    catch e
        println("Max temperature: Apophis ($species_max_temp_Apophis) & Cantera ($species_max_temp_Cantera) for species $(species_Apophis.formula)")
        rethrow(e)
    end
end

test_species_max_temp(gas_Apophis, gas_Cantera) = foreach(_test_species_max_temp, species(gas_Apophis), gas_Cantera.species())

function _test_species_min_temp(species_Apophis, species_Cantera)
    species_min_temp_Apophis = species_Apophis.nasa_polynomial.Tmin
    species_min_temp_Cantera = species_Cantera.thermo.min_temp
    try
        @test species_min_temp_Apophis == species_min_temp_Cantera
    catch e
        println("Min temperature: Apophis ($species_min_temp_Apophis) & Cantera ($species_min_temp_Cantera) for species $(species_Apophis.formula)")
        rethrow(e)
    end
end

test_species_min_temp(gas_Apophis, gas_Cantera) = foreach(_test_species_min_temp, species(gas_Apophis), gas_Cantera.species())

function _test_species_coeffs(species_Apophis, species_Cantera)
    species_coeffs_Apophis = [species_Apophis.nasa_polynomial.A..., species_Apophis.nasa_polynomial.a...]
    species_coeffs_Cantera = species_Cantera.thermo.coeffs[2:end]
    for i in eachindex(species_coeffs_Apophis)
        try
            @test species_coeffs_Apophis[i] == species_coeffs_Cantera[i]
        catch e
            println("Coefficient $i: Apophis ($(species_coeffs_Apophis[i])) & Cantera ($(species_coeffs_Cantera[i])) for species $(species_Apophis.formula)")
            rethrow(e)
        end
    end
end

test_species_coeffs(gas_Apophis, gas_Cantera) = foreach(_test_species_coeffs, species(gas_Apophis), gas_Cantera.species())

# np = species_Apophis.nasa_polynomial
# Tₘ = isequal(np.A, np.a) ? np.Tmax : np.Tₘ ## Because Cantera works like that, if a species has same coefficients for high and low intervals.
# @test Tₘ == species_Cantera.thermo.coeffs[1]

function _test_species_mid_temp(species_Apophis, species_Cantera)
    species_mid_temp_Apophis = species_Apophis.nasa_polynomial.Tₘ
    species_mid_temp_Cantera = species_Cantera.thermo.coeffs[1]
    try
        @test species_mid_temp_Apophis == species_mid_temp_Cantera
    catch e
        error("Mid temperature: Apophis ($species_mid_temp_Apophis) & Cantera ($species_mid_temp_Cantera) for species $(species_Apophis.formula)")
        rethrow(e)
    end
end

test_species_mid_temp(gas_Apophis, gas_Cantera) = foreach(_test_species_mid_temp, species(gas_Apophis), gas_Cantera.species())

########################################################################  Reactions  ########################################################################

function test_total_reactions(gas_Apophis, gas_Cantera)
    total_reactions_Apophis = reactions(gas_Apophis) |> length
    total_reactions_Cantera = gas_Cantera.n_reactions
    @test total_reactions_Apophis == total_reactions_Cantera
end

function _test_reversibility(reaction_Apophis, reaction_Cantera)
    reaction_reversibility_Apophis = reaction_Apophis.isreversible
    reaction_reversibility_Cantera = reaction_Cantera.reversible
    @test reaction_reversibility_Apophis == reaction_reversibility_Cantera
end

test_reversibility(gas_Apophis, gas_Cantera) = foreach(_test_reversibility, reactions(gas_Apophis), gas_Cantera.reactions())

function _test_reaction_type(reaction_Apophis, reaction_Cantera)
    reaction_type_Apophis = adjust_type(reaction_Apophis)
    reaction_type_Cantera = reaction_Cantera.reaction_type
    @test reaction_type_Apophis == reaction_type_Cantera
end

test_reaction_type(gas_Apophis, gas_Cantera) = foreach(_test_reaction_type, reactions(gas_Apophis), gas_Cantera.reactions())

function _test_reaction_reactants(reaction_Apophis, reaction_Cantera)
    for (species, value) in reaction_Apophis.reactants
        reaction_reactant_Apophis = abs(value)
        reaction_reactant_Cantera = reaction_Cantera.reactants[String(species.formula)]
        @test reaction_reactant_Cantera == reaction_reactant_Apophis ## Cantera saves ν as +ve regardless of being for a reactant or product.
    end
end

function _test_reaction_products(reaction_Apophis, reaction_Cantera)
    for (species, value) in reaction_Apophis.products
        reaction_product_Apophis = value
        reaction_product_Cantera = reaction_Cantera.products[String(species.formula)]
        @test reaction_product_Cantera == reaction_product_Apophis
    end
end

test_reaction_reactants(gas_Apophis, gas_Cantera) = foreach(_test_reaction_reactants, reactions(gas_Apophis), gas_Cantera.reactions())
test_reaction_products(gas_Apophis, gas_Cantera) = foreach(_test_reaction_products, reactions(gas_Apophis), gas_Cantera.reactions())

function _test_rate_parameters(reaction_Apophis, reaction_Cantera)
    if reaction_Apophis isa Apophis.FallOffReaction
        @test reaction_Apophis.high_pressure_parameters.A ≈ reaction_Cantera.rate.high_rate.pre_exponential_factor
        @test reaction_Apophis.high_pressure_parameters.β ≈ reaction_Cantera.rate.high_rate.temperature_exponent
        @test reaction_Apophis.high_pressure_parameters.E * 4184 ≈ reaction_Cantera.rate.high_rate.activation_energy ## cal to joule
        
        @test reaction_Apophis.low_pressure_parameters.A ≈ reaction_Cantera.rate.low_rate.pre_exponential_factor
        @test reaction_Apophis.low_pressure_parameters.β ≈ reaction_Cantera.rate.low_rate.temperature_exponent
        @test reaction_Apophis.low_pressure_parameters.E * 4184 ≈ reaction_Cantera.rate.low_rate.activation_energy
    else
        @test reaction_Apophis.forward_rate_parameters.A ≈ reaction_Cantera.rate.pre_exponential_factor
        @test reaction_Apophis.forward_rate_parameters.β ≈ reaction_Cantera.rate.temperature_exponent
        @test reaction_Apophis.forward_rate_parameters.E * 4184 ≈ reaction_Cantera.rate.activation_energy
    end
end

test_rate_parameters(gas_Apophis, gas_Cantera) = foreach(_test_rate_parameters, reactions(gas_Apophis), gas_Cantera.reactions())

function _test_troe_parameters(reaction_Apophis, reaction_Cantera)
    if reaction_Apophis isa Apophis.FallOffReaction && !isnothing(reaction_Apophis.troe_parameters)
        @test reaction_Apophis.troe_parameters.a ≈ reaction_Cantera.rate.falloff_coeffs[1]
        @test reaction_Apophis.troe_parameters.T₃ ≈ reaction_Cantera.rate.falloff_coeffs[2]
        @test reaction_Apophis.troe_parameters.T₁ ≈ reaction_Cantera.rate.falloff_coeffs[3]
        if iszero(reaction_Apophis.troe_parameters.T₂)
            @test length(reaction_Cantera.rate.falloff_coeffs) == 3
        else
            @test reaction_Apophis.troe_parameters.T₂ ≈ reaction_Cantera.rate.falloff_coeffs[4]
        end
    end
end

test_troe_parameters(gas_Apophis, gas_Cantera) = foreach(_test_troe_parameters, reactions(gas_Apophis), gas_Cantera.reactions())

function test_species_reader(gas_Apophis, gas_Cantera)
    @testset "– # Species" test_total_species(gas_Apophis, gas_Cantera)
    @testset "– Name" test_species_name(gas_Apophis, gas_Cantera)
    @testset "– Weight" test_species_molecular_weights(gas_Apophis, gas_Cantera)
    @testset "– Components" test_species_components(gas_Apophis, gas_Cantera)
    @testset "– Max. Temp." test_species_max_temp(gas_Apophis, gas_Cantera)
    @testset "– Min. Temp." test_species_min_temp(gas_Apophis, gas_Cantera)
    @testset "– Mid. Temp." test_species_mid_temp(gas_Apophis, gas_Cantera)
    @testset "– Coefficients" test_species_coeffs(gas_Apophis, gas_Cantera)
end

function test_reactions_reader(gas_Apophis, gas_Cantera)
    @testset "– # Reactions" test_total_reactions(gas_Apophis, gas_Cantera)
    @testset "– Reversibility" test_reversibility(gas_Apophis, gas_Cantera)
    @testset "– Type" test_reaction_type(gas_Apophis, gas_Cantera)
    @testset "– Reactants" test_reaction_reactants(gas_Apophis, gas_Cantera)
    @testset "– Products" test_reaction_products(gas_Apophis, gas_Cantera)
    @testset "– Rate Params." test_rate_parameters(gas_Apophis, gas_Cantera)
    @testset "– Troe Params." test_troe_parameters(gas_Apophis, gas_Cantera)
end

function test_reader(mech::Union{String, Symbol})
    gas_Apophis = Gas(mech)
    gas_Cantera = run_cantera(mech)
    
    # Test the species reader
    @testset "Species" begin
        test_species_reader(gas_Apophis, gas_Cantera)
    end
    
    # Test the reactions reader
    @testset "Reactions" begin
        test_reactions_reader(gas_Apophis, gas_Cantera)
    end
end