ct = pyimport("cantera")
function run_cantera(mech::Symbol)
    path = pkgdir(Apophis, "test/mechanisms/$mech/$mech.yaml")
    gas = ct.Solution(path)
    return gas
end

########################################################################### Reader ##########################################################################################
#############################################################################################################################################################################

test_total_species(gas_Apophis, gas_Cantera) = @test length(species(gas_Apophis)) == length(gas_Cantera.species())
test_species_name(species_Apophis, species_Cantera) = @test species_Apophis.formula == Symbol(species_Cantera.name)
test_species_max_temp(species_Apophis, species_Cantera) = @test species_Apophis.nasa_polynomial.Tmax == species_Cantera.thermo.max_temp
test_species_min_temp(species_Apophis, species_Cantera) = @test species_Apophis.nasa_polynomial.Tmin == species_Cantera.thermo.min_temp
test_species_coeffs(species_Apophis, species_Cantera) = @test [species_Apophis.nasa_polynomial.A..., species_Apophis.nasa_polynomial.a...] == species_Cantera.thermo.coeffs[2:end]

function test_species_mid_temp(species_Apophis, species_Cantera)
    np = species_Apophis.nasa_polynomial
    Tₘ = isequal(np.A, np.a) ? np.Tmax : np.Tₘ ## Because Cantera works like that, if a species has same coefficients for high and low intervals.
    @test Tₘ == species_Cantera.thermo.coeffs[1]
end

function test_species_components(species_Apophis, species_Cantera)
    for (i, v) in species_Apophis.components
        key = String(i) |> lowercase |> uppercasefirst
        @test species_Cantera.composition[key] == v
    end
end

function test_species_molecular_weights(gas_Apophis, gas_Cantera)
    for (mA, mC) in zip(molecular_weights(gas_Apophis), gas_Cantera.molecular_weights)
        @test mA ≈ mC rtol = 0.001
    end
end 

function test_species_reader(gas_Apophis, gas_Cantera)
    OneToK = eachindex(species(gas_Apophis))
    @testset verbose = true "Species" begin
        @testset "No. Species = $(length(OneToK))" test_total_species(gas_Apophis, gas_Cantera)
        @testset "species name" begin for k in OneToK test_species_name(species(gas_Apophis, k), gas_Cantera.species(k-1)) end end 
        @testset "components" begin for k in OneToK test_species_components(species(gas_Apophis, k), gas_Cantera.species(k-1)) end end
        @testset "molecilar weights" test_species_molecular_weights(gas_Apophis, gas_Cantera)
        @testset "max. temperature" begin for k in OneToK test_species_max_temp(species(gas_Apophis, k), gas_Cantera.species(k-1)) end end 
        @testset "min. temperature" begin for k in OneToK test_species_min_temp(species(gas_Apophis, k), gas_Cantera.species(k-1)) end end 
        @testset "mid. temperature" begin for k in OneToK test_species_mid_temp(species(gas_Apophis, k), gas_Cantera.species(k-1)) end end 
        @testset "coefficients" begin for k in OneToK test_species_coeffs(species(gas_Apophis, k), gas_Cantera.species(k-1)) end end 
    end
end

test_total_reactions(gas_Apophis, gas_Cantera) = @test length(reactions(gas_Apophis)) == length(gas_Cantera.reactions())
test_reaction_reversibility(reaction_Apophis, reaction_Cantera) = @test reaction_Apophis.isreversible == reaction_Cantera.reversible

adjust_type(::Apophis.ElementaryReaction) = "reaction"
adjust_type(::Apophis.ThreeBodyReaction) = "three-body"
adjust_type(::Apophis.FallOffReaction) = "falloff"

test_reaction_type(reaction_Apophis, reaction_Cantera) = @test adjust_type(reaction_Apophis) == reaction_Cantera.reaction_type

function test_reaction_components(reaction_Apophis, reaction_Cantera)
    # println(reaction_Apophis, " - ", reaction_Cantera)
    # println(reaction_Apophis.reactants, " - ", reaction_Cantera.reactants)
    # println(reaction_Apophis.products, " - ", reaction_Cantera.products)
    for (i, v) in reaction_Apophis.reactants
        @test reaction_Cantera.reactants[String(i.formula)] == abs(v) ## Cantera saves ν as +ve regardless of being for a reactant or product.
    end
    for (i, v) in reaction_Apophis.products
        @test reaction_Cantera.products[String(i.formula)] == v
    end
end

function test_reaction_rates(reaction_Apophis, reaction_Cantera)
    ∑νᵣ = sum(v -> last(v) |> abs, reaction_Apophis.reactants)# + 1 #??
    unitify = 10.0^3(1-∑νᵣ)
    if reaction_Apophis isa Apophis.FallOffReaction
        #@test reaction_Apophis.high_pressure_parameters.A * unitify *10^-3 ≈ reaction_Cantera.rate.high_rate.pre_exponential_factor
        @test reaction_Apophis.high_pressure_parameters.β == reaction_Cantera.rate.high_rate.temperature_exponent
        @test reaction_Apophis.high_pressure_parameters.E * 4184 ≈ reaction_Cantera.rate.high_rate.activation_energy ## cal to joule
        
        #@test reaction_Apophis.low_pressure_parameters.A * unitify *10^-3 ≈ reaction_Cantera.rate.low_rate.pre_exponential_factor
        @test reaction_Apophis.low_pressure_parameters.β == reaction_Cantera.rate.low_rate.temperature_exponent
        @test reaction_Apophis.low_pressure_parameters.E * 4184 ≈ reaction_Cantera.rate.low_rate.activation_energy
    else reaction_Apophis isa Apophis.ElementaryReaction
        #@test reaction_Apophis.forward_rate_parameters.A * unitify ≈ reaction_Cantera.rate.pre_exponential_factor
        @test reaction_Apophis.forward_rate_parameters.β == reaction_Cantera.rate.temperature_exponent
        @test reaction_Apophis.forward_rate_parameters.E * 4184 ≈ reaction_Cantera.rate.activation_energy
    end
end

function test_reaction_troe(reaction_Apophis, reaction_Cantera)
    if reaction_Apophis isa Apophis.FallOffReaction && !isnothing(reaction_Apophis.troe_parameters)
        #@test reaction_Apophis.high_pressure_parameters.A * unitify *10^-3 ≈ reaction_Cantera.rate.high_rate.pre_exponential_factor
        @test reaction_Apophis.troe_parameters.a ≈ reaction_Cantera.rate.falloff_coeffs[1] ## cal to joule
        @test reaction_Apophis.troe_parameters.T₃ ≈ reaction_Cantera.rate.falloff_coeffs[2] ## cal to joule
        @test reaction_Apophis.troe_parameters.T₁ ≈ reaction_Cantera.rate.falloff_coeffs[3] ## cal to joule
        if iszero(reaction_Apophis.troe_parameters.T₂)
            @test length(reaction_Cantera.rate.falloff_coeffs) == 3
        else
            @test reaction_Apophis.troe_parameters.T₂ ≈ reaction_Cantera.rate.falloff_coeffs[4] ## cal to joule
        end
    end
end

function test_reactions_reader(gas_Apophis, gas_Cantera)
    OneToI = eachindex(reactions(gas_Apophis))
    @testset verbose = true "Reactions" begin 
        @testset "No. Reactions = $(length(OneToI))" test_total_reactions(gas_Apophis, gas_Cantera)
        @testset "reversibility" begin for i in OneToI test_reaction_reversibility(reaction(gas_Apophis, i), gas_Cantera.reaction(i-1)) end end
        @testset "reaction type" begin for i in OneToI test_reaction_type(reaction(gas_Apophis, i), gas_Cantera.reaction(i-1)) end end
        @testset "reactants & products" begin for i in OneToI test_reaction_components(reaction(gas_Apophis, i), gas_Cantera.reaction(i-1)) end end
        #@testset "Arrhenius parameters" begin for i in OneToI test_reaction_rates(reaction(gas_Apophis, i), gas_Cantera.reaction(i-1)) end end
        #@testset "Troe parameters" begin for i in OneToI test_reaction_troe(reaction(gas_Apophis, i), gas_Cantera.reaction(i-1)) end end
    end
end

function test_reader(mech::Symbol)
    @testset verbose = true "Mechanism: $mech" begin
        gas_Apophis = Gas(mech)
        gas_Cantera = run_cantera(mech)
        test_species_reader(gas_Apophis, gas_Cantera)
        test_reactions_reader(gas_Apophis, gas_Cantera)
    end
end

#################################################################### Calculations: values ###################################################################################
#############################################################################################################################################################################

test_enthalpies(species_Apophis, species_Cantera, T) = @test enthalpy(species_Apophis; in=u"J/(kmol*K)") ≈ species_Cantera.thermo.h(T) # partial_molar_test_enthalpies
test_entropies(species_Apophis, species_Cantera, T) = @test entropy(species_Apophis; in=u"J/kmol") ≈ species_Cantera.thermo.s(T)
test_capacities_pressure(species_Apophis, species_Cantera, T) = @test heat_capacity_pressure(species_Apophis; in=u"J/kmol") ≈ species_Cantera.thermo.cp(T)

function test_production_rates(gas_Apophis, gas_Cantera)
    @testset "production rates" begin
        for (i, v) in enumerate(gas_Cantera.net_production_rates)
            @test isapprox(Apophis.production_rates(gas_Apophis, in=u"kmol/m^3/s")[i], v, rtol = 0.005)
        end
    end
end

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