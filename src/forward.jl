function init(mech::Symbol, temperature::K, pressure::K=Pₐ; s...) where {K<:Float64}

    readmechanism(mech)
    gasexpr = :(gas = Gas{$K}(:H2, mechanism))
    eval(gasexpr)

    species = String.(first.(collect(s)))
    fractions = last.(collect(s)) / 100
    ∑fractions = sum(fractions)

    indicies = indexin(species, mechanism.species)
    in(nothing, indicies) && error("one or more of species are not part of $(mech) mechanism")
    ∑fractions == one(K) || @warn "∑Y = $(∑fractions) ≠ 1!"

    W⁻¹ = mechanism.inverse_molecular_weight

    ## integrate with concentrations? ##
    T = gas.initial.temperature = temperature
    P = gas.initial.pressure = pressure

    Y = gas.initial.mass_fractions
    X = gas.initial.molar_fractions
    C = gas.initial.molar_concentration

    for i in eachindex(indicies)
        j = indicies[i]
        Y[j] = fractions[i]
    end

    W̅ = inv(Y ⋅ W⁻¹)
    ρ = gas.initial.density = P * W̅ / (R * T)

    for i in eachindex(mechanism.species)
        YW⁻¹ = Y[i] * W⁻¹[i]
        C[i] = YW⁻¹ * ρ
        X[i] = YW⁻¹ * W̅
    end

    return nothing
end

function init(mech::Symbol, temperature::K, pressure::K, mass_fractions) where {K<:Num}

    readmechanism(mech)
    gasexpr = :(gas = Gas{$K}(:H2, mechanism))
    eval(gasexpr)

    W⁻¹ = mechanism.inverse_molecular_weight

    ## integrate with concentrations? ##
    T = gas.initial.temperature = temperature
    P = gas.initial.pressure = pressure

    Y = gas.initial.mass_fractions = mass_fractions
    X = gas.initial.molar_fractions
    C = gas.initial.molar_concentration

    W̅ = inv(Y ⋅ W⁻¹)
    ρ = gas.initial.density = P * W̅ / (R * T)

    for i in eachindex(mechanism.species)
        YW⁻¹ = Y[i] * W⁻¹[i]
        C[i] = YW⁻¹ * ρ
        X[i] = YW⁻¹ * W̅
    end

    return nothing
end

function polynomials((; current, intermediate, mechanism)::Gas) ## computes nasa polynomials for species entropy, entahlpy, ...

    Tₖ = mechanism.common_temperature
    A = mechanism.upper_temperature_coefficients
    a = mechanism.lower_temperature_coefficients

    coeffs = intermediate.polynomial_coefficients
    a₁, a₂, a₃, a₄, a₅, a₆, a₇ = (view(coeffs, a, :) for a in 1:7) ## multiple assign?

    cₚ = intermediate.heat_capacity_pressure
    cᵥ = intermediate.heat_capacity_volume
    h = intermediate.enthalpy_species
    u = intermediate.internal_energy
    s = intermediate.entropy_species

    T = current.temperature

    for i in eachindex(mechanism.species) ## 500 ns; faster than coeffs[:, i] .= view(A, :, i); 690 ns for GRI3; eachcol?
        offset = 7i - 6  ## 7(i - 1) + 1
        #if T ≥ Tₖ[i]
        copyto!(coeffs, offset, A, offset, 7)
        #else
        #    copyto!(coeffs, offset, a, offset, 7)
        #end
    end

    for i in eachindex(mechanism.species) ## save repeated operations
        cₚ[i] = (a₁[i] + a₂[i] * T + a₃[i] * T^2 + a₄[i] * T^3 + a₅[i] * T^4)R
        cᵥ[i] = cₚ[i] - R

        h[i] = (a₁[i] * T + a₂[i] / 2 * T^2 + a₃[i] / 3 * T^3 + a₄[i] / 4 * T^4 + a₅[i] / 5 * T^5 + a₆[i])R
        u[i] = h[i] - (R * T)

        s[i] = (a₁[i] * log(T) + a₂[i] * T + a₃[i] / 2 * T^2 + a₄[i] / 3 * T^3 + a₅[i] / 4 * T^4 + a₇[i])R
    end

    return nothing
end

function concentrations((; current, intermediate, mechanism)::Gas) ##computes species properties

    W⁻¹ = mechanism.inverse_molecular_weight
    α = mechanism.enhancement_factors
    M = intermediate.total_molar_concentration

    Y = current.mass_fractions
    X = current.molar_fractions
    C = current.molar_concentration
    ρ = current.density

    W̅ = inv(Y ⋅ W⁻¹)
    for i in eachindex(mechanism.species)
        YW⁻¹ = Y[i] * W⁻¹[i]
        C[i] = YW⁻¹ * ρ
        X[i] = YW⁻¹ * W̅
    end

    mul!(M, α, C)
    return nothing
end

function reactionconstants((; current, intermediate, mechanism)::Gas)
    ## optimizable

    T = current.temperature
    RT⁻¹ = inv(R * T)
    RcT⁻¹ = inv(Rc * T)

    nr = length(mechanism.reactions)
    nt = length(mechanism.threebody_reactions) ##optimize
    nf = length(mechanism.falloff_reactions)
    ny = nr - nf

    A, b, E = (view(mechanism.high_pressure_parameters, 1:ny, p) for p in 1:3)
    A∞, b∞, E∞ = (view(mechanism.high_pressure_parameters, ny+1:nr, p) for p in 1:3)
    Aₒ, bₒ, Eₒ = (view(mechanism.low_pressure_parameters, :, p) for p in 1:3)
    Aᵣ, bᵣ, Eᵣ = (view(mechanism.reverse_reaction_parameters, :, p) for p in 1:3)

    ν = mechanism.stoichiometric_reactions
    ∑ν = mechanism.stoichiometric_sum
    a, T₃, T₁, T₂ = (view(mechanism.troe_parameters, :, p) for p in 1:4)

    h = intermediate.enthalpy_species
    s = intermediate.entropy_species
    M = view(intermediate.total_molar_concentration, nt+1:nt+nf)
    H = intermediate.enthalpy_reactions
    S = intermediate.entropy_reactions
    kf = intermediate.forward_rate_constant
    kr = intermediate.reverse_rate_constant

    T = current.temperature

    for k in eachindex(A)
        kf[k] = A[k] * T^b[k] * exp(-E[k] * RcT⁻¹)
    end

    for f in eachindex(mechanism.falloff_reactions)
        kₒ = Aₒ[f] * T^bₒ[f] * exp(-Eₒ[f] * RcT⁻¹)
        k∞ = A∞[f] * T^b∞[f] * exp(-E∞[f] * RcT⁻¹)

        Pᵣ = kₒ * M[f] / k∞
        log10Pᵣ = log10(Pᵣ)

        Fc = (1.0 - a[f]) * exp(-T / T₃[f]) + a[f] * exp(-T / T₁[f]) + exp(-T₂[f] / T)
        log10Fc = log10(Fc)

        c = -0.4 - 0.67 * log10Fc
        n = 0.75 - 1.27 * log10Fc
        F = 10.0^(log10Fc / (1.0 + ((log10Pᵣ + c) / (n - d * (log10Pᵣ + c)))^2.0))

        kf[ny+f] = F * k∞ * Pᵣ / (1.0 + Pᵣ)
    end

    mul!(H, ν, h)
    mul!(S, ν, s)

    for k in mechanism.reversible_equilibrium
        Kp = exp((S[k] * T - H[k]) * RT⁻¹)
        Kc = Kp * (Pₐ * RT⁻¹)^∑ν[k]
        kr[k] = kf[k] / Kc
        kr[k] = Aᵣ[k] * T^bᵣ[k] * exp(-Eᵣ[k] * RcT⁻¹)
    end


    return nothing
end

function rates((; current, intermediate, mechanism)::Gas)

    nt = length(mechanism.threebody_reactions)
    ne = length(mechanism.elementary_reactions)
    threebody_range = ne+1:ne+nt

    W⁻¹ = mechanism.inverse_molecular_weight
    Y = current.mass_fractions
    X = current.molar_fractions
    reactantsindices = mechanism.reactants_indices
    productsindices = mechanism.products_indices
    C = current.molar_concentration
    W = mechanism.molecular_weight

    nur = mechanism.stoichiometric_reactants
    nup = mechanism.stoichiometric_products
    nut = mechanism.stoichiometric_transpose

    cᵥ = intermediate.heat_capacity_volume
    u = intermediate.internal_energy
    kf = intermediate.forward_rate_constant
    kr = intermediate.reverse_rate_constant
    M = view(intermediate.total_molar_concentration, 1:nt)
    q = intermediate.rate_of_progress
    ω̇ = intermediate.production_rate
    Ẏ = intermediate.mass_change_rate
    Ṫ = intermediate.temperature_change_rate

    ρ = current.density

    @inbounds for ci in reactantsindices ## combine into one function
        i, j = Tuple(ci)
        kf[i] *= C[j]^nur[ci]
    end

    @inbounds for ci in productsindices
        i, j = Tuple(ci)
        kr[i] *= C[j]^nup[ci]
    end

    for k in eachindex(q)
        q[k] = kf[k] - kr[k]
        k ∈ threebody_range && (q[k] *= M[k-ne])
    end

    mul!(ω̇, nut, q)

    W̅ = inv(Y ⋅ W⁻¹)
    c̅ᵥ = cᵥ ⋅ X / W̅

    for i in eachindex(Ẏ)
        Ẏ[i] = ω̇[i] * W[i] / ρ
    end

    Ṫ[1] = -(ω̇ ⋅ u) / (ρ * c̅ᵥ) ## worth it?

    return nothing
end

function step!(gas::Gas{K}, Y::AbstractVector{K}, T::K) where {K<:Real} ### AbstractVector allocates +1

    #gas.current = gas.initial
    gas.current.mass_fractions = Y
    gas.current.temperature = T
    gas.current.density = gas.initial.density ## in init?

    polynomials(gas)
    concentrations(gas)
    reactionconstants(gas)
    rates(gas)

    return nothing
end

function IdealGasReactor!(du, u, p, t) #DGL

    (; mechanism, intermediate) = p
    ns = eachindex(mechanism.species)

    Y = view(u, ns)
    T = last(u)
    step!(gas, Y, T)

    Ẏ = intermediate.mass_change_rate
    Ṫ = only(intermediate.temperature_change_rate)

    du[ns] = Ẏ
    du[end] = Ṫ

    return nothing
end

function simulate(t, gas::Gas; maxis=1e5,
    abs::T=1e-10, rel::T=1e-10) where {T<:Real}

    span = (0.0, T(t))
    p = gas
    uₒ = vcat(gas.initial.mass_fractions, gas.initial.temperature)

    ODE = ODEProblem(IdealGasReactor!, uₒ, span, p)
    solution = solve(ODE, CVODE_BDF(), abstol=abs, reltol=rel, maxiters=Int(maxis))

    return solution
end
