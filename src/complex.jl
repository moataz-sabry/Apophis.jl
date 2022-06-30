function interval((; current, intermediate, mechanism)::Gas{K}; forward=true) where {K<:Number}

    Tₖ = mechanism.common_temperature
    A = mechanism.upper_temperature_coefficients
    a = mechanism.lower_temperature_coefficients

    coeffs = intermediate.polynomial_coefficients
    T = only(current.temperature)

    for i in eachindex(mechanism.species) ## 500 ns; faster than coeffs[:, i] .= view(A, :, i); 690 ns for GRI3; eachcol?
        offset = 7i - 6  ## 7(i - 1) + 1
        if real(T) ≥ real(Tₖ[i])
            copyto!(coeffs, offset, A, offset, 7)
        else
            copyto!(coeffs, offset, a, offset, 7)
        end
    end
    return nothing
end

function polynomials((; current, intermediate, mechanism)::Gas{K}; forward=true) where {K<:Number} ## computes nasa polynomials for species entropy, entahlpy, ...

    coeffs = intermediate.polynomial_coefficients
    a₁, a₂, a₃, a₄, a₅, a₆, a₇ = (view(coeffs, a, :) for a in 1:7) ## multiple assign?

    cₚ = first(intermediate.heat_capacity_pressure)
    cᵥ = first(intermediate.heat_capacity_volume)
    h = first(intermediate.enthalpy_species)
    u = first(intermediate.internal_energy)
    s = first(intermediate.entropy_species)

    T = only(current.temperature)

    for i in eachindex(mechanism.species) ## save repeated operations
        cₚ[i] = (a₁[i] + a₂[i]T + a₃[i]T^2 + a₄[i]T^3 + a₅[i]T^4)R
        cᵥ[i] = cₚ[i] - R

        h[i] = (a₁[i]T + a₂[i]T^2 / 2 + a₃[i]T^3 / 3 + a₄[i]T^4 / 4 + a₅[i]T^5 / 5 + a₆[i])R
        u[i] = h[i] - (R * T)

        s[i] = (a₁[i] * log(T) + a₂[i]T + a₃[i]T^2 / 2 + a₄[i]T^3 / 3 + a₅[i]T^4 / 4 + a₇[i])R
    end

    return nothing
end

function concentrations((; current, intermediate, mechanism)::Gas{K}; forward=true) where {K<:Number} ##computes species properties

    W⁻¹ = mechanism.inverse_molecular_weight
    α = mechanism.enhancement_factors
    M = first(intermediate.total_molar_concentrations)

    Y = current.mass_fractions
    X = first(current.molar_fractions)
    C = first(current.molar_concentrations)
    ρ = current.density

    W̅ = 1 / sum(Y .* W⁻¹)
    for i in eachindex(mechanism.species)
        YW⁻¹ = Y[i] * W⁻¹[i]
        C[i] = YW⁻¹ * ρ
        X[i] = YW⁻¹ * W̅
    end

    mul!(M, α, C)
    return nothing
end

function reactionconstants((; current, intermediate, mechanism)::Gas{K}; forward=true) where {K<:Number}
    ## optimizable

    T = only(current.temperature)
    RT⁻¹ = inv(R * T)
    RcT⁻¹ = inv(Rc * T)

    isparameters = mechanism.is_reverse_parameters
    isequilibrium = mechanism.is_reversible_equilibrium
    istroe = mechanism.is_troe_parameters
    nr = length(mechanism.reactions)
    nt = length(mechanism.threebody_reactions) ##optimize
    nf = length(mechanism.falloff_reactions)
    ny = nr - nf

    A, b, E = (view(mechanism.pressure_independant_parameters, :, p) for p in 1:3)
    A∞, b∞, E∞ = (view(mechanism.high_pressure_parameters, :, p) for p in 1:3)
    Aₒ, bₒ, Eₒ = (view(mechanism.low_pressure_parameters, :, p) for p in 1:3)
    Aᵣ, bᵣ, Eᵣ = (view(mechanism.reverse_reaction_parameters, :, p) for p in 1:3)

    ν = mechanism.stoichiometric_reactions
    ∑ν = mechanism.stoichiometric_sum
    a, T₃, T₁, T₂ = (view(mechanism.troe_parameters, :, p) for p in 1:4)

    h = first(intermediate.enthalpy_species)
    s = first(intermediate.entropy_species)
    M = first(intermediate.total_molar_concentrations)
    MT = view(M, nt+1:nt+nf)

    H = first(intermediate.enthalpy_reactions)
    S = first(intermediate.entropy_reactions)
    kf = first(intermediate.forward_rate_constant)
    kr = first(intermediate.reverse_rate_constant)

    for k in eachindex(A)
        kf[k] = A[k] * T^b[k] * exp(-E[k] * RcT⁻¹)

    end

    for f in istroe
        kₒ = Aₒ[f] * T^bₒ[f] * exp(-Eₒ[f] * RcT⁻¹)
        k∞ = A∞[f] * T^b∞[f] * exp(-E∞[f] * RcT⁻¹)

        Pᵣ = kₒ * MT[f] / k∞
        log10Pᵣ = log10(Pᵣ)

        Fc = (1.0 - a[f]) * exp(-T / T₃[f]) + a[f] * exp(-T / T₁[f]) + exp(-T₂[f] / T)
        log10Fc = log10(Fc)

        c = -0.4 - 0.67 * log10Fc
        n = 0.75 - 1.27 * log10Fc
        F = 10.0^(log10Fc / (1.0 + ((log10Pᵣ + c) / (n - d * (log10Pᵣ + c)))^2.0))

        kf[ny+f] = F * k∞ * Pᵣ / (1.0 + Pᵣ)
    end

    for f in filter(∉(istroe), eachindex(mechanism.falloff_reactions))
        kₒ = Aₒ[f] * T^bₒ[f] * exp(-Eₒ[f] * RcT⁻¹)
        k∞ = A∞[f] * T^b∞[f] * exp(-E∞[f] * RcT⁻¹)

        Pᵣ = kₒ * MT[f] / k∞
        kf[ny+f] = k∞ * Pᵣ / (1.0 + Pᵣ)
    end

    mul!(H, ν, h)
    mul!(S, ν, s)

    for k in isequilibrium
        Kp = exp((S[k] * T - H[k]) * RT⁻¹)
        Kc = Kp * (Pₐ * RT⁻¹)^∑ν[k]
        kr[k] = kf[k] / Kc
    end

    for k in isparameters
        kr[k] = Aᵣ[k] * T^bᵣ[k] * exp(-Eᵣ[k] * RcT⁻¹)
    end


    return nothing
end

function rates((; current, intermediate, mechanism)::Gas{K}; forward=true) where {K<:Number}

    nt = length(mechanism.threebody_reactions)
    ne = length(mechanism.elementary_reactions)
    threebody_range = ne+1:ne+nt

    W⁻¹ = mechanism.inverse_molecular_weight
    Y = current.mass_fractions
    X = first(current.molar_fractions)
    C = first(current.molar_concentrations)
    W = mechanism.molecular_weight

    nur = mechanism.stoichiometric_reactants
    nup = mechanism.stoichiometric_products
    nut = mechanism.stoichiometric_transpose
    reactants_indicies = mechanism.reactants_indicies
    products_indicies = mechanism.products_indicies

    cᵥ = first(intermediate.heat_capacity_volume)
    u = first(intermediate.internal_energy)
    kf = first(intermediate.forward_rate_constant)
    kr = first(intermediate.reverse_rate_constant)
    M = first(intermediate.total_molar_concentrations)
    MT = view(M, 1:nt)
    q = first(intermediate.rate_of_progress)
    ω̇ = first(intermediate.production_rate)
    Ẏ = first(intermediate.mass_change_rate)
    Ṫ = first(intermediate.temperature_change_rate)

    ρ = current.density

    for k in eachindex(mechanism.reactions)

        stepf = 1.0
        for i in reactants_indicies[k]
            stepf *= C[i]^nur[k, i]
        end

        stepr = 1.0
        for i in products_indicies[k]
            stepr *= C[i]^nup[k, i]
        end

        k ∈ threebody_range ? q[k] = MT[k-ne] * (kf[k] * stepf - kr[k] * stepr) : q[k] = kf[k] * stepf - kr[k] * stepr

    end

    mul!(ω̇, nut, q)

    W̅ = 1 / sum(Y .* W⁻¹)
    c̅ᵥ = sum(cᵥ .* X) / W̅
    #c̅ᵥ = sum(cᵥ .* W⁻¹ .* Y)
    #println(c̅ᵥ)
    for i in eachindex(Ẏ)
        Ẏ[i] = ω̇[i] * W[i] / ρ
    end

    Ṫ[1] = -sum(ω̇ .* u) / (ρ * c̅ᵥ) ## worth it?

    return nothing
end

function step!(gas::Gas{K}, Y::AbstractVector{K}, T::K; forward=true) where {K<:Number} ### AbstractVector allocates +1

    #gas.current = gas.initial
    gas.current.mass_fractions = Y
    gas.current.temperature[1] = T
    gas.current.density = gas.initial.density ## in init?

    interval(gas; forward=true)
    polynomials(gas; forward=true)
    concentrations(gas; forward=true)
    reactionconstants(gas; forward=true)
    rates(gas; forward=true)

    return nothing
end