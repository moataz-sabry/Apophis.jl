function arrhenius(A::K, b::K, E::K, T::K) where {K<:Float64}
    RcT⁻¹ = inv(Rc * T)
    k = A * T^b * exp(-E * RcT⁻¹)
    dkdA = T^b * exp(-E * RcT⁻¹)
    dkdT = (b + E * RcT⁻¹) * k / T
    return k, dkdT, dkdA
end

function reducedpressure(kₒ::K, M::K, k∞::K) where {K<:Float64}
    k∞⁻¹ = inv(k∞)
    Pᵣ = kₒ * M * k∞⁻¹
    dPᵣdM = kₒ * k∞⁻¹
    dPᵣdkₒ = M * k∞⁻¹
    dPᵣdk∞ = -Pᵣ * k∞⁻¹
    return Pᵣ, dPᵣdM, dPᵣdkₒ, dPᵣdk∞
end

function troerate(k∞::K, Pᵣ::K, F::K) where {K<:Float64}
    invPᵣ = inv(1.0 + Pᵣ)
    kf = F * k∞ * Pᵣ * invPᵣ
    dkfdF = k∞ * Pᵣ / (1.0 + Pᵣ)
    dkfdk∞ = F * Pᵣ * invPᵣ
    dkfdPᵣ = (F * k∞ - kf) * invPᵣ
    return kf, dkfdF, dkfdk∞, dkfdPᵣ
end

function lindemannrate(k∞::K, Pᵣ::K) where {K<:Float64}
    invPᵣ = inv(1.0 + Pᵣ)
    kf = k∞ * Pᵣ * invPᵣ
    dkfdk∞ = Pᵣ * invPᵣ
    dkfdPᵣ = (k∞ - kf) * invPᵣ
    return kf, dkfdk∞, dkfdPᵣ
end

function fcent(a::K, T₃::K, T₁::K, T₂::K, T::K) where {K<:Float64}
    v₁ = exp(-T / T₃)
    v₂ = a * exp(-T / T₁)
    v₃ = exp(-T₂ / T)
    v₄ = (a - 1.0)

    Fc = -v₄ * v₁ + v₂ + v₃
    dFcdT = v₄ * v₁ / T₃ - v₂ / T₁ + v₃ * T₂ / T^2
    return Fc, dFcdT
end

function interval((; current, intermediate, mechanism)::Gas{K}) where {K<:Float64}

    Tₘ = mechanism.common_temperature
    A = mechanism.upper_temperature_coefficients
    a = mechanism.lower_temperature_coefficients

    coeffs = intermediate.polynomial_coefficients
    T = only(current.temperature)

    for i in eachindex(mechanism.species) # 500 ns; faster than coeffs[:, i] .= view(A, :, i); 690 ns for GRI3; eachcol?
        offset = 7i - 6  # 7(i - 1) + 1
        T ≥ Tₘ[i] ? copyto!(coeffs, offset, A, offset, 7) : copyto!(coeffs, offset, a, offset, 7)
    end
    return nothing
end

function polynomials((; current, intermediate, mechanism)::Gas{K}) where {K<:Float64} # computes nasa polynomials for species entropy, entahlpy, ...

    coeffs = intermediate.polynomial_coefficients
    a₁, a₂, a₃, a₄, a₅, a₆, a₇ = (view(coeffs, a, :) for a in 1:7) # multiple assign?

    cₚ, dcₚdT, _ = intermediate.heat_capacity_pressure
    cᵥ, dcᵥdT, _ = intermediate.heat_capacity_volume
    h, dhdT, _ = intermediate.enthalpy_species
    u, dudT, _ = intermediate.internal_energy
    s, dsdT, _ = intermediate.entropy_species

    T = only(current.temperature)
    T², T³, T⁴, T⁵ = T^2, T^3, T^4, T^5 ## necessary?

    for i in eachindex(mechanism.species) # save repeated operations
        cₚ[i] = (a₁[i] + a₂[i]T + a₃[i]T² + a₄[i]T³ + a₅[i]T⁴)R
        dcₚdT[i] = (a₂[i] + 2a₃[i]T + 3a₄[i]T² + 4a₅[i]T³)R

        cᵥ[i] = cₚ[i] - R
        dcᵥdT[i] = dcₚdT[i]

        h[i] = (a₁[i]T + a₂[i]T² / 2 + a₃[i]T³ / 3 + a₄[i]T⁴ / 4 + a₅[i]T⁵ / 5 + a₆[i])R
        dhdT[i] = (a₁[i] + a₂[i]T + a₃[i]T² + a₄[i]T³ + a₅[i]T⁴)R

        u[i] = h[i] - R * T
        dudT[i] = dhdT[i] - R

        s[i] = (a₁[i] * log(T) + a₂[i]T + a₃[i]T² / 2 + a₄[i]T³ / 3 + a₅[i]T⁴ / 4 + a₇[i])R
        dsdT[i] = (a₁[i] / T + a₂[i] + a₃[i]T + a₄[i]T² + a₅[i]T³)R
    end

    return nothing
end

function concentrations((; current, intermediate, mechanism)::Gas{K}) where {K<:Float64} # computes species properties

    W⁻¹ = mechanism.inverse_molecular_weight
    α = mechanism.enhancement_factors
    M, _, _ = intermediate.total_molar_concentrations

    Y = current.mass_fractions
    C, _, dCdY = current.molar_concentrations
    ρ = current.density

    for i in eachindex(mechanism.species)
        C[i] = Y[i] * W⁻¹[i] * ρ
        dCdY[i, i] = W⁻¹[i] * ρ
    end

    mul!(M, α, C) # dMdC = α

    return nothing
end

function forwardrates((; current, intermediate, mechanism)::Gas{K}) where {K<:Float64}
    # optimizable

    T = only(current.temperature)

    istroe = mechanism.is_troe_parameters
    nt = length(mechanism.threebody_reactions) #optimize
    ne = length(mechanism.elementary_reactions)
    offset = ne + nt

    α = mechanism.enhancement_factors
    a, T₃, T₁, T₂ = (view(mechanism.troe_parameters, :, p) for p in 1:4)
    A∞, b∞, E∞ = (view(mechanism.high_pressure_parameters, :, p) for p in 1:3)
    Aₒ, bₒ, Eₒ = (view(mechanism.low_pressure_parameters, :, p) for p in 1:3)
    A, b, E = (view(mechanism.pressure_independant_parameters, :, p) for p in 1:3)

    M, _, _ = intermediate.total_molar_concentrations
    kf, dkfdT, dkfdY, dkfdA = intermediate.forward_rate_constant
    _, _, dCdY = current.molar_concentrations

    for k in axes(mechanism.pressure_independant_parameters, 1)
        kf[k], dkfdT[k], dkfdA[k, k] = arrhenius(A[k], b[k], E[k], T)
    end

    for f in istroe
        kₒ, dkₒdT, dkₒdAₒ = arrhenius(Aₒ[f], bₒ[f], Eₒ[f], T)
        k∞, dk∞dT, dk∞dA∞ = arrhenius(A∞[f], b∞[f], E∞[f], T)

        Pᵣ, dPᵣdM, dPᵣdkₒ, dPᵣdk∞ = reducedpressure(kₒ, M[nt+f], k∞)

        log10Pᵣ = log10(Pᵣ)
        dlgPᵣdPᵣ = inv(log(10) * Pᵣ)

        Fc, dFcdT = fcent(a[f], T₃[f], T₁[f], T₂[f], T)

        log10Fc = log10(Fc)
        dlgFcdFc = inv(log(10) * Fc)

        c = -0.4 - 0.67 * log10Fc
        dcdFc = -0.67 * dlgFcdFc

        n = 0.75 - 1.27 * log10Fc
        dndFc = -1.27 * dlgFcdFc

        ### check 
        v₁ = c + log10Pᵣ
        v₂ = n - d * v₁
        v₃ = 1.0 + v₁^2 / v₂^2
        v₄ = log10Fc / v₃
        v₅ = 10.0^(v₄) * log(10.0)
        v₆ = 2 * log10Fc * v₅ * v₁
        v₇ = v₆ * v₁
        v₈ = v₂^2 * v₃^2
        v₉ = v₂ * v₈
        v₀ = -(v₆ / v₈ + v₇ * d / v₉)

        F = 10.0^v₄

        dFdFc = v₅ / v₃ * dlgFcdFc
        dFdPᵣ = v₀ * dlgPᵣdPᵣ
        dFdc = v₀
        dFdn = v₇ / v₉
        ##check
        dFdT = (dFdFc + dFdc * dcdFc + dFdn * dndFc) * dFcdT + dFdPᵣ * (dPᵣdk∞ * dk∞dT + dPᵣdkₒ * dkₒdT)
        dFdA = dFdPᵣ * (dPᵣdk∞ * dk∞dA∞)# + dPᵣdkₒ * dkₒdAₒ)

        kf[offset+f], dkfdF, dkfdk∞, dkfdPᵣ = troerate(k∞, Pᵣ, F)

        dkfdT[offset+f] = dkfdF * dFdT + dkfdk∞ * dk∞dT + dkfdPᵣ * (dPᵣdk∞ * dk∞dT + dPᵣdkₒ * dkₒdT)
        dkfdA[offset+f, offset+f] = dkfdF * dFdA + dkfdk∞ * dk∞dA∞ + dkfdPᵣ * dPᵣdk∞ * dk∞dA∞# + dPᵣdkₒ * dkₒdAₒ)

        for i in eachindex(mechanism.species)
            dkfdY[offset+f, i] = (dkfdF * dFdPᵣ + dkfdPᵣ) * dPᵣdM * α[nt+f, i] * dCdY[i, i]
        end
    end

    for f in filter(∉(istroe), eachindex(mechanism.falloff_reactions))
        kₒ, dkₒdT, dkₒdAₒ = arrhenius(Aₒ[f], bₒ[f], Eₒ[f], T)
        k∞, dk∞dT, dk∞dA∞ = arrhenius(A∞[f], b∞[f], E∞[f], T)

        Pᵣ, dPᵣdM, dPᵣdkₒ, dPᵣdk∞ = reducedpressure(kₒ, M[nt+f], k∞)
        kf[offset+f], dkfdk∞, dkfdPᵣ = lindemannrate(k∞, Pᵣ)

        dkfdT[offset+f] = dkfdk∞ * dk∞dT + dkfdPᵣ * (dPᵣdk∞ * dk∞dT + dPᵣdkₒ * dkₒdT)
        dkfdA[offset+f, offset+f] = dkfdk∞ * dk∞dA∞ + dkfdPᵣ * dPᵣdk∞ * dk∞dA∞# + dPᵣdkₒ * dkₒdAₒ)

        for i in eachindex(mechanism.species)
            dkfdY[offset+f, i] = dkfdPᵣ * dPᵣdM * α[nt+f, i] * dCdY[i, i]
        end
    end

    return nothing

end

function reverserates((; current, intermediate, mechanism)::Gas{K}) where {K<:Float64}
    # optimizable

    T = only(current.temperature)
    RT⁻¹ = inv(R * T)

    isparameters = mechanism.is_reverse_parameters
    isequilibrium = mechanism.is_reversible_equilibrium
    ν = mechanism.stoichiometric_reactions
    ∑ν = mechanism.stoichiometric_sum

    Aᵣ, bᵣ, Eᵣ = (view(mechanism.reverse_reaction_parameters, :, p) for p in 1:3)

    h, dhdT, _ = intermediate.enthalpy_species
    s, dsdT, _ = intermediate.entropy_species
    H, dHdT, _ = intermediate.enthalpy_reactions
    S, dSdT, _ = intermediate.entropy_reactions
    kf, dkfdT, dkfdY, dkfdA = intermediate.forward_rate_constant
    kr, dkrdT, dkrdY, dkrdA = intermediate.reverse_rate_constant

    mul!(H, ν, h)
    mul!(S, ν, s)
    dHdh = dSds = ν

    mul!(dHdT, dHdh, dhdT) # had to
    mul!(dSdT, dSds, dsdT)

    for k in isequilibrium # isequilibrium

        Kp = exp((S[k] * T - H[k]) * RT⁻¹)
        dKpdS = Kp / R
        dKpdH = -Kp * RT⁻¹
        dKpdT = S[k] * Kp * RT⁻¹ - Kp * (S[k] * T - H[k]) / (R * T^2) ## Matrix Calculus

        Kc = Kp * (Pₐ * RT⁻¹)^∑ν[k]
        dKcdKp = (Pₐ * RT⁻¹)^∑ν[k]
        dKcdT = -Kp * ∑ν[k] * T^-(1.0 + ∑ν[k]) * R^-∑ν[k] * Pₐ^∑ν[k] ##

        kr[k] = kf[k] / Kc
        dkrdkf = inv(Kc)
        dkrdKc = -kf[k] * inv(Kc)^2

        dkrdT[k] = dkrdkf * dkfdT[k] + dkrdKc * (dKcdT + dKcdKp * (dKpdT + dKpdH * dHdT[k] + dKpdS * dSdT[k]))
        dkrdA[k, k] = dkrdkf * dkfdA[k, k]
        for i in eachindex(mechanism.species) ## only falloff 
            dkrdY[k, i] = dkrdkf * dkfdY[k, i]
        end
    end

    for k in isparameters
        kr[k], dkrdT[k] = arrhenius(Aᵣ[k], bᵣ[k], Eᵣ[k], T)
    end

    return nothing
end

function productionrates((; current, intermediate, mechanism)::Gas{K}) where {K<:Float64}

    nt = length(mechanism.threebody_reactions)
    ne = length(mechanism.elementary_reactions)
    threebody_range = ne+1:ne+nt

    α = mechanism.enhancement_factors
    nut = mechanism.stoichiometric_transpose
    nup = mechanism.stoichiometric_products
    nur = mechanism.stoichiometric_reactants
    reactants_indicies = mechanism.reactants_indicies
    products_indicies = mechanism.products_indicies

    kf, dkfdT, dkfdY, dkfdA = intermediate.forward_rate_constant
    kr, dkrdT, dkrdY, dkrdA = intermediate.reverse_rate_constant
    M, _, _ = intermediate.total_molar_concentrations
    MT = view(M, 1:nt)
    q, dqdT, dqdY, dqdA = intermediate.rate_of_progress
    ω̇, _, _ = intermediate.production_rate

    C, _, dCdY = current.molar_concentrations


    for i in eachindex(mechanism.reactions)

        stepf = 1.0
        for k in reactants_indicies[i]
            stepf *= C[k]^nur[i, k]
        end

        stepr = 1.0
        for k in products_indicies[i]
            stepr *= C[k]^nup[i, k]
        end

        q[i] = kf[i] * stepf - kr[i] * stepr
        dqdT[i] = dkfdT[i] * stepf - dkrdT[i] * stepr
        dqdA[i, i] = dkfdA[i, i] * stepf - dkrdA[i, i] * stepr

        for j in eachindex(mechanism.species) ## only falloff
            dqdY[i, j] = stepf * dkfdY[i, j] - stepr * dkrdY[i, j]
        end

        for j in reactants_indicies[i]
            prod = 1.0
            for k in reactants_indicies[i]
                isequal(k, j) ? prod *= nur[i, k] * C[k]^(nur[i, k] - 1) : prod *= C[k]^nur[i, k]
            end
            dqdY[i, j] += prod * kf[i] * dCdY[j, j]
        end

        for j in products_indicies[i]
            prod = 1.0
            for k in products_indicies[i]
                isequal(k, j) ? prod *= nup[i, k] * C[k]^(nup[i, k] - 1) : prod *= C[k]^nup[i, k]
            end
            dqdY[i, j] -= prod * kr[i] * dCdY[j, j]
        end
    end

    for i in filter(∈(threebody_range), eachindex(mechanism.reactions))
        for j in eachindex(mechanism.species)
            dqdY[i, j] *= MT[i-ne]
            dqdY[i, j] += q[i] * α[i-ne, j] * dCdY[j, j]
        end
        q[i] *= MT[i-ne]
        dqdT[i] *= MT[i-ne]
        dqdA[i, i] *= MT[i-ne]
    end

    mul!(ω̇, nut, q)
    return nothing
end

function rhs((; current, intermediate, mechanism)::Gas{K}) where {K<:Float64}

    W⁻¹ = mechanism.inverse_molecular_weight
    W = mechanism.molecular_weight
    nut = mechanism.stoichiometric_transpose

    cᵥ, dcᵥdT, _ = intermediate.heat_capacity_volume
    u, dudT, _ = intermediate.internal_energy
    _, dqdT, dqdY, dqdA = intermediate.rate_of_progress
    ω̇, dω̇dT, dω̇dY, dω̇dA = intermediate.production_rate
    Ẏ, dẎdT, dẎdY, dẎdA = intermediate.mass_change_rate
    Ṫ, dṪdT, dṪdY, dṪdA = intermediate.temperature_change_rate

    ρ = current.density
    Y = current.mass_fractions
    X, _, dXdY = current.molar_fractions

    dω̇dq = nut
    mul!(dω̇dT, dω̇dq, dqdT)
    mul!(dω̇dY, dω̇dq, dqdY)
    mul!(dω̇dA, dω̇dq, dqdA)
    #show(stdout, "text/plain", dω̇dA)

    W̅ = inv(Y ⋅ W⁻¹)
    c̅ᵥ = 0.0
    dc̅ᵥdT = 0.0
    for i in eachindex(mechanism.species)
        c̅ᵥ += cᵥ[i] * W⁻¹[i] * Y[i]
        dc̅ᵥdT += dcᵥdT[i] * W⁻¹[i] * Y[i]
    end

    Ṫ[1] = -(ω̇ ⋅ u) / (ρ * c̅ᵥ) # worth it? #dTdt
    dṪdT[1] = -(Ṫ[1] * dc̅ᵥdT / c̅ᵥ) - (dω̇dT ⋅ u + ω̇ ⋅ dudT) / (ρ * c̅ᵥ)

    mul!(dẎdA, Diagonal(W), dω̇dA, inv(ρ), 0.0) # checked √
    mul!(dṪdA, u', dω̇dA, -inv((ρ * c̅ᵥ)), 0.0) # checked √

    for i in eachindex(mechanism.species)

        X[i] = Y[i] * W⁻¹[i] * W̅
        Ẏ[i] = ω̇[i] * W[i] / ρ #dYdt
        dẎdT[i] = dω̇dT[i] * W[i] / ρ #
        dṪdY[1, i] = -Ṫ[1] * cᵥ[i] * W⁻¹[i] / c̅ᵥ #

        for j in eachindex(mechanism.species)
            dXdY[i, j] = -Y[i] * W⁻¹[i] * W̅^2 * W⁻¹[j]
            dẎdY[i, j] = dω̇dY[i, j] * W[i] / ρ #
            dṪdY[i] -= u[j] * dω̇dY[j, i] / (ρ * c̅ᵥ) #
        end
        dXdY[i, i] += W̅ * W⁻¹[i]
    end
    #show(stdout, "text/plain", dẎdA)
    #show(stdout, "text/plain", dṪdA)
    return nothing
end

function step!(gas::Gas{K}, Y::AbstractVector{K}, T::K) where {K<:Float64} ## AbstractVector allocates +1

    #gas.current = gas.initial
    gas.current.mass_fractions = Y
    gas.current.temperature[1] = T

    interval(gas)
    polynomials(gas)
    concentrations(gas)
    forwardrates(gas)
    reverserates(gas)
    productionrates(gas)
    rhs(gas)

    return nothing
end