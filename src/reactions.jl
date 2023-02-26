struct Arrhenius{N<:Number} ## [A] := M^(∑νᵣ - 1) ⋅ s^-1; where M = cm^3 ⋅ mol^-1, if not otherwise specified in kinetics data file
    A::N
    β::N
    E::N
end

function Arrhenius(itr; order=1, mechunits=chemkin_default_units)
    activation_energy_unit, length_unit, quantity_unit = (mechunits[k] for k in ("activation-energy", "length", "quantity"))
    isempty(itr) && return nothing
    A, β, E = itr
    return Arrhenius(unitfy_rate(A, order, length_unit, quantity_unit), β, unitfy_activation_energy(E, activation_energy_unit))
end

((; A, β, E)::Arrhenius{N})(T::N) where {N<:Number} = A * T^β * exp(-E * inv(Rc * T))
(arrhenius::Arrhenius{N})(::Val{:dg}, T::N) where {N<:Number} = gradient(arrhenius -> arrhenius(T), arrhenius) |> only

function ((; A, β, E)::Arrhenius{N})(::Val{:dT}, T::N) where {N<:Number}
    t = exp(-E * inv(Rc * T))
    kf = A * T^β * t
    dkfdT = kf * (inv(T)β + inv(Rc * T^2)E)
    return kf, dkfdT
end

reduced_pressure(kₒ::N, M::N, k∞::N) where {N<:Number} = kₒ * M * inv(k∞)
reduced_pressure(::Val{:dkₒ}, M::N, k∞::N) where {N<:Number} = M * inv(k∞)
reduced_pressure(::Val{:dM}, kₒ::N, k∞::N) where {N<:Number} = kₒ * inv(k∞)
reduced_pressure(::Val{:dk∞}, kₒ::N, M::N, k∞::N) where {N<:Number} = -kₒ * M * inv(k∞)^2

struct Plog{N<:Number}
    pressures::Vector{N}
    arrhenius::Vector{Arrhenius{N}}
end

Plog(itr...; order=1) = isempty(first(itr)) ? nothing : Plog(first(itr), Arrhenius{Float64}[Arrhenius(p; order) for p in zip(rest(itr, 2)...)])

function isin(p::N, x::Vector{N}) where {N<:Number}
    for i in eachindex(x)
        x[i] ≤ p && x[i+1] > p && return i
    end
end

isout(p::N, x::Vector{N}) where {N<:Number} = p > last(x) ? lastindex(x) : p < first(x) ? firstindex(x) : nothing ## on the assumption that plogs are sorted, To-Do: sort in advance
interpolate(kᵢ::N, kᵢ₊₁::N, Pᵢ::N, Pᵢ₊₁::N, P::N) where {N<:Number} = log(kᵢ) + (log(kᵢ₊₁) - log(kᵢ)) * (log(P) - log(Pᵢ)) / (log(Pᵢ₊₁) - log(Pᵢ))
interpolate(::Val{:dP}, kᵢ::N, kᵢ₊₁::N, Pᵢ::N, Pᵢ₊₁::N, P::N) where {N<:Number} = (log(kᵢ₊₁) - log(kᵢ)) / (log(Pᵢ₊₁) - log(Pᵢ))P

function ((; pressures, arrhenius)::Plog{N})(T::N, P::N) where {N<:Number}
    i = isout(P, pressures)
    isnothing(i) || return arrhenius[i](T)

    i = isin(P, pressures)
    Pᵢ, Pᵢ₊₁ = P[i], P[i+1]
    kᵢ, kᵢ₊₁ = arrhenius[i](T), arrhenius[i+1](T)
    return interpolate(kᵢ, kᵢ₊₁, Pᵢ, Pᵢ₊₁, P)
end

function ((; pressures, arrhenius)::Plog{N})(v::Val{:dT}, T::N, P::N) where {N<:Number}
    i = isout(P, pressures)
    isnothing(i) || return arrhenius[i](v, T)

    i = isin(P, pressures)
    Pᵢ, Pᵢ₊₁ = P[i], P[i+1]
    kᵢ, dkᵢdT, kᵢ₊₁, dkᵢ₊₁dT = arrhenius[i](v, T), arrhenius[i+1](v, T)
    return interpolate(kᵢ, kᵢ₊₁, Pᵢ, Pᵢ₊₁, P), interpolate(dkᵢdT / kᵢ, dkᵢ₊₁dT / kᵢ₊₁, Pᵢ, Pᵢ₊₁, P)
end

function ((; pressures, arrhenius)::Plog{N})(v::Val{:dP}, T::N, P::N) where {N<:Number}
    i = isout(P, pressures)
    isnothing(i) || return arrhenius[i](T), zero(N)

    i = isin(P, pressures)
    Pᵢ, Pᵢ₊₁ = P[i], P[i+1]
    kᵢ, kᵢ₊₁ = arrhenius[i](T), arrhenius[i+1](T)
    return interpolate(kᵢ, kᵢ₊₁, Pᵢ, Pᵢ₊₁, P), interpolate(v, kᵢ, kᵢ₊₁, Pᵢ, Pᵢ₊₁, P)
end

struct Troe{N<:Number}
    a::N
    T₃::N
    T₁::N
    T₂::N
end

Troe() = nothing
Troe(itr) = Troe(itr...)
Troe(a::N, T₃::N, T₁::N) where {N<:Number} = Troe(a, T₃, T₁, zero(N))
Troe(a::N, T₃::N, T₁::N, ::Nothing) where {N<:Number} = Troe(a, T₃, T₁, zero(N))

((; a, T₃, T₁, T₂)::Troe{N})(T::N) where {N<:Number} = (1.0 - a) * exp(-T / T₃) + a * exp(-T / T₁) + (iszero(T₂) ? zero(N) : exp(-T₂ / T))
((; a, T₃, T₁, T₂)::Troe{N})(::Val{:dT}, T::N) where {N<:Number} = (a - one(N)) * exp(-T / T₃) / T₃ - a * exp(-T / T₁) / T₁ + (iszero(T₂) ? zero(N) : exp(-T₂ / T) * T₂ / T^2)
(troe::Troe{N})(::Val{:dg}, T::N) where {N<:Number} = gradient(troe -> troe(T), troe) |> only

function troe_function(Fc::N, Pᵣ::N) where {N<:Number}
    log10Fc = log10(Fc)
    c = -0.4 - 0.67log10Fc
    n = 0.75 - 1.27log10Fc

    t = log10(Pᵣ) + c
    F = 10.0^(log10Fc * inv(one(N) + (inv(n - d * t)t)^2))
    return F
end

function troe_function(::Val{:dT}, Fc::N, Pᵣ::N) where {N<:Number} ## Zygote is somehow slow here
    log10Pᵣ = log10(Pᵣ)
    dlog10PᵣdPᵣ = log(10.0)Pᵣ |> inv
    log10Fc = log10(Fc)
    dlog10FcdFc = log(10.0)Fc |> inv

    c = -0.4 - 0.67log10Fc
    dcdlog10Fc = -0.67
    n = 0.75 - 1.27log10Fc
    dndlog10Fc = -1.27

    t₀ = c + log10Pᵣ
    t₁ = t₀^2
    t₂ = n - (d * t₀)
    t₃ = t₂^2
    t₄ = one(N) + (t₁ / t₃)

    F = 10.0^(log10Fc / t₄)
    t₅ = log(10.0)F
    t₆ = t₄^2
    t₇ = 2.0log10Fc * t₁ * t₅

    dFdlog10Fc = t₅ / t₄
    dFdlog10Pᵣ = dFdc = -(2.0log10Fc * t₀ * t₅ * inv(t₆ * t₃) + d * t₇ * inv(t₂^3 * t₆))
    dFdn = t₇ * inv(t₂^3 * t₄^2)

    dFdFc = (dFdlog10Fc + dFdc * dcdlog10Fc + dFdn * dndlog10Fc) * dlog10FcdFc
    dFdPᵣ = dFdlog10Pᵣ * dlog10PᵣdPᵣ

    return F, dFdFc, dFdPᵣ
end

function troe_function(::Val{:dC}, Fc::N, Pᵣ::N) where {N<:Number}
    log10Pᵣ = log10(Pᵣ)
    dlog10PᵣdPᵣ = log(10.0)Pᵣ |> inv
    log10Fc = log10(Fc)

    c = -0.4 - 0.67log10Fc
    n = 0.75 - 1.27log10Fc

    t₀ = c + log10Pᵣ
    t₁ = t₀^2
    t₂ = n - (d * t₀)
    t₃ = t₂^2
    t₄ = one(N) + (t₁ / t₃)

    F = 10.0^(log10Fc / t₄)
    t₅ = log(10.0)F
    t₆ = t₄^2
    t₇ = 2.0log10Fc * t₁ * t₅

    dFdlog10Pᵣ = -(2.0log10Fc * t₀ * t₅ * inv(t₆ * t₃) + d * t₇ * inv(t₂^3 * t₆))
    dFdPᵣ = dFdlog10Pᵣ * dlog10PᵣdPᵣ
    return F, dFdPᵣ
end

struct ElementaryReaction{N<:Number, S<:AbstractSpecies{N}} <: AbstractReaction{N}
    i::Int
    equation::Symbol
    isreversible::Bool
    reactants::Vector{Pair{S, N}}
    products::Vector{Pair{S, N}}
    reaction_order::N
    forward_rate_parameters::Maybe{Arrhenius{N}}
    reverse_rate_parameters::Maybe{Arrhenius{N}}
    plog_parameters::Maybe{Plog{N}}
    rates::ReactionRates{N}
end

struct ThreeBodyReaction{N<:Number, S<:AbstractSpecies{N}} <: AbstractReaction{N}
    i::Int
    equation::Symbol
    isreversible::Bool
    reactants::Vector{Pair{S, N}}
    products::Vector{Pair{S, N}}
    reaction_order::N
    forward_rate_parameters::Arrhenius{N}
    reverse_rate_parameters::Maybe{Arrhenius{N}}
    enhancement_factors::Vector{Pair{S, N}}
    rates::ReactionRates{N}
end

struct FallOffReaction{N<:Number, S<:AbstractSpecies{N}} <: AbstractReaction{N}
    i::Int
    equation::Symbol
    isreversible::Bool
    reactants::Vector{Pair{S, N}}
    products::Vector{Pair{S, N}}
    reaction_order::N
    high_pressure_parameters::Arrhenius{N}
    low_pressure_parameters::Arrhenius{N}
    troe_parameters::Maybe{Troe{N}}
    reverse_rate_parameters::Maybe{Arrhenius{N}}
    enhancement_factors::Vector{Pair{S, N}}
    rates::ReactionRates{N}
end

const Reaction{N<:Number, S<:AbstractSpecies{N}} = Union{ElementaryReaction{N, S}, ThreeBodyReaction{N, S}, FallOffReaction{N, S}}

total_molar_concentration(C::Vector{N}, α::Vector{Pair{Species{N},N}}) where {N<:Number} = @inbounds sum(C[k] * f for ((; k), f) in α)

forward_rate((; rates, forward_rate_parameters, plog_parameters)::ElementaryReaction{N}, (; T, P)::State{N}) where {N<:Number} = setindex!(rates.kf.val, isnothing(plog_parameters) ? forward_rate_parameters(T) : plog_parameters(T, P), 1)
forward_rate((; rates, forward_rate_parameters)::ThreeBodyReaction{N}, (; T)::State{N}) where {N<:Number} = setindex!(rates.kf.val, forward_rate_parameters(T), 1)

function forward_rate(reaction::FallOffReaction{N}, (; T, C)::State{N}) where {N<:Number}
    (; rates, high_pressure_parameters, low_pressure_parameters, troe_parameters, enhancement_factors) = reaction

    k∞, kₒ = (high_pressure_parameters, low_pressure_parameters)(T)
    M = total_molar_concentration(C, enhancement_factors)
    iszero(M) && (@inbounds rates.kf.val[] = zero(N); return nothing)

    Pᵣ = reduced_pressure(kₒ, M, k∞)
    F = isnothing(troe_parameters) ? one(N) : troe_parameters(T) |> Fix2(troe_function, Pᵣ)
    @inbounds rates.kf.val[] = F * k∞ * Pᵣ * inv(one(N) + Pᵣ)
    return nothing
end

function forward_rate(v::Val{:dT}, (; rates, forward_rate_parameters, plog_parameters)::ElementaryReaction{N}, (; T, P)::State{N}) where {N<:Number}
    @inbounds rates.kf.val[], rates.kf.dT[] = isnothing(plog_parameters) ? forward_rate_parameters(v, T) : plog_parameters(v, T, P)
    return nothing
end

function forward_rate(v::Val{:dT}, (; rates, forward_rate_parameters)::ThreeBodyReaction{N}, (; T)::State{N}) where {N<:Number}
    @inbounds rates.kf.val[], rates.kf.dT[] = forward_rate_parameters(v, T)
    return nothing
end

function forward_rate(v::Val{:dT}, reaction::FallOffReaction{N}, (; T, C)::State{N}) where {N<:Number}
    (; rates, high_pressure_parameters, low_pressure_parameters, troe_parameters, enhancement_factors) = reaction

    (k∞, dk∞dT), (kₒ, dkₒdT) = (high_pressure_parameters, low_pressure_parameters)(v, T)
    M = total_molar_concentration(C, enhancement_factors)
    iszero(M) && (@inbounds rates.kf.val[] = rates.kf.dT[] = zero(N); return nothing)

    Pᵣ = reduced_pressure(kₒ, M, k∞)
    t = inv(one(N) + Pᵣ)

    dPᵣdkₒ = reduced_pressure(Val(:dkₒ), M, k∞)
    dPᵣdk∞ = reduced_pressure(Val(:dk∞), kₒ, M, k∞)

    dPᵣdT = dPᵣdkₒ * dkₒdT + dPᵣdk∞ * dk∞dT

    if isnothing(troe_parameters)
        @inbounds rates.kf.val[] = k∞ * Pᵣ * t
        @inbounds rates.kf.dT[] = t * (dk∞dT * Pᵣ + dPᵣdT * k∞ * t)
    else
        Fc = troe_parameters(T)
        dFcdT = troe_parameters(v, T)

        F, dFdFc, dFdPᵣ = troe_function(v, Fc, Pᵣ)
        dFdT = dFdFc * dFcdT + dFdPᵣ * dPᵣdT

        @inbounds rates.kf.val[] = F * k∞ * Pᵣ * t
        @inbounds rates.kf.dT[] = t * (F * (dk∞dT * Pᵣ + dPᵣdT * k∞ * t) + dFdT * k∞ * Pᵣ)
    end
    return nothing
end

forward_rate(::Val{:dC}, reaction::ElementaryReaction{N}, state::State{N}) where {N<:Number} = forward_rate(reaction, state)
forward_rate(::Val{:dC}, reaction::ThreeBodyReaction{N}, state::State{N}) where {N<:Number} = forward_rate(reaction, state)

function forward_rate(v::Val{:dC}, reaction::FallOffReaction{N}, (; T, C)::State{N}) where {N<:Number}
    (; rates, high_pressure_parameters, low_pressure_parameters, troe_parameters, enhancement_factors) = reaction

    k∞, kₒ = (high_pressure_parameters, low_pressure_parameters)(T)
    M = total_molar_concentration(C, enhancement_factors)
    iszero(M) && (@inbounds rates.kf.dC .= rates.kf.val[] = zero(N); return nothing)

    Pᵣ = reduced_pressure(kₒ, M, k∞)
    t = inv(one(N) + Pᵣ)

    dPᵣdM = reduced_pressure(Val(:dM), kₒ, k∞)
    F, dFdPᵣ = isnothing(troe_parameters) ? (one(N), zero(N)) : troe_parameters(T) |> Fc -> troe_function(v, Fc, Pᵣ)

    @inbounds rates.kf.val[] = F * k∞ * Pᵣ * t
    dkfdPᵣ = F * k∞ * t^2
    dkfdF = k∞ * Pᵣ * t
    for ((; k), dMdCₖ) in enhancement_factors
        @inbounds rates.kf.dC[k] = (dkfdPᵣ + dkfdF * dFdPᵣ) * dPᵣdM * dMdCₖ
    end
    return nothing
end

function forward_rate(v::Val{:dP}, (; rates, forward_rate_parameters, plog_parameters)::ElementaryReaction{N}, (; T, P)::State{N}) where {N<:Number}
    rates.kf.val[], rates.kf.dT[] = isnothing(plog_parameters) ? (forward_rate_parameters(T), zero(N)) : plog_parameters(v, T, P)
    return nothing
end

change_enthalpy((; reactants, products)::AbstractReaction{<:Number}, i::Int=1) =
    @inbounds sum(getfield(s.thermo.h, i)[] * ν for (s, ν) in flatten((reactants, products)))

change_entropy((; reactants, products)::AbstractReaction{<:Number}, i::Int=1) =
    @inbounds sum(getfield(s.thermo.s, i)[] * ν for (s, ν) in flatten((reactants, products)))

function equilibrium_constants(reaction::AbstractReaction{N}, T::N) where {N<:Number}
    ∑ν = reaction.reaction_order
    t = inv(R * T)
    ∆H = change_enthalpy(reaction)
    ∆S = change_entropy(reaction)

    Kp = exp((∆S * T - ∆H) * t)
    Kc = Kp * (Pa * t)^∑ν
    return Kc
end

function equilibrium_constants(::Val{:dT}, reaction::AbstractReaction{N}, T::N) where {N<:Number}
    ∆H, d∆HdT = (change_enthalpy(reaction, f) for f in OneTo(2))
    ∆S, d∆SdT = (change_entropy(reaction, f) for f in OneTo(2))

    ∑ν = reaction.reaction_order
    t₀ = inv(R * T)
    t₁ = ∆S * T - ∆H

    Kp = exp(t₁ * t₀)
    ∂Kp∂∆S = Kp / R
    ∂Kp∂∆H = -Kp * t₀
    ∂Kp∂T = (∆S - t₁ / T) * Kp * t₀

    Kc = Kp * (Pa * t₀)^∑ν
    ∂Kc∂Kp = (Pa * t₀)^∑ν
    ∂Kc∂T = -Kp * ∑ν * (Pa * t₀)^∑ν / T

    dKcdT = ∂Kc∂T + (∂Kp∂T + ∂Kp∂∆H * d∆HdT + ∂Kp∂∆S * d∆SdT) * ∂Kc∂Kp
    return Kc, dKcdT
end

function reverse_rate(reaction::AbstractReaction{N}, (; T)::State{N}) where {N<:Number}
    (; rates, isreversible, reverse_rate_parameters) = reaction
    isreversible || return nothing
    isnothing(reverse_rate_parameters) || (@inbounds rates.kr.val[] = reverse_rate_parameters(T);
        return nothing
    )

    Kc = equilibrium_constants(reaction, T)
    @inbounds rates.kr.val[] = rates.kf.val[] / Kc
    return nothing
end

function reverse_rate(v::Val{:dT}, reaction::AbstractReaction{N}, (; T)::State{N}) where {N<:Number}
    (; rates, isreversible, reverse_rate_parameters) = reaction
    (; kf, kr) = rates

    isreversible || return nothing
    isnothing(reverse_rate_parameters) || ((kr.val[], kr.dT[]) = reverse_rate_parameters(v, T);
        return nothing
    )

    Kc, dKcdT = equilibrium_constants(v, reaction, T)
    t = inv(Kc)

    ∂kr∂kf = t
    @inbounds ∂kr∂Kc = -kf.val[] * t^2
    @inbounds kr.val[] = kf.val[] * t
    @inbounds kr.dT[] = ∂kr∂kf * kf.dT[] + ∂kr∂Kc * dKcdT
    return nothing
end

reverse_rate(::Val{:dC}, reaction::ElementaryReaction{N}, state::State{N}) where {N<:Number} = reverse_rate(reaction, state)
reverse_rate(::Val{:dC}, reaction::ThreeBodyReaction{N}, state::State{N}) where {N<:Number} = reverse_rate(reaction, state)

function reverse_rate(::Val{:dP}, reaction::AbstractReaction{N}, (; T)::State{N}) where {N<:Number}
    (; rates, isreversible, reverse_rate_parameters) = reaction
    isreversible || return nothing
    isnothing(reverse_rate_parameters) || (@inbounds rates.kr.val[] = reverse_rate_parameters(T);
        return nothing
    )
    
    Kc = equilibrium_constants(reaction, T)
    @inbounds rates.kr.dP[] = rates.kf.dP[] / Kc
    return nothing
end

function reverse_rate(::Val{:dC}, reaction::FallOffReaction{N}, (; T)::State{N}) where {N<:Number}
    (; rates, isreversible, reverse_rate_parameters) = reaction
    (; kf, kr) = rates

    isreversible || return nothing
    isnothing(reverse_rate_parameters) || (@inbounds kr.val[] = reverse_rate_parameters(T);
        return nothing
    )
    
    Kc = equilibrium_constants(reaction, T)
    t = inv(Kc)

    dkrdkf = t
    @inbounds kr.val[] = kf.val[] * t
    map!(kr.dC, kf.dC) do dkfdCₖ
        dkrdCₖ = dkrdkf * dkfdCₖ
        return dkrdCₖ
    end
    return nothing
end

step(part::Vector{Pair{Species{N},N}}, C::Vector{N}) where {N<:Number} = @inbounds prod(C[k]^abs(ν) for ((; k), ν) in part) ## C in small values like in units of mol/cm^3 may lead to slightly different results due to this operation
step(part::Vector{Pair{Species{N},N}}, C::Vector{N}, j::Int) where {N<:Number} = @inbounds in(j, k for ((; k), ν) in part) ? prod(k ≠ j ? C[k]^abs(ν) : abs(ν) * C[k]^(abs(ν) - 1) for ((; k), ν) in part) : zero(N)

function progress_rate(reaction::AbstractReaction{N}, (; C)::State{N}) where {N<:Number}
    (; kf, kr, q) = reaction.rates
    ∏ᴵ = step(reaction.reactants, C)
    ∏ᴵᴵ = reaction.isreversible ? step(reaction.products, C) : zero(N)

    M = reaction isa ThreeBodyReaction ? total_molar_concentration(C, reaction.enhancement_factors) : one(N)
    @inbounds q.val[] = M * (kf.val[] * ∏ᴵ - kr.val[] * ∏ᴵᴵ)
    return nothing
end

function progress_rate(::Val{:dT}, reaction::AbstractReaction{N}, (; C)::State{N}) where {N<:Number}
    (; kf, kr, q) = reaction.rates
    ∏ᴵ = step(reaction.reactants, C)
    ∏ᴵᴵ = reaction.isreversible ? step(reaction.products, C) : zero(N)

    M = reaction isa ThreeBodyReaction ? total_molar_concentration(C, reaction.enhancement_factors) : one(N)
    @inbounds q.val[] = M * (kf.val[] * ∏ᴵ - kr.val[] * ∏ᴵᴵ)
    @inbounds q.dT[] = M * (kf.dT[] * ∏ᴵ - kr.dT[] * ∏ᴵᴵ)
    return nothing
end

function progress_rate(::Val{:dC}, reaction::ElementaryReaction{N}, (; C)::State{N}) where {N<:Number}
    (; kf, kr, q) = reaction.rates
    ∏ᴵ = step(reaction.reactants, C)
    ∏ᴵᴵ = reaction.isreversible ? step(reaction.products, C) : zero(N)
    for k in eachindex(q.dC)
        d∏ᴵdCₖ = step(reaction.reactants, C, k)
        d∏ᴵᴵdCₖ = step(reaction.products, C, k)
        @inbounds q.dC[k] = (kf.val[] * d∏ᴵdCₖ) - (kr.val[] * d∏ᴵᴵdCₖ)
    end
    @inbounds q.val[] = kf.val[] * ∏ᴵ - kr.val[] * ∏ᴵᴵ
    return reaction
end

function progress_rate(::Val{:dC}, reaction::ThreeBodyReaction{N}, (; C)::State{N}) where {N<:Number}
    (; kf, kr, q) = reaction.rates
    M = total_molar_concentration(C, reaction.enhancement_factors)
    ∏ᴵ = step(reaction.reactants, C)
    ∏ᴵᴵ = reaction.isreversible ? step(reaction.products, C) : zero(N)

    for k in eachindex(q.dC)
        d∏ᴵdCₖ = step(reaction.reactants, C, k)
        d∏ᴵᴵdCₖ = step(reaction.products, C, k)
        @inbounds q.dC[k] = M * (kf.val[] * d∏ᴵdCₖ - kr.val[] * d∏ᴵᴵdCₖ)
    end

    for ((; k), dMdCₖ) in reaction.enhancement_factors
        @inbounds q.dC[k] += dMdCₖ * (kf.val[] * ∏ᴵ - kr.val[] * ∏ᴵᴵ)
    end
    @inbounds q.val[] = M * (kf.val[] * ∏ᴵ - kr.val[] * ∏ᴵᴵ)
    return reaction
end

function progress_rate(::Val{:dC}, reaction::FallOffReaction{N}, (; C)::State{N}) where {N<:Number}
    (; kf, kr, q) = reaction.rates
    ∏ᴵ = step(reaction.reactants, C)
    ∏ᴵᴵ = reaction.isreversible ? step(reaction.products, C) : zero(N)

    for k in eachindex(q.dC)
        d∏ᴵdCₖ = step(reaction.reactants, C, k)
        d∏ᴵᴵdCₖ = step(reaction.products, C, k) ## get it out of the loop somehow in case of not reversible
        @inbounds q.dC[k] = (kf.val[] * d∏ᴵdCₖ - kr.val[] * d∏ᴵᴵdCₖ) + (kf.dC[k] * ∏ᴵ - kr.dC[k] * ∏ᴵᴵ)
    end
    @inbounds q.val[] = kf.val[] * ∏ᴵ - kr.val[] * ∏ᴵᴵ
    return reaction
end

function progress_rate(::Val{:dP}, reaction::AbstractReaction{N}, (; C)::State{N}) where {N<:Number}
    (; kf, kr, q) = reaction.rates
    ∏ᴵ = step(reaction.reactants, C)
    ∏ᴵᴵ = reaction.isreversible ? step(reaction.products, C) : zero(N)

    M = reaction isa ThreeBodyReaction ? total_molar_concentration(C, reaction.enhancement_factors) : one(N)
    @inbounds q.dP[] = M * (kf.dP[] * ∏ᴵ - kr.dP[] * ∏ᴵᴵ)
    return nothing
end

_update_reaction_rates(reaction::AbstractReaction{N}, state::State{N}) where {N<:Number} = (forward_rate, reverse_rate, progress_rate)(reaction, state)
_update_reaction_rates(v::Val, reaction::AbstractReaction{N}, state::State{N}) where {N<:Number} = (forward_rate, reverse_rate, progress_rate)(v, reaction, state)