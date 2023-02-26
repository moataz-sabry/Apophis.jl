function _kinetics_sensitivity(reaction::AbstractReaction{N}, (; T, C)::State{N}) where {N<:Number}
    ∏ᴵ = step(reaction.reactants, C)
    M = reaction isa ThreeBodyReaction ? total_molar_concentration(C, reaction.enhancement_factors) : one(N)

    if reaction.isreversible && isnothing(reaction.reverse_rate_parameters)
        Kc = equilibrium_constants(reaction, T)
        ∏ᴵᴵ = step(reaction.products, C)
        dkrdkf = inv(Kc)
    else
        ∏ᴵᴵ = dkrdkf = zero(N)
    end

    dqdkf = M * (∏ᴵ - dkrdkf * ∏ᴵᴵ)
    dqdkr = -M * ∏ᴵᴵ
    return dqdkf, dqdkr
end

_dkfdA((; forward_rate_parameters)::AbstractReaction{N}, (; T)::State{N}) where {N<:Number} = forward_rate_parameters(Val(:dg), T)[:A]

function _dkfdA((; high_pressure_parameters, low_pressure_parameters, enhancement_factors, troe_parameters)::FallOffReaction{N}, (; T, C)::State{N}) where {N<:Number} ## Add dkfdAₒ later
    k∞, kₒ = (high_pressure_parameters, low_pressure_parameters)(T)
    M = total_molar_concentration(C, enhancement_factors)
    iszero(M) && return zero(N)

    Pᵣ = reduced_pressure(kₒ, M, k∞)
    t = inv(one(N) + Pᵣ)

    dk∞dA∞ = high_pressure_parameters(Val(:dg), T)[:A]
    dPᵣdk∞ = reduced_pressure(Val(:dk∞), kₒ, M, k∞)

    if isnothing(troe_parameters)
        dkfdk∞ = Pᵣ * t
        dkfdPᵣ = k∞ * t^2
        dkfdA∞ = dkfdk∞ * dk∞dA∞ + dkfdPᵣ * dPᵣdk∞ * dk∞dA∞
    else
        Fc = troe_parameters(T)
        F, dFdPᵣ = troe_function(Val(:dC), Fc, Pᵣ)

        dkfdF = k∞ * Pᵣ * t
        dkfdk∞ = F * Pᵣ * t
        dkfdPᵣ = F * k∞ * t^2

        dFdA = dFdPᵣ * dPᵣdk∞ * dk∞dA∞
        dkfdA∞ = dkfdF * dFdA + dkfdk∞ * dk∞dA∞ + dkfdPᵣ * dPᵣdk∞ * dk∞dA∞
    end
    return dkfdA∞
end

function dω̇dA((; mechanism, state)::Gas{N}) where {N<:Number}
    dω̇dA = zeros(N, length(mechanism.species), length(mechanism.reactions))
    for (i, reaction) in enumerate(mechanism.reactions)
        dqdA = _kinetics_sensitivity(reaction, state)[1] * _dkfdA(reaction, state)
        for ((; k), ν) in flatten((reaction.reactants, reaction.products))
            @inbounds dω̇dA[k, i] = ν * dqdA
        end
    end
    return dω̇dA
end