function _kinetics_sensitivity(reaction::Reaction{N}, (; T, C)::State{N}) where {N<:Number}
    ∏ᴵ = step(reaction.reactants, C)
    M = reaction isa ThreeBodyReaction ? total_molar_concentration(C, reaction.enhancement_factors) : one(N)

    if reaction.isreversible
        Kc = equilibrium_constants(reaction, T)
        ∏ᴵᴵ = step(reaction.products, C)
        dkrdkf = inv(Kc)
    else
        Kc = ∏ᴵᴵ = dkrdkf = zero(N)
    end
    
    dqdkf = M * (∏ᴵ - dkrdkf * ∏ᴵᴵ)
    dqdkr = -M * ∏ᴵᴵ
    return dqdkf, dqdkr
end

_kinetics_sensitivity((; mechanism, state)::Gas{<:Number}, d::Int = 1) = @inbounds map(reaction -> _kinetics_sensitivity(reaction, state)[d], mechanism.reactions) |> Diagonal
kinetics_sensitivity(gas::Gas{N}, d::Int = 1) where {N<:Number} = stoichiometry_matrix(gas) * _kinetics_sensitivity(gas, d)

_dkfdA((; forward_rate_parameters)::Union{ElementaryReaction{N}, ThreeBodyReaction{N}}, (; T)::State{N}) where {N<:Number} = forward_rate_parameters(Val(:dg), T) |> first

function _dkfdA((; high_pressure_parameters, low_pressure_parameters, enhancement_factors, troe_parameters)::FallOffReaction{N}, (; T, C)::State{N}) where {N<:Real} ## Add dkfdAₒ later
    k∞, kₒ = (high_pressure_parameters, low_pressure_parameters)(T)
    M = total_molar_concentration(C, enhancement_factors)
    Pᵣ = reduced_pressure(kₒ, M, k∞)
    t = inv(one(N) + Pᵣ)
    
    dk∞dA∞ = high_pressure_parameters(Val(:dg), T) |> first
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

dkfdA((; mechanism, state)::Gas{<:Number}) = map(reaction -> _dkfdA(reaction, state), mechanism.reactions) |> Diagonal
dω̇dA(gas::Gas{<:Number}) = stoichiometry_matrix(gas) * dkfdA(gas)