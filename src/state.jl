mutable struct State{N<:Number}
    T::N
    P::N
    ρ::N
    Y::Vector{N}
    X::Vector{N}
    C::Vector{N}
end

State(mechanism::NamedTuple{(:species, :reactions), Tuple{Vector{Species{N}}, Vector{R}}}) where {N<:Number, R} = State(N(Tᵣ), N(Pa), zero(N), (zeros(N, length(mechanism.species)) for _ in OneTo(3))...)

XW2W̅(X::Vector{N}, species::Vector{Species{N}}) where {N<:Number} = sum(x * s.weight for (x, s) in zip(X, species))
YW2W̅(Y::Vector{N}, species::Vector{Species{N}}) where {N<:Number} = sum(y * inv(s.weight) for (y, s) in zip(Y, species)) |> inv
CW2W̅(C::Vector{N}, species::Vector{Species{N}}) where {N<:Number} = sum(c * s.weight for (c, s) in zip(C, species)) / sum(C)

W̅YW2Y!(Y::Vector{N}, X::Vector{N}, species::Vector{Species{N}}, W̅::N) where {N<:Number} = map!((x, s) -> x * s.weight / W̅, Y, X, species)
W̅YW2X!(X::Vector{N}, Y::Vector{N}, species::Vector{Species{N}}, W̅::N) where {N<:Number} = map!((y, s) -> y * W̅ * inv(s.weight), X, Y, species)

function C2X!(X::Vector{N}, C::Vector{N}) where {N<:Number}
    ∑C = sum(C)
    map!(c -> c / ∑C, X, C)
    return nothing
end

function CW2Y!(Y::Vector{N}, C::Vector{N}, species::Vector{Species{N}}) where {N<:Number}
    ∑CW = sum(c * s.weight for (c, s) in zip(C, species))
    map!((c, s) -> c * s.weight, Y, C, species) / ∑CW
    return nothing
end

Cρ2Y!(Y::Vector{N}, C::Vector{N}, species::Vector{Species{N}}, ρ::N) where {N<:Number} = map!((c, s) -> c * s.weight, Y, C, species) / ρ
YWρ2C!(C::Vector{N}, Y::Vector{N}, species::Vector{Species{N}}, ρ::N) where {N<:Number} = map!((y, s) -> y * ρ * inv(s.weight), C, Y, species)
XW̅ρ2C!(C::Vector{N}, X::Vector{N}, W̅::N, ρ::N) where {N<:Number} = map!(x -> x * ρ / W̅, C, X)

CT2P(C::Vector{N}, T::N) where {N<:Number} = sum(C) * R * T
ρW̅T2P(ρ::N, W̅::N, T::N) where {N<:Number} = ρ * R * T / W̅

CW2ρ(C::Vector{N}, species::Vector{Species{N}}) where {N<:Number} = sum(c * s.weight for (c, s) in zip(C, species))
PW̅T2ρ(P::N, W̅::N, T::N) where {N<:Number} = P * W̅ / (R * T)

struct Gas{N<:Number, R<:AbstractReaction{N}}
    mechanism::NamedTuple{(:species, :reactions), Tuple{Vector{Species{N}}, Vector{R}}}
    state::State{N}
end

"""
Main function to create a `Gas` object.
"""
function Gas(name::Union{Nothing, Symbol,String} = nothing; 
    kinetics_path = "", 
    thermo_path = "", 
    transport_path = "",
    as::Type{<:Number} = Float64, init...)

    mechanism = read_mechanism(name; kinetics_path = kinetics_path, thermo_path = thermo_path, transport_path = transport_path, as = as)
    gas = Gas(mechanism, State(mechanism))

    set!(gas; init...)
    return gas
end

Gas(mechanism::NamedTuple{(:species, :reactions), Tuple{Vector{Species{N}}, Vector{R}}}) where {N<:Number, R<:AbstractReaction{N}} = Gas{N, R}(mechanism, State(mechanism))

function TPY!(gas::Gas{N}, T::N, P::N, Y::AbstractVector{N}) where {N<:Number}
    (; mechanism, state) = gas
    state.T = T
    state.P = P
    map!(identity, state.Y, Y)
    W̅ = YW2W̅(state.Y, mechanism.species)
    W̅YW2X!(state.X, state.Y, mechanism.species, W̅)
    state.ρ = PW̅T2ρ(P, W̅, T)
    XW̅ρ2C!(state.C, state.X, W̅, state.ρ) 
    return gas
end

function TPX!(gas::Gas{N}, T::N, P::N, X::AbstractVector{N}) where {N<:Number}
    (; mechanism, state) = gas
    state.T = T
    state.P = P
    map!(identity, state.X, X)
    W̅ = XW2W̅(state.X, mechanism.species)
    W̅YW2Y!(state.Y, state.X, mechanism.species, W̅)
    state.ρ = PW̅T2ρ(P, W̅, T)
    XW̅ρ2C!(state.C, X, W̅, state.ρ) 
    return gas
end

function TPC!(gas::Gas{N}, T::N, P::N, C::AbstractVector{N}) where {N<:Number}
    (; mechanism, state) = gas
    state.T = T
    state.P = P
    map!(identity, state.C, C)
    state.ρ = CW2ρ(state.C, mechanism.species)
    C2X!(state.X, state.C)
    Cρ2Y!(state.Y, state.C, mechanism.species, state.ρ)
    return gas
end

function TρY!(gas::Gas{N}, T::N, ρ::N, Y::AbstractVector{N}) where {N<:Number}
    (; mechanism, state) = gas
    state.T = T
    state.ρ = ρ
    map!(identity, state.Y, Y)
    W̅ = YW2W̅(state.Y, mechanism.species)
    W̅YW2X!(state.X, state.Y, mechanism.species, W̅)
    state.P = ρW̅T2P(ρ, W̅, T)
    XW̅ρ2C!(state.C, state.X, W̅, ρ) 
    return gas
end

function TρX!(gas::Gas{N}, T::N, ρ::N, X::AbstractVector{N}) where {N<:Number}
    (; mechanism, state) = gas
    state.T = T
    state.ρ = ρ
    map!(identity, state.X, X)
    W̅ = XW2W̅(state.X, mechanism.species)
    W̅YW2Y!(state.Y, state.X, mechanism.species, W̅)
    state.P = ρW̅T2P(ρ, W̅, T)
    XW̅ρ2C!(state.C, X, W̅, state.ρ) 
    return gas
end

function TρC!(gas::Gas{N}, T::N, ρ::N, C::AbstractVector{N}) where {N<:Number}
    (; mechanism, state) = gas
    state.T = T
    state.ρ = ρ
    map!(identity, state.C, C)
    state.P = CT2P(state.C, state.T)
    C2X!(state.X, state.C)
    Cρ2Y!(state.Y, state.C, mechanism.species, state.ρ) 
    return gas
end