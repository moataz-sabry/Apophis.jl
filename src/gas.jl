const Mechanism{N<:Number, S<:AbstractSpecies{N}, R<:AbstractReaction{N}} = NamedTuple{
    tuple(:species, :reactions, :stoichiometry_matrix),
    Tuple{Vector{S}, Vector{R}, SparseMatrixCSC{N, Int}}
}

mutable struct State{N<:Number}
    T::N
    P::N
    ρ::N
    Y::Vector{N}
    X::Vector{N}
    C::Vector{N}
end

State(mechanism::Mechanism{N}) where {N<:Number} = State(N(Tᵣ), N(Pa), zero(N), (zeros(N, length(mechanism.species)) for _ in OneTo(3))...)

XW2W̅(X::Vector{N}, W::AbstractVector{N}) where {N<:Number} = sum(x * w for (x, w) in zip(X, W))
YW2W̅(Y::Vector{N}, W::AbstractVector{N}) where {N<:Number} = sum(y / w for (y, w) in zip(Y, W)) |> inv
CW2W̅(C::Vector{N}, W::AbstractVector{N}) where {N<:Number} = sum(c * w for (c, w) in zip(C, W)) / sum(C)

W̅YW2Y!(Y::Vector{N}, X::Vector{N}, W::AbstractVector{N}, W̅::N) where {N<:Number} = map!((x, w) -> x * w / W̅, Y, X, W)
W̅YW2X!(X::Vector{N}, Y::Vector{N}, W::AbstractVector{N}, W̅::N) where {N<:Number} = map!((y, w) -> y * W̅ / w, X, Y, W)

function C2X!(X::Vector{N}, C::Vector{N}) where {N<:Number}
    ∑C = sum(C)
    map!(c -> c / ∑C, X, C)
    return nothing
end

function CW2Y!(Y::Vector{N}, C::Vector{N}, W::AbstractVector{N}) where {N<:Number}
    ∑CW = sum(c * w for (c, w) in zip(C, W))
    map!((c, w) -> c * w / ∑CW, Y, C, W)
    return nothing
end

Cρ2Y!(Y::Vector{N}, C::Vector{N}, W::AbstractVector{N}, ρ::N) where {N<:Number} = map!((c, w) -> c * w / ρ, Y, C, W)
YWρ2C!(C::Vector{N}, Y::Vector{N}, W::AbstractVector{N}, ρ::N) where {N<:Number} = map!((y, w) -> y * ρ / w, C, Y, W)
XW̅ρ2C!(C::Vector{N}, X::Vector{N}, W̅::N, ρ::N) where {N<:Number} = map!(x -> x * ρ / W̅, C, X)

CT2P(C::Vector{N}, T::N) where {N<:Number} = sum(C) * R * T
ρW̅T2P(ρ::N, W̅::N, T::N) where {N<:Number} = ρ * R * T / W̅

CW2ρ(C::Vector{N}, W::AbstractVector{N}) where {N<:Number} = sum(c * w for (c, w) in zip(C, W))
PW̅T2ρ(P::N, W̅::N, T::N) where {N<:Number} = P * W̅ / (R * T)

struct Gas{N<:Number, S<:AbstractSpecies{N}, R<:AbstractReaction{N}}
    mechanism::Mechanism{N, S, R}
    state::State{N}
end

"""
Main function to create a `Gas` object.
"""
function Gas(name::Union{Nothing, Symbol, String} = nothing;
    kinetics_path::Maybe{String} = nothing,
    thermo_path::Maybe{String} = nothing,
    transport_path::Maybe{String} = nothing,
    as::Type{<:Number} = Float64,
    init...)

    mechanism = read_mechanism(name; kinetics_path, thermo_path, transport_path, as)
    gas = Gas(mechanism, State(mechanism))
    set!(gas; init...)
    return gas
end

Gas(mechanism::Mechanism{N}) where {N<:Number} = Gas(mechanism, State(mechanism))

function TPY!(gas::Gas{N}, T::N, P::N, Y::AbstractVector{N}) where {N<:Number}
    state = gas.state
    W = molecular_weights(gas)
    state.T = T
    state.P = P
    map!(identity, state.Y, Y)
    W̅ = YW2W̅(state.Y, W)
    W̅YW2X!(state.X, state.Y, W, W̅)
    state.ρ = PW̅T2ρ(P, W̅, T)
    XW̅ρ2C!(state.C, state.X, W̅, state.ρ) 
    return gas
end

function TPX!(gas::Gas{N}, T::N, P::N, X::AbstractVector{N}) where {N<:Number}
    state = gas.state
    W = molecular_weights(gas)
    state.T = T
    state.P = P
    map!(identity, state.X, X)
    W̅ = XW2W̅(state.X, W)
    W̅YW2Y!(state.Y, state.X, W, W̅)
    state.ρ = PW̅T2ρ(P, W̅, T)
    XW̅ρ2C!(state.C, X, W̅, state.ρ) 
    return gas
end

function TPC!(gas::Gas{N}, T::N, P::N, C::AbstractVector{N}) where {N<:Number}
    state = gas.state
    W = molecular_weights(gas)
    state.T = T
    state.P = P
    map!(identity, state.C, C)
    state.ρ = CW2ρ(state.C, W)
    C2X!(state.X, state.C)
    Cρ2Y!(state.Y, state.C, W, state.ρ)
    return gas
end

function TρY!(gas::Gas{N}, T::N, ρ::N, Y::AbstractVector{N}) where {N<:Number}
    state = gas.state
    W = molecular_weights(gas)
    state.T = T
    state.ρ = ρ
    map!(identity, state.Y, Y)
    W̅ = YW2W̅(state.Y, W)
    W̅YW2X!(state.X, state.Y, W, W̅)
    state.P = ρW̅T2P(ρ, W̅, T)
    XW̅ρ2C!(state.C, state.X, W̅, ρ) 
    return gas
end

function TρX!(gas::Gas{N}, T::N, ρ::N, X::AbstractVector{N}) where {N<:Number}
    state = gas.state
    W = molecular_weights(gas)
    state.T = T
    state.ρ = ρ
    map!(identity, state.X, X)
    W̅ = XW2W̅(state.X, W)
    W̅YW2Y!(state.Y, state.X, W, W̅)
    state.P = ρW̅T2P(ρ, W̅, T)
    XW̅ρ2C!(state.C, X, W̅, state.ρ) 
    return gas
end

function TρC!(gas::Gas{N}, T::N, ρ::N, C::AbstractVector{N}) where {N<:Number}
    state = gas.state
    W = molecular_weights(gas)
    state.T = T
    state.ρ = ρ
    map!(identity, state.C, C)
    state.P = CT2P(state.C, state.T)
    C2X!(state.X, state.C)
    Cρ2Y!(state.Y, state.C, W, state.ρ) 
    return gas
end